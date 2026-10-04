"""Count-level WGBS simulation with immutable CpG and DMR truth."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.special import expit, logit, betaincinv, ndtr

from wgbs_v3.provenance import digest

SIMULATOR_VERSION = 4


@dataclass(frozen=True)
class Scenario:
    seed: int
    n_per_group: int = 5
    coverage: int = 20
    delta: float = 0.2
    dispersion_multiplier: float = 1.0
    regions: int | None = None
    correlated: bool = False
    mode: str = "dmr"  # dmr, dml, null
    dml_count: int = 1000
    balanced_coverage: bool = False

    def check(self):
        if self.n_per_group < 2 or self.coverage not in (5, 10, 15, 20, 30):
            raise ValueError("Samples must be >=2; coverage must be one of 5,10,15,20,30")
        if not 0 <= self.delta < 0.5 or self.dispersion_multiplier <= 0:
            raise ValueError("Invalid effect or dispersion")
        if self.mode not in ("dmr", "dml", "null"):
            raise ValueError("Unknown simulation mode")
        if self.mode != "null" and self.delta == 0:
            raise ValueError("Signal scenarios require a positive effect")
        if self.regions is not None and self.regions < 1:
            raise ValueError("regions must be positive")
        if self.dml_count < 1:
            raise ValueError("dml_count must be positive")


def _rng(seed: int, chromosome: str, stream: str):
    key = int.from_bytes(hashlib.sha256(f"{seed}/{chromosome}/{stream}".encode()).digest()[:8], "little")
    return np.random.default_rng(key)


def _region_shift(values: np.ndarray, delta: float, direction: int) -> np.ndarray:
    target = np.clip(values.mean() + direction * delta, 0.01, 0.99)
    if abs(target - values.mean()) < delta - 1e-6:
        raise ValueError("Target effect exceeds methylation range")
    base = logit(np.clip(values, 0.01, 0.99))
    low, high = -40.0, 40.0
    for _ in range(60):
        mid = (low + high) / 2
        if expit(base + mid).mean() < target:
            low = mid
        else:
            high = mid
    return expit(base + (low + high) / 2)


def place_regions(pos: np.ndarray, mu: np.ndarray, target: int, seed: int, chrom: str):
    """Place nonoverlapping truth regions on successive, sufficiently dense CpGs."""
    rng = _rng(seed, chrom, "placement")
    occupied = np.zeros(len(pos), dtype=bool)
    regions = []
    tries = 0
    while len(regions) < target and tries < max(1000, target * 500):
        tries += 1
        n = max(5, int(round(rng.gamma(4, 5))))
        if n + 2 >= len(pos):
            break
        start = int(rng.integers(1, len(pos) - n - 1))
        end = start + n
        gaps = np.diff(pos[start:end])
        if occupied[start - 1:end + 1].any() or gaps.max(initial=0) > 500:
            continue
        if (pos[end - 1] - pos[start] + 1) / n > 100:
            continue
        direction = 1 if mu[start:end].mean() <= 0.5 else -1
        occupied[start:end] = True
        regions.append((start, end, direction))
    if len(regions) < target:
        raise ValueError(f"Only {len(regions)}/{target} eligible regions on {chrom}")
    return sorted(regions)


def _chrom_catalog(files: list[str], chrom: str):
    paths = [f for f in files if Path(f).name.startswith(chrom + "_")]
    # Avoid chr1 matching chr10; the underscore above is significant.
    if not paths:
        return None
    table = pa.concat_tables([pq.read_table(f) for f in sorted(paths)])
    return table.to_pandas()


def _counts(mu: np.ndarray, tau: np.ndarray, coverage: np.ndarray, rng,
            correlated: bool, pos, length_bp: float | None, count_rng=None):
    alpha = mu / tau
    beta = (1 - mu) / tau
    if correlated:
        if length_bp is None:
            raise ValueError("Correlation scenario requires a fitted BLUEPRINT decay length")
        z = rng.normal(size=len(mu))
        for i in range(1, len(mu)):
            correlation = np.exp(-min(pos[i] - pos[i - 1], 10_000) / length_bp)
            z[i] = correlation * z[i - 1] + np.sqrt(1 - correlation**2) * z[i]
        probability = betaincinv(alpha, beta, np.clip(ndtr(z), 1e-7, 1 - 1e-7))
    else:
        probability = rng.beta(alpha, beta)
    return (rng if count_rng is None else count_rng).binomial(coverage, probability)


def coverage_parameters(frame, balanced=False):
    """Expected depth after dropout has mean 30, before sample scaling."""
    depth = frame.coverage_mean.to_numpy(dtype=float)
    if np.any(~np.isfinite(depth)) or np.any(depth < 0) or not np.any(depth > 0):
        raise ValueError("Invalid empirical depth library")
    dropout = np.clip(frame.zero_fraction.to_numpy(dtype=float) - np.exp(-depth), 0, .95)
    if balanced:
        dropout = np.zeros(len(depth))
    normalizer = np.mean(depth * (1 - dropout))
    return 30 * depth / normalizer, dropout


def _thin(meth, cov, target, rng):
    # Each next tier is a nested Bernoulli subsample of the preceding tier.
    p = target / 30 if target == 20 else {15: 15 / 20, 10: 10 / 15, 5: 5 / 10}[target]
    m = rng.binomial(meth, p)
    u = rng.binomial(cov - meth, p)
    return m, m + u


def simulate(calibration: Path, out: Path, config: Scenario, *, chromosomes: list[str] | None = None,
             max_sites: int | None = None) -> dict:
    config.check()
    manifest = json.loads((calibration / "manifest.json").read_text())
    if manifest["release_ready"]:
        proof = calibration / "reference_validation.json"
        validation = json.loads(proof.read_text()) if proof.exists() else {}
        if not validation.get("passed") or validation.get("calibration_sha256") != digest(calibration / "manifest.json"):
            raise ValueError("Release simulation requires a passed GRCh38 CpG validation")
    out.mkdir(parents=True, exist_ok=False)
    chroms = chromosomes or manifest.get("scope", {}).get("chromosomes") or [f"chr{i}" for i in range(1, 23)]
    files = manifest["parquet_files"]
    per_chrom = {}
    for chrom in chroms:
        frame = _chrom_catalog(files, chrom)
        if frame is None:
            raise ValueError(f"Calibration missing requested chromosome {chrom}; no silent subset simulation")
        if frame is not None:
            if max_sites is not None:
                frame = frame.iloc[:max_sites]
            per_chrom[chrom] = frame
    total = sum(len(x) for x in per_chrom.values())
    if total == 0:
        raise ValueError("No calibrated CpGs for selected chromosomes")
    if config.mode == "dml" and total > 100_000:
        selected_global = np.sort(_rng(config.seed, "all", "dml_panel").choice(total, 100_000, replace=False))
        reduced = {}
        offset = 0
        for chrom, frame in per_chrom.items():
            indices = selected_global[(selected_global >= offset) & (selected_global < offset + len(frame))] - offset
            if len(indices):
                reduced[chrom] = frame.iloc[indices]
            offset += len(frame)
        per_chrom = reduced
        total = 100_000
    if config.mode == "dml" and total < config.dml_count:
        raise ValueError("DML panel has fewer loci than positives")
    if config.mode == "dmr":
        region_total = config.regions if config.regions is not None else max(1, round(total / 26_883_210 * 3000))
        fractions = np.array([len(x) / total * region_total for x in per_chrom.values()])
        quotas = np.floor(fractions).astype(int)
        for idx in np.argsort(-(fractions - quotas))[:region_total - quotas.sum()]:
            quotas[idx] += 1
        region_counts = dict(zip(per_chrom, quotas))
    else:
        region_counts = {c: 0 for c in per_chrom}
    dml_global = set(_rng(config.seed, "all", "dml_placement").choice(total, config.dml_count, replace=False).tolist()) if config.mode == "dml" else set()
    tau_library = np.asarray(manifest["tau_by_mean_bin"])
    spatial_model = manifest["spatial_correlation"]
    spatial_length = spatial_model["fitted_length_bp"]
    if config.correlated and not spatial_model.get("validated", False):
        raise ValueError("Spatial correlation scenario is unavailable until decay validation passes")
    depth_library = np.asarray(manifest["sample_mean_depth"], dtype=float)
    depth_library /= depth_library.mean()
    truth_dir = out / "truth_parts"
    truth_dir.mkdir()
    eligibility_dir = out / "eligibility_parts"
    eligibility_dir.mkdir()
    count_dir = out / "counts"
    count_dir.mkdir()
    sample_names = [f"control_{i+1}" for i in range(config.n_per_group)] + [f"treatment_{i+1}" for i in range(config.n_per_group)]
    # One sample-level draw persists across chromosomes. RNG keys are independent
    # of effect placement, biological sampling, coverage and thinning.
    library_factors = _rng(config.seed, "all", "library").choice(depth_library, size=len(sample_names), replace=True)
    library_factors = np.ones(len(sample_names)) if config.balanced_coverage else library_factors / library_factors.mean()
    streams = {name: io.TextIOWrapper(gzip.GzipFile(filename="", mode="wb", compresslevel=3,
               fileobj=(count_dir / f"{name}.cov.gz").open("wb"), mtime=0)) for name in sample_names}
    (out / "samples.tsv").write_text("sample_id\tgroup\tpath\n" + "".join(
        f"{name}\t{'control' if name.startswith('control') else 'treatment'}\t{(count_dir / (name + '.cov.gz')).resolve()}\n"
        for name in sample_names))
    regions = []
    covered = 0
    positives = 0
    dml_offset = 0
    try:
        for chrom, frame in per_chrom.items():
            print(f"Simulating {config.mode} {chrom}: {len(frame):,} CpGs", flush=True)
            pos = frame.pos.to_numpy(dtype=np.int64)
            mu = frame.mu.to_numpy(dtype=float)
            bin_idx = np.clip(np.floor(mu * 100).astype(int) - 1, 0, 98)
            tau = np.clip(tau_library[bin_idx] * config.dispersion_multiplier, 1e-4, 4)
            treatment = mu.copy()
            ids = np.zeros(len(pos), dtype=np.int32)
            if config.mode == "dmr":
                placed = place_regions(pos, mu, region_counts[chrom], config.seed, chrom)
                for j, (lo, hi, direction) in enumerate(placed, start=1):
                    key = len(regions) + 1
                    treatment[lo:hi] = _region_shift(mu[lo:hi], config.delta, direction)
                    ids[lo:hi] = key
                    regions.append({"region_id": key, "chrom": chrom, "start": int(pos[lo]),
                                    "end": int(pos[hi-1]), "n_cpg": hi - lo, "direction": direction,
                                    "mean_delta": float(np.mean(treatment[lo:hi] - mu[lo:hi])),
                                    "min_abs_delta": float(np.min(np.abs(treatment[lo:hi] - mu[lo:hi]))),
                                    "max_abs_delta": float(np.max(np.abs(treatment[lo:hi] - mu[lo:hi])))})
            elif config.mode == "dml":
                selected = np.array([i for i in range(len(pos)) if dml_offset + i in dml_global], dtype=int)
                direction = np.where(mu[selected] <= 0.5, 1, -1)
                treatment[selected] = mu[selected] + direction * config.delta
                ids[selected] = -1
            dml_offset += len(pos)
            effect = treatment - mu
            if np.any((ids == 0) & (effect != 0)) or np.any((ids != 0) & (effect == 0)):
                raise AssertionError("Truth and latent effect disagree")
            truth = pa.table({"chrom": [chrom] * len(pos), "pos": pos,
                              "control_mu": mu.astype("float32"),
                              "treatment_mu": treatment.astype("float32"),
                              "signed_delta": effect.astype("float32"),
                              "tau": tau.astype("float32"), "is_dml": ids != 0,
                              "region_id": ids})
            pq.write_table(truth, truth_dir / f"{chrom}.parquet", compression="zstd")
            positives += int(np.count_nonzero(ids))
            base_depth, dropout = coverage_parameters(frame, config.balanced_coverage)
            counts = np.zeros((len(pos), len(sample_names)), dtype=np.int32)
            meth = np.zeros_like(counts)
            for col, name in enumerate(sample_names):
                coverage_rng = _rng(config.seed, chrom, "coverage_" + name)
                expected = base_depth * library_factors[col]
                counts[:, col] = coverage_rng.poisson(expected)
                counts[coverage_rng.random(len(pos)) < dropout, col] = 0
                group_mu = mu if name.startswith("control") else treatment
                meth[:, col] = _counts(group_mu, tau, counts[:, col],
                                       _rng(config.seed, chrom, "biology_" + name),
                                       config.correlated, pos, spatial_length,
                                       count_rng=_rng(config.seed, chrom, "sampling_" + name))
                if config.coverage != 30:
                    thin_rng = _rng(config.seed, chrom, "thin_" + name)
                    for tier in (20, 15, 10, 5):
                        meth[:, col], counts[:, col] = _thin(meth[:, col], counts[:, col], tier, thin_rng)
                        if config.coverage == tier:
                            break
            eligible = counts.min(axis=1) >= 1
            covered += int(eligible.sum())
            pq.write_table(pa.table({"chrom": [chrom] * len(pos), "pos": pos,
                                     "eligible": eligible}),
                           eligibility_dir / f"{chrom}.parquet", compression="zstd")
            for col, name in enumerate(sample_names):
                valid = np.flatnonzero(counts[:, col] > 0)
                stream = streams[name]
                for i in valid:
                    m, n = int(meth[i, col]), int(counts[i, col])
                    stream.write(f"{chrom}\t{pos[i]}\t{pos[i]+1}\t{100*m/n:.6f}\t{m}\t{n-m}\n")
    finally:
        for stream in streams.values():
            stream.close()
    import pandas as pd
    pd.DataFrame(regions, columns=["region_id", "chrom", "start", "end", "n_cpg", "direction", "mean_delta", "min_abs_delta", "max_abs_delta"]).to_csv(out / "dmr_truth.tsv", sep="\t", index=False)
    result = {"simulator_version": SIMULATOR_VERSION,
              "simulator_sha256": digest(Path(__file__)),
              "chromosomes": list(per_chrom), "max_sites_per_chromosome": max_sites,
              "library_factors": dict(zip(sample_names, library_factors.tolist())),
              "rng_contract": "independent placement/library/coverage/biology/sampling/thinning streams",
              "scenario": asdict(config), "calibration_sha256": digest(calibration / "manifest.json"),
              "calibration_release_ready": manifest["release_ready"], "n_sites": total,
              "n_positive_cpg": positives, "n_regions": len(regions),
              "n_common_eligible": covered, "coordinate_system": "GRCh38 1-based CpG; Bismark end=pos+1"}
    (out / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
