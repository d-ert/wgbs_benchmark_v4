"""Chunked calibration from the 206 BLUEPRINT WGBS methylomes."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from wgbs_v3.provenance import digest
from .reference import ReferenceGenome


def estimate_parameters(methylated: np.ndarray, coverage: np.ndarray, minimum: int = 20):
    """Moment fit of beta-binomial mean and rho, correcting binomial sampling noise.

    A beta-binomial has Var(M/N) = mu(1-mu)[rho + (1-rho)/N].
    Estimates are site specific. The fitted mean-bin median tau drives simulation.
    """
    ok = coverage >= 5
    count = ok.sum(axis=1)
    beta = np.divide(methylated, coverage, out=np.zeros_like(methylated, dtype=float), where=ok)
    mu = np.divide(beta.sum(axis=1), count, out=np.zeros(len(count)), where=count > 0)
    mu = np.clip(mu, 0.01, 0.99)
    sq = np.divide((((beta - mu[:, None]) * ok) ** 2).sum(axis=1), np.maximum(count - 1, 1))
    inverse_n = np.divide(np.divide(ok, coverage, out=np.zeros_like(beta), where=ok).sum(axis=1),
                          count, out=np.zeros(len(count)), where=count > 0)
    normalized = sq / (mu * (1 - mu))
    rho = np.clip((normalized - inverse_n) / np.maximum(1 - inverse_n, 1e-8), 1e-4, 0.8)
    tau = rho / (1 - rho)
    valid = (count >= minimum) & np.isfinite(tau)
    return mu, tau, count, valid


def _matrix_chunks(path: Path, chunk_rows: int):
    with gzip.open(path, "rt") as stream:
        for frame in pd.read_csv(stream, sep="\t", chunksize=chunk_rows):
            yield frame


def _read_intervals(handle, chrom: str, lo: int, hi: int, coords: np.ndarray):
    result = np.zeros(len(coords), dtype=np.float32)
    spans = handle.intervals(chrom, max(0, lo - 1), hi + 1) or ()
    if not spans:
        return result
    starts = np.fromiter((s + 1 for s, e, v in spans if e - s == 2), dtype=np.int64)
    vals = np.fromiter((v for s, e, v in spans if e - s == 2), dtype=np.float32)
    idx = np.searchsorted(coords, starts)
    good = idx < len(coords)
    good[good] &= coords[idx[good]] == starts[good]
    result[idx[good]] = vals[good]
    return result


def calibrate(matrix: Path, out: Path, samples: list[dict], source_dir: Path,
              minimum: int = 20, chunk_rows: int = 250_000, chromosomes: list[str] | None = None,
              smoke_local_matrix: bool = False, max_sites: int | None = None,
              reference_fasta: Path | None = None):
    import pyBigWig

    if not smoke_local_matrix and len(samples) != 206:
        raise ValueError("Release calibration requires all 206 Table S1 samples")
    if not smoke_local_matrix:
        verified_file = source_dir / "verified_sources.json"
        if not verified_file.is_file():
            raise ValueError("Verify all 206 source file pairs before release calibration")
        verified = json.loads(verified_file.read_text())
        if len(verified["verified"]) != 206 or {r["sample_id"] for r in verified["verified"]} != {r["sample_id"] for r in samples}:
            raise ValueError("Verified source manifest does not match Table S1")
    if reference_fasta is None or not reference_fasta.is_file():
        raise ValueError("A GRCh38 reference FASTA is required for calibration")
    reference = ReferenceGenome(reference_fasta)
    out.mkdir(parents=True, exist_ok=True)
    sample_ids = [r["sample_id"] for r in samples]
    handles = {}
    if not smoke_local_matrix:
        for sample_id in sample_ids:
            pair = []
            for kind in ("coverage", "call"):
                path = source_dir / f"{sample_id}.{kind}.bw"
                if not path.is_file():
                    raise FileNotFoundError(f"Missing required 206-sample source: {path}")
                pair.append(pyBigWig.open(str(path)))
            handles[sample_id] = pair

    hist = np.zeros((99, 400), dtype=np.int64)
    counts = {"all": 0, "valid": 0, "non_reference_cpg": 0,
              "reconstruction_disagreements": 0}
    output_files = []
    sample_depth_sum = np.zeros(len(sample_ids), dtype=np.float64)
    spatial_edges = np.array([0, 25, 50, 100, 200, 500, 1000, 2000], dtype=int)
    spatial_sum = np.zeros(len(spatial_edges) - 1, dtype=float)
    spatial_count = np.zeros(len(spatial_edges) - 1, dtype=np.int64)
    chroms = chromosomes or [f"chr{i}" for i in range(1, 23)]
    chunk_no = {}
    try:
        for frame in _matrix_chunks(matrix, chunk_rows):
            if max_sites is not None:
                remaining = max_sites - counts["all"]
                if remaining <= 0:
                    break
                frame = frame.iloc[:remaining]
            # The ten-sample matrix supplies a reference CpG coordinate catalogue.
            for chrom, sub in frame.groupby("chr", sort=False):
                if chrom not in chroms:
                    continue
                sub = sub.sort_values("pos")
                pos = sub["pos"].to_numpy(dtype=np.int64)
                if len(pos) == 0 or np.any(np.diff(pos) <= 0):
                    raise ValueError(f"Unsorted or duplicate CpGs on {chrom}")
                reference_ok = reference.is_cpg(chrom, pos)
                counts["non_reference_cpg"] += int((~reference_ok).sum())
                sub = sub.iloc[np.flatnonzero(reference_ok)]
                pos = pos[reference_ok]
                if len(pos) == 0:
                    continue
                if smoke_local_matrix:
                    cov = sub[[f"Cov{i}" for i in range(1, 11)]].to_numpy(dtype=np.int32)
                    meth = sub[[f"M{i}" for i in range(1, 11)]].to_numpy(dtype=np.int32)
                else:
                    cov = np.zeros((len(pos), len(sample_ids)), dtype=np.int32)
                    meth = np.zeros_like(cov)
                    for j, sid in enumerate(sample_ids):
                        cov_values = _read_intervals(handles[sid][0], chrom, int(pos[0]), int(pos[-1]), pos)
                        call_values = _read_intervals(handles[sid][1], chrom, int(pos[0]), int(pos[-1]), pos)
                        cov[:, j] = np.rint(cov_values).astype(np.int32)
                        meth[:, j] = np.rint(call_values * cov[:, j]).astype(np.int32)
                        bad = (cov_values > 0) & (np.abs(cov_values - np.rint(cov_values)) > 1e-4)
                        counts["reconstruction_disagreements"] += int(bad.sum())
                if np.any(meth < 0) or np.any(meth > cov):
                    raise ValueError(f"Impossible methylated counts at {chrom}")
                mu, tau, n_obs, valid = estimate_parameters(meth, cov, minimum)
                sample_depth_sum += cov.sum(axis=0)
                if len(pos) > 1:
                    observed = cov >= 5
                    fractions = np.divide(meth, cov, out=np.zeros_like(meth, dtype=float), where=observed)
                    residual = (fractions - mu[:, None]) * observed
                    site_sd = np.sqrt(np.divide((residual**2).sum(axis=1), np.maximum(observed.sum(axis=1) - 1, 1)))
                    residual /= np.maximum(site_sd[:, None], 1e-3)
                    gap = np.diff(pos)
                    bucket = np.searchsorted(spatial_edges, gap, side="right") - 1
                    pair_valid = observed[:-1] & observed[1:]
                    pair_product = np.clip(residual[:-1] * residual[1:], -16, 16) * pair_valid
                    pair_count = pair_valid.sum(axis=1)
                    pair_sum = pair_product.sum(axis=1)
                    for k in range(len(spatial_count)):
                        mask = (bucket == k) & valid[:-1] & valid[1:]
                        spatial_sum[k] += float(pair_sum[mask].sum())
                        spatial_count[k] += int(pair_count[mask].sum())
                indices = np.clip(np.floor(mu[valid] * 100).astype(int) - 1, 0, 98)
                log_tau = np.clip(np.log10(tau[valid]), -4, np.log10(4))
                tau_bins = np.minimum(np.floor((log_tau + 4) / (np.log10(4) + 4) * 400).astype(int), 399)
                np.add.at(hist, (indices, tau_bins), 1)
                counts["all"] += len(pos)
                counts["valid"] += int(valid.sum())
                chunk_no[chrom] = chunk_no.get(chrom, 0) + 1
                target = out / f"{chrom}_{chunk_no[chrom]:05d}.parquet"
                pq.write_table(pa.table({"chrom": [chrom] * int(valid.sum()),
                                         "pos": pos[valid], "mu": mu[valid].astype("float32"),
                                         "tau_fit": tau[valid].astype("float32"),
                                         "coverage_mean": cov.mean(axis=1)[valid].astype("float32"),
                                         "zero_fraction": (cov == 0).mean(axis=1)[valid].astype("float32"),
                                         "n_observed": n_obs[valid].astype("int16")}), target,
                               compression="zstd")
                output_files.append(str(target))
                print(f"Calibrated {chrom} chunk {chunk_no[chrom]}: {counts['valid']:,} eligible CpGs so far", flush=True)
            if max_sites is not None and counts["all"] >= max_sites:
                break
    finally:
        for pair in handles.values():
            for handle in pair:
                handle.close()
    np.savez(out / "sufficient_statistics.npz", hist=hist, depth_sum=sample_depth_sum,
             spatial_sum=spatial_sum, spatial_count=spatial_count)
    tau_edges = np.linspace(-4, np.log10(4), 401)
    tau_by_bin = []
    for row in hist:
        if row.sum() == 0:
            tau_by_bin.append(None)
        else:
            k = int(np.searchsorted(np.cumsum(row), (row.sum() + 1) / 2))
            tau_by_bin.append(float(10 ** ((tau_edges[k] + tau_edges[k + 1]) / 2)))
    # Empty centile bins inherit the nearest observed bin; record those bins.
    observed = np.flatnonzero([v is not None for v in tau_by_bin])
    if len(observed) == 0:
        raise ValueError("Calibration produced no eligible CpGs")
    empty = []
    for i, value in enumerate(tau_by_bin):
        if value is None:
            empty.append(i + 1)
            tau_by_bin[i] = tau_by_bin[observed[np.argmin(np.abs(observed - i))]]
    correlations = np.divide(spatial_sum, spatial_count, out=np.zeros_like(spatial_sum), where=spatial_count > 0)
    mid = np.sqrt(np.maximum(spatial_edges[:-1], 1) * spatial_edges[1:])
    use = (spatial_count >= 1000) & (correlations > .01) & (correlations < .99)
    spatial_length = None
    if use.sum() >= 2:
        slope = np.polyfit(mid[use], np.log(correlations[use]), 1,
                           w=np.sqrt(spatial_count[use]))[0]
        if slope < 0:
            spatial_length = float(np.clip(-1 / slope, 25, 2000))
    manifest = {"release_ready": not smoke_local_matrix and len(samples) == 206
                and max_sites is None,
    "scope": {"chromosomes": chroms, "max_sites": max_sites,
                          "kind": "provisional" if smoke_local_matrix or max_sites is not None else
                          ("all_autosomes" if chromosomes is None else "chromosome_subset")},
                "method": "site moment beta-binomial fit; median tau in 99 mean bins",
                "sample_count": len(samples), "minimum_observed": minimum,
                "reference_fasta": str(reference_fasta),
                "reference_sha256": digest(reference_fasta),
                "verified_sources_sha256": None if smoke_local_matrix else digest(source_dir / "verified_sources.json"),
                "sample_mean_depth": (sample_depth_sum / max(counts["all"], 1)).tolist(),
                "matrix_sha256": digest(matrix), "source_dir": str(source_dir),
                "counts": counts, "tau_by_mean_bin": tau_by_bin, "empty_mean_bins": empty,
                "spatial_correlation": {"gap_edges_bp": spatial_edges.tolist(),
                                        "mean_residual_correlation": correlations.tolist(),
                                        "pairs": spatial_count.tolist(),
                                        "fitted_length_bp": spatial_length,
                                        "validated": False},
                "parquet_files": output_files}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
