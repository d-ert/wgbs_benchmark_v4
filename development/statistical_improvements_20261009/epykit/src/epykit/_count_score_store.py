"""Experimental count-score store: bounded prior training and streamed fitting."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

from . import _cache, _count_blocks, _count_score, _dispersion_prior
from ._dmc_store import DMCStore, _chrom_filename

ALGORITHM_REVISION = "bb-score-store-v1-complete-family"


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def _coordinate_hash(positions, chrom, seed):
    salt = int.from_bytes(hashlib.blake2b(str(chrom).encode(), digest_size=8).digest(), "little")
    with np.errstate(over="ignore"):
        z = (
            np.asarray(positions, dtype=np.uint64)
            ^ np.uint64(salt)
            ^ np.uint64(seed & ((1 << 64) - 1))
        )
        z = z + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def _training_coordinates(md, cfg, samples, chromosomes, intersect):
    cap = max(1, (cfg.score_prior_max_groups + 1) // 2)
    selected = []
    for chrom in chromosomes:
        sites = _count_blocks.canonical_sites(md.store, chrom, samples, intersect=intersect)
        pos = sites["pos"].to_numpy()
        if not len(pos):
            continue
        keys = _coordinate_hash(pos, chrom, cfg.score_prior_seed)
        indices = np.argsort(keys, kind="stable")[:cap]
        selected.extend((int(keys[i]), str(chrom), int(pos[i])) for i in indices)
        selected = sorted(selected)[:cap]
    grouped = {}
    for _, chrom, pos in selected:
        grouped.setdefault(chrom, []).append(pos)
    return {c: np.array(sorted(p), dtype=np.int64) for c, p in grouped.items()}


def _fallback_prior(reason):
    weights = np.full(len(_dispersion_prior.RHO_GRID), 1 / len(_dispersion_prior.RHO_GRID)).tolist()
    return dict(
        rho_grid=_dispersion_prior.RHO_GRID.tolist(),
        mean_bin_edges=_dispersion_prior.MEAN_BIN_EDGES.tolist(),
        rho_weights_by_bin=[weights] * 4,
        pooled_rho_weights=weights,
        algorithm_revision=_dispersion_prior.ALGORITHM_REVISION,
        n_training_groups=0,
        training_status="unavailable: " + reason,
        limitations="No usable training groups; uniform grid fallback, experimental only",
    )


def _train(md, cfg, samples, chromosomes, intersect, x, n_case, block_size):
    selected = _training_coordinates(md, cfg, samples, chromosomes, intersect)
    ms, ns = [], []
    for chrom, positions in selected.items():
        for _, m, n in _count_blocks.iter_count_blocks(
            md.store, chrom, samples, positions, block_size=block_size
        ):
            ms.append(m)
            ns.append(n)
    if not ms:
        return _fallback_prior("empty coordinate family")
    m, n = np.vstack(ms), np.vstack(ns)
    has_groups = ((n[:, :n_case] > 0).sum(axis=1) >= 2) | ((n[:, n_case:] > 0).sum(axis=1) >= 2)
    if not has_groups.any():
        prior = _fallback_prior("no group with two covered replicates")
    else:
        # Integrate a common group mean for simple designs, condition on
        # per-sample fitted means for nuisance designs, recording approximation.
        simple = np.allclose(x[:n_case], x[0]) and np.allclose(x[n_case:], x[n_case])
        try:
            prior = _dispersion_prior.fit_dispersion_prior(
                m,
                n,
                n_case,
                design=None if simple else x,
                max_training=cfg.score_prior_max_groups,
                seed=cfg.score_prior_seed,
            )
        except ValueError as exc:
            if not str(exc).startswith("Dispersion training requires two covered replicates"):
                raise
            prior = _fallback_prior("no usable conditional-regression training group")
    prior.update(
        samples_case=samples[:n_case],
        samples_control=samples[n_case:],
        selected_coordinates={c: p.tolist() for c, p in selected.items()},
        coordinate_sampling="seeded coordinate hash; independent of counts, group effects and truth",
        design=x.tolist(),
    )
    return prior


def _frame(chrom, positions, strands, m, n, x, c, fit, cfg, n_case):
    covered = n > 0
    n_a, n_b = covered[:, :n_case].sum(axis=1), covered[:, n_case:].sum(axis=1)
    eligible = (n_a >= max(2, cfg.min_samples_treatment)) & (n_b >= max(2, cfg.min_samples_control))
    p = np.where(eligible, fit["pvalue"], np.nan)
    beta = np.divide(m, n, out=np.zeros(m.shape), where=covered)

    def descriptive_mean(lo, hi):
        support = covered[:, lo:hi].sum(axis=1)
        return np.divide(
            beta[:, lo:hi].sum(axis=1), support, out=np.full(len(m), np.nan), where=support > 0
        )

    a, b = descriptive_mean(0, n_case), descriptive_mean(n_case, n.shape[1])
    meth_a, meth_b = m[:, :n_case].sum(axis=1), m[:, n_case:].sum(axis=1)
    unmeth_a, unmeth_b = n[:, :n_case].sum(axis=1) - meth_a, n[:, n_case:].sum(axis=1) - meth_b
    odds = np.log2((meth_a + 0.5) / (unmeth_a + 0.5)) - np.log2((meth_b + 0.5) / (unmeth_b + 0.5))
    data = dict(
        chrom=[chrom] * len(m),
        pos=positions,
        strand=strands,
        n_case=n_a,
        n_control=n_b,
        mean_beta_case=a,
        mean_beta_control=b,
        pvalue=p,
        pvalue_for_correction=np.nan_to_num(p, nan=1.0),
        log2_odds_ratio_pooled=odds,
        log2_odds_ratio=np.full(len(m), np.nan),
        coef_treatment_log2=fit["coefficients"] @ c / np.log(2.0),
        meth_diff=fit["effect"],
        meth_diff_raw=a - b,
        meth_diff_ci_lo=fit["ci_lo"],
        meth_diff_ci_hi=fit["ci_hi"],
        score_stat=fit["statistic"],
        score_variance=fit["variance"],
        score=fit["score"],
        score_valid_samples=fit["valid_counts"],
        score_converged=fit["converged"],
        score_full_converged=fit["full_converged"],
        score_estimable=eligible & fit["estimable"],
        interval_available=eligible & fit["interval_available"],
    )
    for field in (
        "rho",
        "rho_case",
        "rho_control",
        "rho_uncertainty",
        "rho_case_uncertainty",
        "rho_control_uncertainty",
        "dispersion_iterations",
        "dispersion_converged",
        "dispersion_fallback",
    ):
        data[field] = fit[field]
    # Canonical schema stays consistent across engine stores.
    from .dmc import _EMPTY_SCHEMA

    return pl.DataFrame(data).with_columns(
        [pl.col(name).cast(dtype) for name, dtype in _EMPTY_SCHEMA.items()]
    )


def _cache_valid(path, fingerprint):
    try:
        store = DMCStore.open(path)
        manifest = store.manifest
        if manifest.get("fingerprint") != fingerprint or not manifest.get("complete"):
            return None
        if _sha(store.path / manifest["prior_file"]) != manifest["prior_sha256"]:
            return None
        for entry in manifest["chroms"]:
            p = store.path / _chrom_filename(entry["name"])
            if _sha(p) != entry["sha256"]:
                return None
        return store
    except (OSError, ValueError, KeyError):
        return None


def run_count_score_store(md, cfg, *, design, contrast, out_dir=None):
    """Own bounded observed-count training, fitting and exact global correction.

    Design rows follow md.treatment_ids + md.control_ids. Explicit out_dir is
    for isolated direct callers/benchmarks; otherwise use a content-keyed cache.
    Model work blocks also cap dense-information memory at approximately 128MiB.
    Existing correction currently uses global vectors, never a genome-wide
    sample matrix. The planned disk-backed correction can replace that step.
    """
    cfg.validate()
    cfg.validate_resolved()
    samples = list(md.treatment_ids) + list(md.control_ids)
    n_case = len(md.treatment_ids)
    if min(n_case, len(md.control_ids)) < 2 or len(set(samples)) != len(samples):
        raise ValueError("bb_score requires two unique biological replicates per group")
    if len(samples) != md.obs.height:
        raise ValueError(
            "bb_score supports the declared two groups only; extra samples are unassigned"
        )
    x, c = np.asarray(design, dtype=float), np.asarray(contrast, dtype=float)
    if x.ndim != 2 or x.shape[0] != len(samples) or x.shape[1] >= len(samples):
        raise ValueError("bb_score design needs sample-aligned rows and residual replication")
    if c.ndim == 2 and c.shape[0] == 1:
        c = c[0]
    if c.shape != (x.shape[1],) or not np.isfinite(c).all() or not np.any(c):
        raise ValueError("bb_score supports one finite nonzero scalar contrast")
    if not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("bb_score design must be finite and full rank")
    from .dmc import resolve_dmc_chromosomes

    chromosomes = resolve_dmc_chromosomes(
        md.store, cfg.chromosomes, canonical_only=cfg.canonical_only
    )
    intersect = md.uns.get("unite", {}).get("type") == "intersect"
    inputs = []
    for chrom in chromosomes:
        for sample in samples:
            for p in _count_blocks.sample_parts(md.store, chrom, sample):
                stat = p.stat()
                inputs.append((str(p.resolve()), stat.st_size, stat.st_mtime_ns))
    params = dict(
        algorithm_revision=ALGORITHM_REVISION,
        score_revision=_count_score.ALGORITHM_REVISION,
        prior_revision=_dispersion_prior.ALGORITHM_REVISION,
        blocks_revision=_count_blocks.ALGORITHM_REVISION,
        tail_revision="normal-plugin-v1",
        samples=samples,
        n_case=n_case,
        design=x.tolist(),
        contrast=c.tolist(),
        inputs=inputs,
        chromosomes=chromosomes,
        intersect=intersect,
        min_samples_case=cfg.min_samples_treatment,
        min_samples_control=cfg.min_samples_control,
        fdr_method=cfg.fdr_method,
        **cfg.score_options(),
    )
    fingerprint = hashlib.sha256(
        json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    root = (
        Path(out_dir)
        if out_dir is not None
        else Path(md.analysis_root or Path(md.store).parent) / ".cache/dmc/bb_score" / fingerprint
    )
    cached = _cache_valid(root, fingerprint)
    if cached is not None:
        cached._manifest["_cache_hit"] = True  # transient outcome provenance; not written
        return cached
    if out_dir is not None and root.exists() and any(root.iterdir()):
        raise ValueError(
            "bb_score output directory is nonempty and has no valid matching completed cache"
        )
    root.mkdir(parents=True, exist_ok=True)
    block_size = min(
        cfg.score_block_size, max(1, 134217728 // (8 * (8 * x.shape[1] ** 2 + 12 * len(samples))))
    )
    prior = _train(md, cfg, samples, chromosomes, intersect, x, n_case, block_size)
    prior_path = root / "dispersion_prior.json"
    _cache.write_json(prior_path, prior)
    prior_hash = _sha(prior_path)
    entries, total, untestable = [], 0, 0
    for chrom in chromosomes:
        sites = _count_blocks.canonical_sites(md.store, chrom, samples, intersect=intersect)
        if not len(sites):
            continue
        pos = sites["pos"].to_numpy()
        strands = sites["strand"].to_numpy()
        target = root / _chrom_filename(chrom)
        temporary = target.with_suffix(".parquet.tmp")
        writer, offset = None, 0
        try:
            for positions, m, n in _count_blocks.iter_count_blocks(
                md.store, chrom, samples, pos, block_size=block_size
            ):
                fit = _count_score.test_counts(
                    m,
                    n,
                    x,
                    c,
                    prior=prior,
                    n_case=n_case,
                    dispersion_mode=cfg.score_dispersion,
                    intervals=cfg.score_ci,
                )
                frame = _frame(
                    chrom,
                    positions,
                    strands[offset : offset + len(m)],
                    m,
                    n,
                    x,
                    c,
                    fit,
                    cfg,
                    n_case,
                )
                table = frame.to_arrow()
                if writer is None:
                    writer = pq.ParquetWriter(temporary, table.schema, compression="zstd")
                writer.write_table(table, row_group_size=block_size)
                untestable += int(frame["pvalue"].is_nan().sum())
                offset += len(m)
        finally:
            if writer is not None:
                writer.close()
        temporary.replace(target)
        entries.append(
            dict(name=chrom, n_sites=offset, prior_file=prior_path.name, prior_sha256=prior_hash)
        )
        total += offset
    manifest = dict(
        test="bb_score",
        algorithm_revision=ALGORITHM_REVISION,
        fingerprint=fingerprint,
        parameters=params,
        chroms=entries,
        total_sites=total,
        hypothesis_count=total,
        untestable_sites=untestable,
        complete=False,
        prior_file=prior_path.name,
        prior_sha256=prior_hash,
        actual_block_size=block_size,
        inference_status="experimental",
        missing_hypothesis_policy="retain p=1 for correction, NaN statistic status",
    )
    store = DMCStore(path=root, test="bb_score", _manifest=manifest)
    store.write_manifest(manifest.copy())
    if total:
        from .dmc import apply_multiple_testing_correction

        store = apply_multiple_testing_correction(
            store, method=cfg.fdr_method, pvalue_col="pvalue_for_correction"
        )
    completed = store.manifest
    for entry in completed["chroms"]:
        entry["sha256"] = _sha(root / _chrom_filename(entry["name"]))
    completed["complete"] = True
    store.write_manifest(completed)
    return store
