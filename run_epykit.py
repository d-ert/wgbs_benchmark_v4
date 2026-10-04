#!/usr/bin/env python3
"""Run the local epykit lr/chain_merge pair on a v3 sample sheet."""

from __future__ import annotations

import json
import os
import hashlib
import sys
import time
from pathlib import Path

import pandas as pd
import polars as pl


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python run_epykit.py DATASET OUTPUT_DIR")
    dataset, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    profile = os.environ.get("WGBS_EPYKIT_PROFILE", "baseline")
    if profile not in {"baseline", "historical_site_smooth", "permutation_max_t", "permutation_region"}:
        raise ValueError(f"Unknown epykit profile: {profile}")
    dmc_kwargs = dict(test="lr", allow_n1=False, dispersion="eb", reference="adaptive",
                      fdr_method="fdr_bh", backend="sequential", smoothing=False)
    dispersion_override = os.environ.get("WGBS_EPYKIT_DISPERSION")
    if dispersion_override:
        if dispersion_override not in {"eb", "site", "bb_site", "bb_eb"}:
            raise ValueError("Unsupported benchmark dispersion profile")
        dmc_kwargs["dispersion"] = dispersion_override
    dmr_kwargs = dict(method="chain_merge", alpha=.05, min_abs_meth_diff=.1,
                      minlen_bp=50, min_cpgs=5, dis_merge_bp=500, pct_sig=.5,
                      min_mean_qvalue=.05, use_q_for_sig=False, empirical_fdr=False)
    if profile == "historical_site_smooth":
        dmc_kwargs.update(dispersion="site", smoothing=True, smoothing_span_bp=500)
        dmr_kwargs.update(alpha=1e-5, min_abs_meth_diff=0., min_cpgs=3, dis_merge_bp=100)
    if profile.startswith("permutation_"):
        dmr_kwargs.update(empirical_fdr=True, fdr_method=profile.removeprefix("permutation_"),
                          n_perm=int(os.environ.get("WGBS_EPYKIT_PERMUTATIONS", "100")),
                          perm_seed=20260930, perm_n_jobs=1, min_mean_qvalue=None)
    if (out / "store").exists():
        raise ValueError("Use a fresh output directory to avoid reusing an old profile cache")
    manifest_path = dataset / "manifest.json"
    dataset_manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    strand_collapsed = dataset_manifest.get("kind") == "real" and dataset_manifest.get("study") == "GSE64177"
    (out / "runner_config.json").write_text(json.dumps({
        "profile": profile, "dmc": dmc_kwargs, "dmr": dmr_kwargs,
        "coverage_contract": "coverage >=1 in every sample; no high-coverage trimming",
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "coordinate_system": "1-based inclusive regions; native epykit end minus one",
        "assembly": dataset_manifest.get("assembly", "hg38"),
        "input_strand_collapsed": strand_collapsed,
        "region_qvalue_note": "permutation output requires empirical validation" if dmr_kwargs["empirical_fdr"] else "asymptotic ranking only; not demonstrated region FDR control"
    }, indent=2) + "\n")
    sheet = pd.read_csv(dataset / "samples.tsv", sep="\t")
    if set(sheet.group) != {"control", "treatment"}:
        raise ValueError("Invalid sample groups")
    os.environ.setdefault("NUMBA_CACHE_DIR", str((out / "numba_cache").resolve()))
    import epykit as ep
    if strand_collapsed:
        # GEO supplies one count row per CpG coordinate. The public reader
        # does not expose convert_sample(merge_strands=False), so constrain
        # this temporary adapter to the GEO real-data manifest.
        import epykit.convert as conversion
        original_convert = conversion.convert_sample
        def convert_collapsed(*args, **kwargs):
            kwargs["merge_strands"] = False
            return original_convert(*args, **kwargs)
        conversion.convert_sample = convert_collapsed
    t0 = time.monotonic()
    sample_csv = out / "samples.csv"
    sheet.to_csv(sample_csv, index=False)
    assembly = dataset_manifest.get("assembly", "hg38")
    md = ep.read_bismark(str(sample_csv), treatment_group="treatment",
                         control_group="control", assembly=assembly,
                         store_dir=str(out / "store"))
    ep.pp.set_unite_type(md, type="intersect")
    t_load = time.monotonic()
    ep.tl.dmc(md, **dmc_kwargs)
    dmc = md.dmc
    if dmc is None:
        raise ValueError("epykit returned no active DMC result")
    dmc.select([c for c in ["chrom", "pos", "pvalue", "qvalue", "meth_diff"]
                if c in dmc.columns]).write_csv(out / "dml.tsv", separator="\t")
    t_dml = time.monotonic()
    ep.tl.dmr(md, **dmr_kwargs)
    (out / "native_metadata.json").write_text(json.dumps({"dmc": md.uns.get("dmc"),
        "dmr": md.uns.get("dmr_params"), "epykit_version": getattr(ep, "__version__", "unknown")},
        indent=2, default=str) + "\n")
    dmr = md.uns["dmr"]
    if isinstance(dmr, dict):
        dmr = dmr["frame"] if "frame" in dmr else next(iter(dmr.values()))
    dmr = dmr.with_columns((pl.col("end") - 1).alias("end"))
    dmr.write_csv(out / "dmr_candidates.tsv", separator="\t")
    significance = "empirical_qvalue" if dmr_kwargs["empirical_fdr"] else ("combined_qvalue" if "combined_qvalue" in dmr.columns else "qvalue")
    if significance in dmr.columns:
        dmr = dmr.filter(pl.col(significance).is_not_null() & (pl.col(significance) <= 0.05))
    dmr.write_csv(out / "dmr.tsv", separator="\t")
    (out / "timing.json").write_text(json.dumps({"read_filter_s": t_load - t0,
                                                 "dml_s": t_dml - t_load,
                                                 "dmr_s": time.monotonic() - t_dml}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
