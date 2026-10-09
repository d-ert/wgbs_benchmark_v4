"""Independent bounded review checks; never modifies callers or benchmark data."""
from pathlib import Path
import ast
import csv
import hashlib
import importlib.util
import json
import math
import sys

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "vendor310"))
import polars as pl


def rows(path):
    with open(ROOT / path) as handle:
        return list(csv.DictReader(handle))


def load_function(path, name):
    tree = ast.parse(path.read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    env = {"np": np, "sp_stats": stats}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), env)
    return env[name]


result = {}
source = ROOT / "results/default_search_by_20261007/source_snapshot/epykit"
combine = load_function(source / "dmr.py", "_stouffer_combine_signed")
result["zero_p_probe"] = {
    "p_zero_and_half": combine(np.array([0., .5]), np.ones(2)),
    "p_tiny_and_half": combine(np.array([np.finfo(float).tiny, .5]), np.ones(2)),
    "all_zero_is_nan": bool(np.isnan(combine(np.zeros(2), np.ones(2)))),
}

for filename in ["dmc.py", "dmr.py"]:
    baseline = ROOT / "results/genome_baseline/source_snapshot/epykit" / filename
    current = Path("/scratch/epykit3/src/epykit") / filename
    result[f"current_matches_baseline_{filename}"] = hashlib.sha256(baseline.read_bytes()).hexdigest() == hashlib.sha256(current.read_bytes()).hexdigest()

data = rows("results/beta_binomial_20261007_full/comparison/all_tool_region_comparison.csv")
region_checks = 0
for row in data:
    if row["dataset"] != "genome_autosomes_signal":
        continue
    calls, truth, matched = (float(row[k]) for k in ["called_regions", "true_regions", "matched_regions"])
    if calls and truth:
        for key, expected in [("region_precision", matched/calls), ("region_recall", matched/truth), ("region_f1", 2*matched/(calls+truth))]:
            assert math.isclose(float(row[key]), expected, abs_tol=1e-12)
        region_checks += 1
result["region_rows_arithmetic_checked"] = region_checks

result["tails_at_statistic_30"] = {"F_1_8": float(stats.f.sf(30, 1, 8)), "chi2_1": float(stats.chi2.sf(30, 1))}
result["separate_22_chromosome_tests_global_null_fdr_example"] = 1 - .95**22
result["two_sided_mirror_invariant_exact_minima"] = {"six_pairs": 2/64, "five_plus_five": 2/math.comb(10,5)}
depth = np.array([5., 5., 5., 5., 100.])
rho = .1
design = 1+(depth-1)*rho
pooled_variance_factor = np.sum(depth*design)/depth.sum()**2
optimal_variance_factor = 1/np.sum(depth/design)
result["known_rho_inverse_variance_example"] = {
    "depths": depth.tolist(), "rho": rho,
    "pooled_to_optimal_variance_ratio": float(pooled_variance_factor/optimal_variance_factor),
    "effective_read_weights": (depth/design).tolist(),
    "note": "Exact variance comparison for independent replicate proportions at a common mean and known rho; not a fitted-method benchmark.",
}

# Smoothing diagnostic: independent, homogeneous CpGs with the same depth.
# Averaging counts leaves pseudo-depth N but reduces variance by K. This is
# valid underdispersion induced by preprocessing, not biological underdispersion.
rng = np.random.default_rng(20261009)
k, n, rho, mu = 10, 20, .1, .5
latent = rng.beta(mu*(1-rho)/rho, (1-mu)*(1-rho)/rho, (20000, 10, k))
averaged_m = rng.binomial(n, latent).mean(axis=2)
pearson = np.zeros(len(averaged_m))
for group in [averaged_m[:, :5], averaged_m[:, 5:]]:
    fitted = group.mean(axis=1)/n
    pearson += (((group-n*fitted[:, None])**2)/n).sum(axis=1)/(fitted*(1-fitted))
raw_phi = pearson/8
result["average_count_smoothing_floor_example"] = {
    "K": k, "N": n, "rho": rho,
    "theoretical_pseudo_count_scale": (1+(n-1)*rho)/k,
    "empirical_mean_raw_scale": float(raw_phi.mean()),
    "fraction_raw_scale_below_one": float((raw_phi < 1).mean()),
    "variance_inflation_from_floor_against_known_model": k/(1+(n-1)*rho),
    "scope": "Controlled homogeneous independent-CpG model; not an estimate of the effect in GSE263850.",
}

# Recompute complete-family BY directly, without calling adjust_selected.
spec = importlib.util.spec_from_file_location("region_search", source / "_region_search.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
base = ROOT / "results/default_search_by_20261007/simulated/genome_autosomes_signal/epykit/store/.cache/dmc/lr"
paths = sorted(base.glob("chrom=*.parquet"))
family = sum(module.count_intervals(pl.read_parquet(p, columns=["pos"])["pos"].to_numpy(), 5, 50, 500) for p in paths)
cache_paths = list(base.glob(".dmr_chain_merge*.parquet"))
assert len(cache_paths) == 1
cache = pl.read_parquet(cache_paths[0])
p = cache["combined_pvalue"].to_numpy()
order = np.argsort(p)
# Euler-Maclaurin harmonic approximation is essentially machine-accurate at this family size.
harmonic = math.log(family)+np.euler_gamma+1/(2*family)-1/(12*family**2)
q_sorted = np.minimum(1, np.minimum.accumulate((p[order]*family*harmonic/np.arange(1,len(p)+1))[::-1])[::-1])
q = np.empty_like(p)
q[order] = q_sorted
assert np.allclose(q, cache["combined_qvalue"].to_numpy(), rtol=1e-11, atol=1e-14)
result["signal_search_correction"] = {"chromosomes": len(paths), "family": family, "harmonic": harmonic, "candidates": len(p), "calls": int((q <= .05).sum())}

real_base = ROOT / "results/GSE263850_paper_20261009/real/epykit"
real_caches = list((real_base / "store/.cache/dmc").rglob(".dmr_chain_merge*.parquet"))
assert len(real_caches) == 1
real = pl.read_parquet(real_caches[0])
result["GSE263850_candidates"] = {"before_q": len(real), "after_q": int((real["combined_qvalue"] <= .05).sum()), "config": json.loads((real_base / "runner_config.json").read_text())}

output = Path(__file__).with_name("checks.json")
output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
print(json.dumps(result, indent=2, allow_nan=False))
