"""Caller-independent DML and DMR evaluation against immutable latent truth."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import average_precision_score, roc_auc_score, matthews_corrcoef


def truth_cpg(dataset: Path) -> pd.DataFrame:
    parts = sorted((dataset / "truth_parts").glob("*.parquet"))
    if not parts:
        raise ValueError("Missing CpG truth")
    frame = pd.concat((pq.read_table(p).to_pandas() for p in parts), ignore_index=True)
    if frame.duplicated(["chrom", "pos"]).any():
        raise ValueError("Duplicate chromosome/CpG truth key")
    if any(np.any(np.diff(g.pos.to_numpy()) <= 0) for _, g in frame.groupby("chrom", sort=False)):
        raise ValueError("Truth must be strictly sorted within chromosomes")
    return frame


def common_eligibility(dataset: Path, truth: pd.DataFrame, verify_files: bool = False) -> np.ndarray:
    cached = sorted((dataset / "eligibility_parts").glob("*.parquet"))
    if cached and not verify_files:
        frame = pd.concat((pq.read_table(p).to_pandas() for p in cached), ignore_index=True)
        if not frame[["chrom", "pos"]].equals(truth[["chrom", "pos"]]):
            raise ValueError("Cached eligibility keys differ from truth")
        return frame.eligible.to_numpy(dtype=bool)
    sample_sheet = pd.read_csv(dataset / "samples.tsv", sep="\t")
    index = pd.MultiIndex.from_frame(truth[["chrom", "pos"]])
    retain = np.ones(len(index), dtype=bool)
    for row in sample_sheet.itertuples(index=False):
        part = pd.read_csv(row.path, sep="\t", header=None, usecols=[0, 1, 4, 5],
                           names=["chrom", "pos", "m", "u"], compression="gzip")
        values = part[["pos", "m", "u"]].to_numpy(dtype=float)
        if not np.isfinite(values).all() or np.any(values != np.floor(values)) or (part[['m', 'u']] < 0).any().any():
            raise ValueError(f"Noninteger or invalid counts: {row.sample_id}")
        keys = pd.MultiIndex.from_frame(part[["chrom", "pos"]])
        if keys.has_duplicates:
            raise ValueError(f"Duplicate sample CpG: {row.sample_id}")
        if not keys.isin(index).all():
            raise ValueError(f"Sample CpG outside truth: {row.sample_id}")
        retain &= index.isin(keys[(part.m + part.u).to_numpy() >= 1])
    if cached:
        expected = common_eligibility(dataset, truth, verify_files=False)
        if not np.array_equal(retain, expected):
            raise ValueError("Cached eligibility differs from raw Bismark files")
    return retain


def load_calls(path: Path, kind: str) -> pd.DataFrame:
    frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")
    if "chr" in frame.columns and "chrom" not in frame.columns:
        frame = frame.rename(columns={"chr": "chrom"})
    if kind == "dml":
        if "start" in frame.columns and "pos" not in frame.columns:
            frame = frame.rename(columns={"start": "pos"})
        required = {"chrom", "pos", "pvalue"}
        key = ["chrom", "pos"]
    else:
        required = {"chrom", "start", "end"}
        key = ["chrom", "start", "end"]
    if not required.issubset(frame.columns) or frame[key].isna().any().any():
        raise ValueError(f"Malformed {kind} calls: {path}")
    if frame.duplicated(key).any():
        raise ValueError(f"Duplicate caller results: {path}")
    return frame


def dml_metrics(truth: pd.DataFrame, eligible: np.ndarray, calls: pd.DataFrame,
                q_threshold: float = 0.05) -> dict:
    keys = pd.MultiIndex.from_frame(truth[["chrom", "pos"]])
    ck = pd.MultiIndex.from_frame(calls[["chrom", "pos"]])
    if ck.has_duplicates:
        raise ValueError("Duplicate DML calls")
    if not ck.isin(keys).all():
        raise ValueError("Caller reported CpG outside simulation universe")
    idx = keys.get_indexer(ck)
    if not eligible[idx].all():
        raise ValueError("Caller reported CpG outside common coverage mask")
    p = np.ones(len(truth), dtype=float)
    q = np.ones(len(truth), dtype=float)
    p[idx] = calls.pvalue.fillna(1).to_numpy(dtype=float)
    if "qvalue" in calls:
        q[idx] = calls.qvalue.fillna(1).to_numpy(dtype=float)
    if ((p < 0) | (p > 1) | (q < 0) | (q > 1)).any():
        raise ValueError("Invalid p or q values")
    y = truth.is_dml.to_numpy(dtype=bool)
    pred = (q <= q_threshold) & eligible
    pos = y & eligible
    neg = ~y & eligible
    tp, fp = int(np.sum(pred & pos)), int(np.sum(pred & neg))
    fn, tn = int(np.sum(~pred & pos)), int(np.sum(~pred & neg))
    valid = eligible & np.isfinite(p)
    scored = np.unique(p[valid]).size > 1 and y[valid].any() and (~y[valid]).any()
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": tp / (tp + fp) if tp + fp else None,
            "recall_eligible": tp / (tp + fn) if tp + fn else None,
            "recall_full_truth": tp / int(y.sum()) if y.sum() else None,
            "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.,
            "mcc": float(matthews_corrcoef(y[eligible], pred[eligible])) if eligible.any() else None,
            "false_discovery_proportion": fp / (tp + fp) if tp + fp else 0.,
            "average_precision": float(average_precision_score(y[valid], -p[valid])) if scored else None,
            "auroc": float(roc_auc_score(y[valid], -p[valid])) if scored else None,
            "n_eligible": int(eligible.sum()), "n_truth": int(y.sum()),
            "n_scored": int((p[eligible] < 1).sum())}


def _mask_from_regions(truth: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    result = np.zeros(len(truth), dtype=bool)
    for chrom, group in truth.groupby("chrom", sort=False):
        pos = group.pos.to_numpy()
        indices = group.index.to_numpy()
        for call in regions.loc[regions.chrom == chrom].itertuples():
            lo = np.searchsorted(pos, call.start)
            hi = np.searchsorted(pos, call.end, side="right")
            result[indices[lo:hi]] = True
    return result


def dmr_metrics(truth: pd.DataFrame, eligible: np.ndarray, calls: pd.DataFrame,
                overlap: float = 0.5) -> dict:
    """Reciprocal CpG overlap plus maximum-cardinality weighted matching."""
    if calls.empty:
        calls = pd.DataFrame(columns=["chrom", "start", "end"])
    calls = calls.reset_index(drop=True)
    if calls.duplicated(["chrom", "start", "end"]).any():
        raise ValueError("Duplicate DMR calls")
    coords = calls[["start", "end"]].to_numpy(dtype=float)
    if not np.isfinite(coords).all() or np.any(coords != np.floor(coords)) or np.any(coords < 1):
        raise ValueError("DMR coordinates must be 1-based integer inclusive intervals")
    if (calls.end < calls.start).any():
        raise ValueError("Invalid region boundaries")
    unknown = set(calls.chrom) - set(truth.chrom)
    if unknown:
        raise ValueError(f"Caller reported unknown chromosomes: {unknown}")
    predicted_mask = _mask_from_regions(truth, calls)
    y = truth.region_id.to_numpy() > 0
    tp_cpg = int(np.sum(predicted_mask & y & eligible))
    fp_cpg = int(np.sum(predicted_mask & ~y & eligible))
    fn_cpg = int(np.sum(~predicted_mask & y & eligible))
    true_regions = truth.loc[truth.region_id > 0].groupby("region_id", sort=False).agg(
        chrom=("chrom", "first"), start=("pos", "min"), end=("pos", "max"),
        direction=("signed_delta", lambda x: np.sign(x.mean())))
    graph = nx.Graph()
    truth_to_calls: dict[int, set[int]] = {}
    call_to_truth: dict[int, set[int]] = {}
    calls_without_eligible_cpg = 0
    for tid in true_regions.index:
        graph.add_node(("t", int(tid)), bipartite=0)
    for j in range(len(calls)):
        graph.add_node(("p", j), bipartite=1)
    for chrom, group in truth.groupby("chrom", sort=False):
        pos = group.pos.to_numpy()
        ids = group.region_id.to_numpy()
        chrom_eligible = eligible[group.index.to_numpy()]
        true_sizes = np.bincount(ids[chrom_eligible].clip(min=0))
        selected = calls.loc[calls.chrom == chrom]
        for j, call in selected.iterrows():
            lo = np.searchsorted(pos, call.start)
            hi = np.searchsorted(pos, call.end, side="right")
            pred_size = int(chrom_eligible[lo:hi].sum())
            if pred_size == 0:
                calls_without_eligible_cpg += 1
                continue
            sub_ids = ids[lo:hi][chrom_eligible[lo:hi]]
            t_ids, counts = np.unique(sub_ids[sub_ids > 0], return_counts=True)
            for tid, inter in zip(t_ids, counts):
                truth_to_calls.setdefault(int(tid), set()).add(j)
                call_to_truth.setdefault(j, set()).add(int(tid))
                true_size = int(true_sizes[tid])
                if inter / pred_size >= overlap and inter / true_size >= overlap:
                    graph.add_edge(("t", int(tid)), ("p", j), weight=float(inter / (pred_size + true_size - inter)))
    matching = nx.algorithms.matching.max_weight_matching(graph, maxcardinality=True, weight="weight")
    matched = len(matching)
    boundary_errors = []
    direction_correct = []
    for a, b in matching:
        truth_node, call_node = (a, b) if a[0] == "t" else (b, a)
        tid, j = truth_node[1], call_node[1]
        expected = true_regions.loc[tid]
        observed = calls.iloc[j]
        boundary_errors.append((abs(int(expected.start) - int(observed.start)),
                                abs(int(expected.end) - int(observed.end))))
        predicted_direction = observed.get("direction", None)
        if predicted_direction is None or pd.isna(predicted_direction):
            if "mean_meth_diff" in calls:
                predicted_direction = np.sign(observed.mean_meth_diff)
            elif "meth_diff" in calls:
                predicted_direction = np.sign(observed.meth_diff)
        if predicted_direction is not None and not pd.isna(predicted_direction):
            direction_correct.append(int(np.sign(float(predicted_direction)) == expected.direction))
    t_count = len(true_regions)
    p_count = len(calls)
    return {"matched_regions": matched, "true_regions": t_count, "called_regions": p_count,
            "region_precision": matched / p_count if p_count else None,
            "region_recall": matched / t_count if t_count else None,
            "region_f1": 2 * matched / (p_count + t_count) if p_count + t_count else None,
            "matched_direction_accuracy": float(np.mean(direction_correct)) if direction_correct else None,
            "start_boundary_mae_bp": float(np.mean([x[0] for x in boundary_errors])) if boundary_errors else None,
            "end_boundary_mae_bp": float(np.mean([x[1] for x in boundary_errors])) if boundary_errors else None,
            "fragmented_true_regions": sum(len(c) > 1 for c in truth_to_calls.values()),
            "merged_called_regions": sum(len(t) > 1 for t in call_to_truth.values()),
            "calls_without_eligible_cpg": calls_without_eligible_cpg,
            "calls_without_positive_eligible_cpg": p_count - len(call_to_truth),
            "strict_unmatched_fraction": (p_count - matched) / p_count if p_count else None,
            "null_region_fdp": (p_count - len(call_to_truth)) / p_count if p_count else 0.,
            "cpg_tp": tp_cpg, "cpg_fp": fp_cpg, "cpg_fn": fn_cpg,
            "cpg_precision": tp_cpg / (tp_cpg + fp_cpg) if tp_cpg + fp_cpg else None,
            "cpg_recall": tp_cpg / (tp_cpg + fn_cpg) if tp_cpg + fn_cpg else None,
            "cpg_recall_full_truth": tp_cpg / int(y.sum()) if y.sum() else None}


def score(dataset: Path, dml: Path | None, dmr: Path | None, out: Path) -> dict:
    truth = truth_cpg(dataset)
    mask = common_eligibility(dataset, truth)
    result = {"dataset": str(dataset), "eligible_count": int(mask.sum())}
    if dml is not None and dml.exists():
        result["dml"] = dml_metrics(truth, mask, load_calls(dml, "dml"))
    if dmr is not None and dmr.exists():
        regions = load_calls(dmr, "dmr")
        result["dmr"] = dmr_metrics(truth, mask, regions)
        result["dmr_80"] = dmr_metrics(truth, mask, regions, overlap=0.8)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    return result
