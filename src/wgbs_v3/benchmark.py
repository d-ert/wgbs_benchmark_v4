"""Sequential six-tool benchmark with one frozen input and score contract."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import signal
import fcntl
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import psutil

from .score import score
from .provenance import snapshot, write_json

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ("epykit", "DSS", "methylKit", "BSmooth", "dmrseq", "DMRcate")


def run_monitored(command: list[str], log: Path, env: dict, memory_gib: float = 48,
                  timeout_seconds: float = 86400) -> dict:
    start = time.monotonic()
    with log.open("wb") as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        peak = 0
        termination_reason = None
        parent = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                family = [parent] + parent.children(recursive=True)
                rss = sum(p.memory_info().rss for p in family if p.is_running())
                peak = max(peak, rss)
                if rss > memory_gib * 2**30 or time.monotonic() - start > timeout_seconds:
                    termination_reason = "memory_limit" if rss > memory_gib * 2**30 else "timeout"
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    break
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(.2)
        code = process.wait()
    return {"exit_code": code, "wall_seconds": time.monotonic() - start,
            "peak_process_tree_rss_bytes": peak, "command": command,
            "termination_reason": termination_reason,
            "memory_limit_gib": memory_gib, "timeout_seconds": timeout_seconds,
            "memory_measure": "sum of sampled process-tree RSS; shared pages can be counted repeatedly"}


def benchmark(dataset: Path, results: Path, tools: tuple[str, ...] = TOOLS,
              allow_provisional: bool = False, parallel: bool = False,
              ep_source: Path | None = None, overrides: dict | None = None) -> dict:
    meta = json.loads((dataset / "manifest.json").read_text())
    if not meta["calibration_release_ready"] and not allow_provisional:
        raise ValueError("Benchmark release requires 206-sample calibration")
    results.mkdir(parents=True, exist_ok=True)
    if (results / "manifest.json").exists() or any((results / t).exists() for t in tools):
        raise ValueError("Refusing to overwrite a previous run; use a fresh output directory")
    if parallel:
        raise ValueError("Timed comparisons must be sequential")
    if len(set(tools)) != len(tools) or not set(tools).issubset(TOOLS):
        raise ValueError("Unknown or duplicate tools")
    # One timed caller across all study invocations on this workspace.
    timing_lock = (ROOT / "results/.timing.lock").open("a")
    fcntl.flock(timing_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    ep_source = (ep_source or Path(os.environ.get("WGBS_EPYKIT_SOURCE", "/scratch/epykit3/src"))).resolve()
    snapshot(ROOT, ep_source, dataset, results)
    env = dict(os.environ)
    env["PYTHONPATH"] = ":".join([str(ROOT / "vendor310"), str(ROOT / "src"),
                                     str(ep_source), env.get("PYTHONPATH", "")])
    env.update({k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "POLARS_MAX_THREADS")})
    env.update(overrides or {})
    env["MPLCONFIGDIR"] = str(ROOT / "data/mplcache")
    output = {"dataset": str(dataset), "provisional": allow_provisional,
              "tools": {}, "coordinate_system": "1-based inclusive",
              "dataset_manifest": meta, "epykit_source": str(ep_source),
              "thread_budget": 1, "statistical_units": {
                  "epykit": "native CpG tests; region inference depends on recorded profile",
                  "DSS": "CpG BH; heuristic regions", "methylKit": "CpG BH and fixed-window BH",
                  "BSmooth": "t-statistic heuristic; no nominal FDR", "dmrseq": "native region q",
                  "DMRcate": "CpG BH gate and smoothed region statistic"}}
    write_json(results / "manifest.json", output)
    def run_one(tool: str):
        if tool not in TOOLS:
            raise ValueError(f"Unknown tool: {tool}")
        target = results / tool
        target.mkdir(exist_ok=True)
        if tool == "epykit":
            cmd = [sys.executable, str(ROOT / "run_epykit.py"), str(dataset), str(target)]
        else:
            cmd = ["Rscript", str(ROOT / "run_r_tool.R"), tool, str(dataset), str(target)]
        monitoring = run_monitored(cmd, target / "run.log", env)
        record = {"monitoring": monitoring, "status": "failed" if monitoring["exit_code"] else "ok"}
        record["configuration_overrides"] = {key: value for key, value in env.items()
                                              if key.startswith("WGBS_EPYKIT_") or
                                              key.startswith("WGBS_DSS_") or
                                              key.startswith("WGBS_BSMOOTH_") or
                                              key.startswith("WGBS_DMRCATE_") or
                                              key.startswith("WGBS_METHYLKIT_")}
        if monitoring["exit_code"] == 0:
            dml = target / "dml.tsv"
            dmr = target / "dmr.tsv"
            try:
                record["score"] = score(dataset, dml if dml.exists() else None,
                                        dmr if dmr.exists() else None, target / "score.json")
            except Exception as exc:
                record["status"] = "scoring_failed"
                record["scoring_error"] = f"{type(exc).__name__}: {exc}"
        return tool, record

    if parallel:
        with ThreadPoolExecutor(max_workers=len(tools)) as pool:
            futures = [pool.submit(run_one, tool) for tool in tools]
            for future in as_completed(futures):
                tool, record = future.result()
                output["tools"][tool] = record
                write_json(results / "manifest.json", output)
    else:
        for tool in tools:
            name, record = run_one(tool)
            output["tools"][name] = record
            write_json(results / "manifest.json", output)
    timing_lock.close()
    return output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", type=Path)
    ap.add_argument("results", type=Path)
    ap.add_argument("--tools", default=",".join(TOOLS))
    ap.add_argument("--allow-provisional", action="store_true")
    ap.add_argument("--parallel", action="store_true")
    ap.add_argument("--epykit-source", type=Path)
    args = ap.parse_args()
    result = benchmark(args.dataset, args.results, tuple(args.tools.split(",")), args.allow_provisional,
                       args.parallel, ep_source=args.epykit_source)
    if any(r["status"] != "ok" for r in result["tools"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
