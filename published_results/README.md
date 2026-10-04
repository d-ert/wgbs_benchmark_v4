# Compact benchmark results

This directory holds small, selected outputs so the repository includes evidence
alongside the benchmark code. The original output paths are preserved under
`published_results/`; the full `results/` and `archive/` directories remain local.

- `archive/v3/results/chr1_3_reviewed/`: reviewed chr1–3 comparison, report,
  score curves, and per-tool scores.
- `archive/v3/results/calibration_comparison_v5/`: calibration comparison and
  report.
- `archive/v3/results/publication_workflow/simulated/simulation_metrics.csv`:
  historical simulation summary.
- `results/genome_baseline/simulated/`: available per-tool score files from
  the whole-genome run. These are an **interim snapshot**: the source run was
  marked `running` when copied on 2026-10-04. They do not represent a completed
  all-tool benchmark. Refresh these files only after reviewing run status.

The large input counts, per-CpG DML tables, caches, and source snapshots are
excluded. To reproduce a run, supply the inputs described in the root README
and run the benchmark locally.
