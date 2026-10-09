#!/usr/bin/env bash
cd /scratch/wgbs_benchmark_v4 || exit 1
python3 real_benchmark.py run --name GSE263850_paper_20261009 --profile paper_aligned > development/gse263850_20261009/full.log 2>&1
run_status=$?
printf '%s\n' "$run_status" > development/gse263850_20261009/full.exit
exit "$run_status"
