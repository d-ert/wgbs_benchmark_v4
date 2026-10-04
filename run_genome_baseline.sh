#!/usr/bin/env bash
cd /scratch/wgbs_benchmark_v4 || exit 1
python3 -u benchmark.py baseline --whole-genome --name genome_baseline >> results/genome_baseline.tmux.log 2>&1
status=$?
printf '\nBenchmark exited with status %s at %s\n' "$status" "$(date -Is)" >> results/genome_baseline.tmux.log
exit "$status"
