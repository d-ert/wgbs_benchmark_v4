#!/usr/bin/env bash
set -u
cd /scratch/wgbs_benchmark_v4 || exit 1
python3 development/beta_binomial_20261007/run.py > development/beta_binomial_20261007/benchmark_full.log 2>&1
epy_job_exit=$?
printf '%s\n' "$epy_job_exit" > development/beta_binomial_20261007/benchmark_full.exit
if [ "$epy_job_exit" -eq 0 ]; then
    python3 development/beta_binomial_20261007/finish.py > development/beta_binomial_20261007/analysis_full.log 2>&1
    epy_job_exit=$?
    printf '%s\n' "$epy_job_exit" > development/beta_binomial_20261007/analysis_full.exit
fi
printf '%s\n' "$epy_job_exit" > development/beta_binomial_20261007/job_full.exit
exit "$epy_job_exit"
