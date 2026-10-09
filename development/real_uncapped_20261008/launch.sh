#!/usr/bin/env bash
cd /scratch/wgbs_benchmark_v4 || exit 1
python3 development/real_uncapped_20261008/run.py > development/real_uncapped_20261008/job.log 2>&1
job_status=$?
printf '%s\n' "$job_status" > development/real_uncapped_20261008/job.exit
exit "$job_status"
