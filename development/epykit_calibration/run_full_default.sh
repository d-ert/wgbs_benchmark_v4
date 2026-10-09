#!/usr/bin/env bash
cd /scratch/wgbs_benchmark_v4 || exit 1
export PYTHONPATH=/scratch/wgbs_benchmark_v4/vendor310:/scratch/wgbs_benchmark_v4/src:/scratch/epykit-calibration/src
export MPLCONFIGDIR=/scratch/wgbs_benchmark_v4/data/mplcache
export POLARS_MAX_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python3 -u development/epykit_calibration/scripts/run_full_default_benchmark.py >> development/epykit_calibration/results/full_default.log 2>&1
status=$?
printf '\nExited with status %s at %s\n' "$status" "$(date -Is)" >> development/epykit_calibration/results/full_default.log
exit "$status"
