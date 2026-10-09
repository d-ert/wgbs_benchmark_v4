"""Run comparisons only after the complete three-dataset benchmark succeeds."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time

HERE = Path(__file__).parent.resolve()
ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT/'results/beta_binomial_20261007'
env = dict(os.environ)
env['PYTHONPATH'] = ':'.join([str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)])
env['MPLCONFIGDIR'] = '/tmp/epykit_bb_analysis_mpl'
env['OPENBLAS_NUM_THREADS'] = '1'
while True:
    path = RUN/'run.json'
    if path.exists():
        record = json.loads(path.read_text())
        if record['status']=='failed':
            raise RuntimeError(f'Benchmark failed: {record.get("error")}')
        if record['status']=='complete':
            break
    time.sleep(30)
print('Complete benchmark confirmed; analyzing all tool comparisons...',flush=True)
subprocess.run([sys.executable,str(HERE/'plot_validation.py')],env=env,check=True)
subprocess.run([sys.executable,str(HERE/'analyze.py')],env=env,check=True)
subprocess.run(['Rscript',str(HERE/'validate_paper.R')],env=env,check=True)
subprocess.run([sys.executable,str(HERE/'report.py')],env=env,check=True)
print(f'Comparison artifacts ready: {RUN}/comparison/report.html',flush=True)
