"""Run the four pending real-data callers sequentially, without resource caps."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime,timezone

ROOT = Path('/scratch/wgbs_benchmark_v4')
HERE = Path(__file__).resolve().parent
OUT = ROOT/'results/real_uncapped_20261008'
DATA = ROOT/'data/real'
sys.path[:0] = [str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)]
from wgbs_v3.benchmark import run_monitored
from wgbs_v3.provenance import dataset_inventory,write_json

TOOLS = ['DMRcate','BSmooth','dmrseq','methylKit']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()

def main():
    if OUT.exists():raise ValueError('Preserving existing results; select a fresh run name')
    assert json.loads((HERE/'smoke_verified.json').read_text())['status']=='verified'
    env = dict(os.environ)
    prefixes = ('WGBS_EPYKIT_','WGBS_DSS_','WGBS_BSMOOTH_','WGBS_DMRCATE_','WGBS_METHYLKIT_')
    assert not [k for k in env if k.startswith(prefixes)],'Unexpected statistical overrides'
    env['PYTHONPATH'] = ':'.join([str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)])
    env['MPLCONFIGDIR'] = '/tmp/wgbs_real_uncapped_mpl'
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','POLARS_MAX_THREADS']:
        env[key] = '1'
    with (ROOT/'results/.workflow.lock').open('a') as workflow,(ROOT/'results/.timing.lock').open('a') as timing:
        fcntl.flock(workflow,fcntl.LOCK_EX|fcntl.LOCK_NB)
        fcntl.flock(timing,fcntl.LOCK_EX|fcntl.LOCK_NB)
        OUT.mkdir()
        snap = OUT/'source_snapshot';snap.mkdir()
        for name in ['run.py','compare.py','run_r_tool_original.R','run_r_tool_uncapped.R','runner_changes.json','smoke_verified.json','validate_paper.R']:
            shutil.copy2(HERE/name,snap/name)
        shutil.copy2(ROOT/'analysis/annotate_real_dmrs.R',snap/'annotate_real_dmrs.R')
        record = dict(status='running',stage='input_verification',started_at=now(),dataset=str(DATA),
            study='GSE64177',assembly='hg19',samples=12,paired_model=False,tools=TOOLS,
            memory_limit_gib=None,timeout_seconds=None,thread_budget=1,
            reused={'epykit':str(ROOT/'results/default_search_by_20261007/real/epykit'),
                    'DSS':str(ROOT/'results/genome_baseline/real/DSS')},
            source_sha256={p.name:sha(p) for p in snap.iterdir() if p.is_file()})
        write_json(OUT/'run.json',record)
        expected = json.loads((ROOT/'results/genome_baseline/run.json').read_text())['contract']['real']
        assert dataset_inventory(DATA)==expected,'Real input changed'
        record['input_inventory'] = expected
        record['stage'] = 'benchmarking'
        manifests = {}
        real = OUT/'real';real.mkdir()
        for tool in TOOLS:
            target = real/tool;target.mkdir()
            record.update(active_tool=tool,updated_at=now())
            write_json(OUT/'run.json',record)
            print(f'Starting {tool}: all 12 samples, no memory or timeout cap',flush=True)
            command = ['Rscript',str(snap/'run_r_tool_uncapped.R'),tool,str(DATA),str(target)]
            try:
                monitor = run_monitored(command,target/'run.log',env,memory_gib=float('inf'),timeout_seconds=float('inf'))
                monitor.update(memory_limit_gib=None,timeout_seconds=None)
                complete = monitor['exit_code']==0 and all((target/name).exists() for name in ['dmr.tsv','timing.tsv','session_info.txt'])
                entry = dict(status='ok' if complete else 'failed',monitoring=monitor,assembly='hg19',
                             profile='recorded defaults; methylKit unused-object cleanup only')
            except Exception as error:
                entry = dict(status='failed',error=f'{type(error).__name__}: {error}')
            manifests[tool] = entry
            write_json(real/'run_manifest.json',manifests)
            print(f'{tool}: {entry["status"]}',flush=True)
        assert dataset_inventory(DATA)==expected,'Real input changed during the run'
        assert all(sha(snap/name)==digest for name,digest in record['source_sha256'].items())
        failed = [t for t,r in manifests.items() if r['status']!='ok']
        record.update(stage='comparison',active_tool=None,failed_tools=failed,updated_at=now())
        write_json(OUT/'run.json',record)
        subprocess.run([sys.executable,str(snap/'compare.py')],env=env,check=True)
        record.update(status='complete' if not failed else 'complete_with_failures',stage='complete',finished_at=now())
        write_json(OUT/'run.json',record)
        if failed:raise RuntimeError(f'Failed tools require follow-up: {failed}; others continued')
        print(f'Complete: {OUT}',flush=True)

if __name__=='__main__':
    try:main()
    except BaseException as error:
        path = OUT/'run.json'
        if path.exists():
            record=json.loads(path.read_text())
            record.update(status='failed',error=f'{type(error).__name__}: {error}',updated_at=now())
            write_json(path,record)
        raise
