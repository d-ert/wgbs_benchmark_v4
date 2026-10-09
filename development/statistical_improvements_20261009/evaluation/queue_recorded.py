"""Bounded concurrent scheduling of authorized recorded-data evaluations."""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from wgbs_v3.benchmark import run_monitored
from wgbs_v3.provenance import write_json


def main():
    p=argparse.ArgumentParser()
    p.add_argument('run',type=Path)
    a=p.parse_args()
    root=Path.cwd()
    run=a.run.resolve()
    plan=json.loads((run/'plan.json').read_text())['jobs']
    snap=run/'source_snapshot'
    env=dict(os.environ)
    env['PYTHONPATH']=':'.join([str(root/'vendor310'),str(snap),str(root/'src')])
    env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS','POLARS_MAX_THREADS','NUMBA_NUM_THREADS']})
    env.update(MPLCONFIGDIR='/tmp/epykit-real-benchmark-mpl',PYTHONDONTWRITEBYTECODE='1')
    pending=plan[2:]
    records={}
    live={}
    external={x['key']:x for x in plan[:2]}

    def evaluate(job):
        out=Path(job['output'])
        out.mkdir(parents=True,exist_ok=True)
        if job['legacy_script']:
            cmd=[sys.executable,str(snap/'run_epykit_legacy.py'),job['dataset'],str(out)]
        else:
            cmd=[sys.executable,str(snap/'run_recorded.py'),job['dataset'],str(out),
                '--engine',job['engine'],'--config',job['config'],'--source',str(snap)]
            if job['paired']:
                cmd.append('--paired')
        write_json(out/'job.json',dict(status='running',command=cmd,job=job,concurrent=True))
        print('START',job['key'],flush=True)
        monitor=run_monitored(cmd,out/'run.log',env,memory_gib=float('inf'),timeout_seconds=float('inf'))
        monitor.update(memory_limit_gib=None,timeout_seconds=None)
        record=dict(status='ok' if monitor['exit_code']==0 and (out/'dmr.tsv').is_file() else 'failed',
                    monitoring=monitor,job=job,concurrent=True)
        write_json(out/'job.json',record)
        if record['status']=='ok':
            record['scoring']=score_job(job)
            if record['scoring']['exit_code']:
                record['status']='scoring_failed'
            write_json(out/'job.json',record)
        print('END',job['key'],record['status'],flush=True)
        return record

    def score_job(job):
        out=Path(job['output'])
        cmd=[sys.executable,str(snap/'score_recorded.py'),job['kind'],job['dataset'],str(out)]
        if job['kind']=='real':
            paper=(run/'input_references/GSE64177_paper.tsv' if job['dataset'].endswith('/data/real')
                   else root/'results/GSE263850_paper_20261009/input/paper_dmrs.tsv')
            cmd.extend(['--paper',str(paper)])
        with (out/'scoring.log').open('w') as f:
            result=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
        return dict(exit_code=result.returncode,command=cmd)

    with ThreadPoolExecutor(max_workers=2) as pool:
        while external or pending or live:
            for key,job in list(external.items()):
                out=Path(job['output'])
                if key=='signal/updated_legacy':
                    entry=json.loads((run/'run.json').read_text()).get('jobs',{}).get(key)
                else:
                    entry=json.loads((out/'job.json').read_text()) if (out/'job.json').exists() else None
                    progress=out/'progress.json'
                    if (not entry or entry.get('status')=='running') and progress.exists():
                        finished=json.loads(progress.read_text())
                        if finished.get('stage')=='complete':
                            entry=dict(status='ok',monitoring=dict(exit_code=None,
                                wall_seconds=sum(finished['phases'].values()),
                                peak_process_tree_rss_bytes=None,
                                timing_quality='Recovered package phase sum; wrapper record unavailable',
                                memory_limit_gib=None,timeout_seconds=None),
                                completion_evidence='Caller wrote final complete marker and all outputs; initial wrapper cannot serialize infinite caps',
                                job=job,concurrent=True)
                if entry and entry.get('status') in ('ok','failed'):
                    print('EXTERNAL FINISHED',key,entry['status'],flush=True)
                    external.pop(key)
                    if entry['status']=='ok':
                        entry['scoring']=score_job(job)
                        if entry['scoring']['exit_code']:
                            entry['status']='scoring_failed'
                    entry['job']=job
                    records[key]=entry
                    write_json(out/'job.json',entry)
            for future,key in list(live.items()):
                if future.done():
                    try:records[key]=future.result()
                    except Exception as e:records[key]=dict(status='controller_failed',error=str(e))
                    live.pop(future)
            while pending and len(live)+len(external)<2:
                job=pending.pop(0)
                live[pool.submit(evaluate,job)]=job['key']
            state=dict(status='running' if external or live or pending else 'complete',
                controller_pid=os.getpid(),completed=records,
                active=list(external)+list(live.values()),pending=[x['key'] for x in pending],
                concurrency=2,threads_per_caller=1,timing_context='Concurrent older methylKit workflow and up to two new epykit jobs')
            write_json(run/'queue_status.json',state)
            if external or live or pending:
                time.sleep(10)
    print('QUEUE COMPLETE',len(records),flush=True)


if __name__=='__main__':
    main()
