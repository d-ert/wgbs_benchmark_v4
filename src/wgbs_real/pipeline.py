import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime,timezone

from wgbs_v3.benchmark import run_monitored
from wgbs_v3.provenance import source_inventory,write_json
from .dataset import prepare,digest
from .runners import build_r_runner,build_epykit_runner,PAPER_SETTINGS,replace_once
from .comparison import compare

ROOT=Path(__file__).resolve().parents[2]
TOOLS=['DSS','epykit','DMRcate','BSmooth','methylKit','dmrseq']
def now():return datetime.now(timezone.utc).isoformat()

def annotation_runner(dataset,target):
    text=(ROOT/'analysis/annotate_real_dmrs.R').read_text()
    text=text.replace('  library(TxDb.Hsapiens.UCSC.hg19.knownGene)\n','')
    start=text.index('txdb <- TxDb.')
    end=text.index('prom <- promoters(',start)
    block=f'''ref <- fread({json.dumps(str(Path(dataset)/'genes.tsv'))})
gene_gr <- GRanges(ref$chrom, IRanges(ref$start, ref$end), strand=ref$strand)
names(gene_gr) <- ref$gene_id
gene_id <- ref$gene_id
gene_symbol <- setNames(ref$symbol, ref$gene_id)
'''
    text=text[:start]+block+text[end:]
    text=text.replace('hg19','hg38; supplied RefGene')
    Path(target).write_text(text)

def run(raw,out,ep_source,profile='paper_aligned',tools=None,max_sites=None):
    raw,out,ep_source=Path(raw).resolve(),Path(out).resolve(),Path(ep_source).resolve()
    tools=tools or TOOLS
    if not set(tools)<=set(TOOLS) or len(set(tools))!=len(tools):raise ValueError('Unknown or duplicate tool')
    if out.exists():raise ValueError('Use a fresh run name; previous results are preserved')
    if profile not in ['paper_aligned','native']:raise ValueError('Unknown profile')
    if not (ep_source/'epykit/_region_search.py').exists():raise ValueError('Expected the region-fix epykit source')
    env=dict(os.environ)
    for key in env:
        if key.startswith(('WGBS_EPYKIT_','WGBS_DSS_','WGBS_BSMOOTH_','WGBS_DMRCATE_','WGBS_METHYLKIT_')):
            raise ValueError('Remove inherited runner overrides')
    with (ROOT/'results/.workflow.lock').open('a') as workflow,(ROOT/'results/.timing.lock').open('a') as timing:
        fcntl.flock(workflow,fcntl.LOCK_EX|fcntl.LOCK_NB);fcntl.flock(timing,fcntl.LOCK_EX|fcntl.LOCK_NB)
        out.mkdir(parents=True)
        dataset=out/'input'
        record=dict(status='running',stage='preparing',started_at=now(),study='GSE263850',assembly='hg38',
            dataset=str(dataset),source_dataset=str(raw),tools=tools,profile=profile,smoke_sites=max_sites,
            memory_limit_gib=None,timeout_seconds=None,thread_budget=1,active_tool=None,
            settings={t:PAPER_SETTINGS[t] for t in tools} if profile=='paper_aligned' else 'recorded native settings')
        write_json(out/'run.json',record)
        try:
            prepare(raw,dataset,max_sites=max_sites)
            print('Dataset prepared and verified',flush=True)
            snap=out/'source_snapshot';snap.mkdir()
            shutil.copytree(Path(__file__).parent,snap/'wgbs_real',ignore=shutil.ignore_patterns('__pycache__'))
            shutil.copytree(ep_source/'epykit',snap/'epykit',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
            build_r_runner(ROOT/'run_r_tool.R',snap/'run_r_profile.R')
            build_epykit_runner(ROOT/'run_epykit.py',snap/'run_epykit_profile.py')
            annotation_runner(dataset,snap/'annotate_hg38.R')
            record['package_sources']=source_inventory(snap/'epykit')
            record['source_inventory']={str(p.relative_to(snap)):digest(p) for p in snap.rglob('*') if p.is_file()}
            record['input_manifest_sha256']=digest(dataset/'manifest.json')
            env['PYTHONPATH']=':'.join([str(ROOT/'vendor310'),str(snap),str(ROOT/'src'),str(ROOT)])
            env['MPLCONFIGDIR']='/tmp/wgbs_gse263850_mpl'
            for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','POLARS_MAX_THREADS','NUMBA_NUM_THREADS']:env[key]='1'
            statuses={};real=out/'real';real.mkdir()
            record['stage']='benchmarking'
            for tool in tools:
                target=real/tool;target.mkdir()
                record.update(active_tool=tool,updated_at=now());write_json(out/'run.json',record)
                print(f'Starting {tool} ({profile}), no memory/time cap',flush=True)
                command=[sys.executable,'-m','wgbs_real.run_tool',tool,str(dataset),str(target),str(snap),profile]
                monitor=run_monitored(command,target/'run.log',env,memory_gib=float('inf'),timeout_seconds=float('inf'))
                monitor.update(memory_limit_gib=None,timeout_seconds=None)
                complete=monitor['exit_code']==0 and (target/'dmr.tsv').is_file()
                statuses[tool]=dict(status='ok' if complete else 'failed',monitoring=monitor,profile=profile)
                write_json(real/'run_manifest.json',statuses)
                print(f'{tool}: {statuses[tool]["status"]}',flush=True)
                compare(out,annotate=False)
            for sample in json.loads((dataset/'manifest.json').read_text())['samples']:
                assert digest(Path(sample['path']))==sample['sha256'],'Count input changed'
            assert digest(dataset/'manifest.json')==record['input_manifest_sha256']
            assert all(digest(snap/path)==value for path,value in record['source_inventory'].items())
            compare(out,annotate=True)
            failed=[t for t,s in statuses.items() if s['status']!='ok']
            record.update(status='complete' if not failed else 'complete_with_failures',stage='complete',
                active_tool=None,failed_tools=failed,finished_at=now())
            write_json(out/'run.json',record)
            if failed:raise RuntimeError(f'Failed tools need follow-up: {failed}')
        except BaseException as error:
            record.update(status='failed',error=f'{type(error).__name__}: {error}',updated_at=now())
            write_json(out/'run.json',record)
            raise
