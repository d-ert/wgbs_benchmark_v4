"""Run frozen six-tool simulation/real-data screens and collect reviewable outputs.

The real dataset has no truth. Its outputs are native DMR lists and descriptive
cross-tool concordance, never simulated precision/FDR. Every comparative row
uses the existing runner defaults; exploratory epykit runs are kept separate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from wgbs_v3.benchmark import TOOLS,benchmark,run_monitored
from wgbs_v3.provenance import write_json

CONFIG_PREFIXES=('WGBS_EPYKIT_','WGBS_DSS_','WGBS_BSMOOTH_','WGBS_DMRCATE_','WGBS_METHYLKIT_')
DONOR_RE=re.compile(r'_DC(\d+)_(MTB|NI)_5mC\.cov$')


def environment(ep_source:Path) -> dict:
    leaked=[k for k in os.environ if k.startswith(CONFIG_PREFIXES)]
    if leaked:raise ValueError(f'Remove runner overrides before default comparison: {leaked}')
    env=dict(os.environ)
    env['PYTHONPATH']=':'.join(map(str,(ROOT/'vendor310',ROOT/'src',ep_source)))
    env['MPLCONFIGDIR']=str(ROOT/'data/mplcache')
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','POLARS_MAX_THREADS','NUMBA_NUM_THREADS'):
        env[k]='1'
    return env


def real_samples(cov_dir:Path):
    samples=[]
    for p in sorted(cov_dir.glob('*.cov')):
        match=DONOR_RE.search(p.name)
        if not match:raise ValueError(f'Unexpected real sample name: {p.name}')
        donor,condition=match.groups()
        samples.append(dict(sample_id=f'DC{donor}_{condition}',group='treatment' if condition=='MTB' else 'control',
                            donor=f'DC{donor}',path=p.resolve()))
    if len(samples)!=12 or len({s['donor'] for s in samples})!=6:
        raise ValueError('Expected six NI/MTB donor pairs, twelve coverage files')
    for donor in {s['donor'] for s in samples}:
        if {s['group'] for s in samples if s['donor']==donor}!={'control','treatment'}:
            raise ValueError(f'Incomplete pair: {donor}')
    return samples


def prepare_real(cov_dir:Path,dataset:Path,max_sites:int|None=None):
    """Normalize decimal-comma percent for both methylKit and epykit readers."""
    samples=real_samples(cov_dir)
    if (dataset/'manifest.json').exists():
        prior=json.loads((dataset/'manifest.json').read_text())
        if prior.get('source_dir')!=str(cov_dir.resolve()) or prior.get('max_sites')!=max_sites:
            raise ValueError('Existing real input has different source or site limit')
        if not all(Path(s['path']).exists() for s in pd.read_csv(dataset/'samples.tsv',sep='\t').to_dict('records')):
            raise ValueError('Prepared sample is missing; use a new output directory')
        return dataset
    dataset.mkdir(parents=True,exist_ok=True)
    unexpected={p.name for p in dataset.iterdir()}-{'counts','samples.tsv','samples.tsv.partial'}
    if unexpected:raise ValueError(f'Unexpected files in real input directory: {sorted(unexpected)}')
    count_dir=dataset/'counts';count_dir.mkdir(exist_ok=True)
    expected={sample['path'].name for sample in samples}
    unexpected={p.name for p in count_dir.iterdir()}-(expected|{n+'.partial' for n in expected})
    if unexpected:raise ValueError(f'Unexpected real count files: {sorted(unexpected)}')
    output=[];audit=[]
    for sample in samples:
        source=sample['path'];dest=count_dir/source.name
        temporary=dest.with_name(dest.name+'.partial')
        hasher=hashlib.sha256();source_hash=hashlib.sha256();rows=0
        with source.open('rb') as read,temporary.open('wb') as write:
            for line in read:
                if max_sites is not None and rows>=max_sites:break
                source_hash.update(line)
                fields=line.rstrip(b'\r\n').split(b'\t')
                if len(fields)!=6 or fields[1]!=fields[2]:raise ValueError(f'Bad coverage row in {source} near {rows+1}')
                if rows<1000:
                    int(fields[1]);int(fields[4]);int(fields[5])
                fields[3]=fields[3].replace(b',',b'.')
                normalized=b'\t'.join(fields)+b'\n'
                write.write(normalized);hasher.update(normalized);rows+=1
        if not rows:raise ValueError(f'Empty coverage file: {source}')
        temporary.replace(dest)
        print(f'Prepared real sample {sample["sample_id"]}: {rows:,} sites',flush=True)
        output.append({**sample,'path':str(dest.resolve())})
        audit.append(dict(source=str(source),normalized=str(dest),rows=rows,
                          source_sha256=source_hash.hexdigest(),normalized_sha256=hasher.hexdigest(),
                          source_hash_scope='first N rows only' if max_sites is not None else 'complete file'))
    sample_sheet=dataset/'samples.tsv'
    pd.DataFrame(output).to_csv(sample_sheet.with_suffix('.tsv.partial'),sep='\t',index=False)
    sample_sheet.with_suffix('.tsv.partial').replace(sample_sheet)
    write_json(dataset/'manifest.json',dict(kind='real',assembly='hg19',study='GSE64177',
        source_dir=str(cov_dir.resolve()),max_sites=max_sites,samples=audit,
        coordinate_system='1-based inclusive',paired_donors=6,
        percent_normalization='decimal comma converted to decimal point; counts unchanged'))
    return dataset


def flatten(record,prefix=''):
    result={}
    for k,v in record.items():
        key=f'{prefix}{k}'
        if isinstance(v,dict):result.update(flatten(v,key+'_'))
        elif not isinstance(v,(list,tuple)):result[key]=v
    return result


def simulated(datasets:list[Path],out:Path,tools:tuple[str,...],ep_source:Path):
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for dataset in datasets:
        target=out/dataset.name
        manifest_path=target/'manifest.json'
        if manifest_path.exists():
            result=json.loads(manifest_path.read_text())
            if Path(result['dataset']).resolve()!=dataset.resolve():raise ValueError('Existing result uses different dataset')
            if not set(tools).issubset(result['tools']):raise ValueError('Existing result lacks requested tools; use new output')
        else:
            result=benchmark(dataset,target,tools,allow_provisional=False,ep_source=ep_source)
        for tool in tools:
            rec=result['tools'][tool]
            metrics=flatten(rec.get('score',{}))
            metrics.pop('dataset',None)
            row=dict(dataset=dataset.name,tool=tool,status=rec['status'],profile='default')
            row.update(metrics)
            row.update(wall_seconds=rec['monitoring']['wall_seconds'],
                       peak_rss_bytes=rec['monitoring']['peak_process_tree_rss_bytes'])
            rows.append(row)
    table=pd.DataFrame(rows)
    table.to_csv(out/'simulation_metrics.csv',index=False)
    if not table.status.eq('ok').all():raise RuntimeError('Some simulated runs failed; see manifests and logs')
    return table


def run_real(dataset:Path,out:Path,tools:tuple[str,...],ep_source:Path,extra_profiles:list[str]):
    out.mkdir(parents=True,exist_ok=True)
    env=environment(ep_source)
    manifest_path=out/'run_manifest.json'
    records=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    provenance_path=out/'source_provenance.json'
    source_provenance=dict(real_input_manifest_sha256=hashlib.sha256((dataset/'manifest.json').read_bytes()).hexdigest(),
        runners={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                 (ROOT/'run_epykit.py',ROOT/'run_r_tool.R',ROOT/'analysis/annotate_real_dmrs.R',Path(__file__))},
        epykit_source=str(ep_source.resolve()),assembly='hg19',tools=list(tools))
    if provenance_path.exists() and json.loads(provenance_path.read_text())!=source_provenance:
        raise ValueError('Real runner/input sources changed; use a new output directory')
    write_json(provenance_path,source_provenance)
    for tool in tools:
        target=out/tool
        if tool in records and records[tool]['status']=='ok' and (target/'dmr.tsv').exists():continue
        if target.exists():raise ValueError(f'Partial real run at {target}; inspect it and use a fresh output')
        target.mkdir()
        cmd=[sys.executable,str(ROOT/'run_epykit.py'),str(dataset),str(target)] if tool=='epykit' else [
            'Rscript',str(ROOT/'run_r_tool.R'),tool,str(dataset),str(target)]
        record=run_monitored(cmd,target/'run.log',env)
        records[tool]=dict(status='ok' if record['exit_code']==0 and (target/'dmr.tsv').exists() else 'failed',
                           monitoring=record,profile='default',assembly='hg19')
        write_json(manifest_path,records)
        if records[tool]['status']!='ok':raise RuntimeError(f'{tool} failed; inspect {target}/run.log')
    if len(tools) >= 2:
        subprocess.run(['Rscript',str(ROOT/'analysis/annotate_real_dmrs.R'),str(out),str(out/'concordance'),
                        ','.join(tools)],env=env,check=True)
