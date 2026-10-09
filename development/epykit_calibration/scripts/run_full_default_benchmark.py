"""Run updated default epykit on full v4 inputs and compare frozen tool results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
DEV=Path(__file__).resolve().parents[1]
REPO=Path('/scratch/epykit-calibration')
sys.path[:0]=[str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)]
import pandas as pd
from analysis.run_publication_workflow import simulated,run_real
from wgbs_v3.provenance import source_inventory,write_json


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def real_paper_comparison(call_sets,out):
    source=ROOT/'data/real/paper_reported_dmrs.xlsx'
    paper=pd.read_excel(source,header=2).dropna(subset=['Start','End','Direction'])
    paper=paper.rename(columns={'Chr':'chrom','Start':'start','End':'end','Direction':'direction'})
    rows=[]
    for name,path in call_sets.items():
        calls=pd.read_csv(path,sep='\t')
        if 'mean_meth_diff' in calls:
            calls['direction']=calls.mean_meth_diff.map(lambda x:'hyper' if x>0 else 'hypo')
        else:calls['direction']=calls.direction.map({1:'hyper',-1:'hypo'})
        for cutoff in [0,.5,.8]:
            ci=set();pi=set()
            for i,c in calls.iterrows():
                p=paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)&(paper.direction==c.direction)]
                for j,r in p.iterrows():
                    overlap=min(c.end,r.end)-max(c.start,r.start)+1
                    if overlap/(c.end-c.start+1)>=cutoff and overlap/(r.end-r.start+1)>=cutoff:ci.add(i);pi.add(j)
            rows.append(dict(tool=name,reciprocal_threshold=cutoff,calls=len(calls),calls_overlapping_paper=len(ci),caller_overlap_fraction=len(ci)/len(calls) if len(calls) else None,paper_regions_covered=len(pi),paper_total=len(paper)))
    pd.DataFrame(rows).to_csv(out/'real_paper_comparison.csv',index=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=DEV/'results/full_default_position_screen_20261006')
    args=ap.parse_args();out=args.out.resolve()
    if out.exists():raise ValueError(f'Use a new output directory: {out}')
    out.mkdir(parents=True)
    before=ROOT/'results/genome_baseline'
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    frozen=out/'source_snapshot';frozen.mkdir()
    shutil.copytree(REPO/'src/epykit',frozen/'epykit',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    shutil.copy2(Path(__file__),frozen/'run_full_default_benchmark.py')
    source_hashes=source_inventory(frozen/'epykit')
    record=dict(status='running',source_commit=commit,source_hashes=source_hashes,baseline=str(before),
                configuration='Original v4 lr/eb/chain_merge settings; empirical_fdr=False; updated default position-window screen',
                real_design='Original unpaired model retained for this comparison',
                tools_reused='All six original simulation tools; only DSS is available as a completed real-data competitor')
    write_json(out/'run.json',record)
    try:
        datasets=[ROOT/'data/genome_autosomes_signal',ROOT/'data/genome_autosomes_null']
        print('Running full signal and null simulations with default epykit...',flush=True)
        simulated(datasets,out/'simulated',('epykit',),frozen)
        print('Running full real data with default epykit...',flush=True)
        run_real(ROOT/'data/real',out/'real',('epykit',),frozen,[])
        dml_checks=[]
        for dataset in ['genome_autosomes_signal','genome_autosomes_null']:
            old=before/'simulated'/dataset/'epykit'
            new=out/'simulated'/dataset/'epykit'
            dml_checks.append(dict(dataset=dataset,dml_identical=digest(old/'dml.tsv')==digest(new/'dml.tsv')))
            config=json.loads((new/'runner_config.json').read_text())
            old_config=json.loads((old/'runner_config.json').read_text())
            assert config['dmc']==old_config['dmc'] and config['dmr']==old_config['dmr']
            assert not config['dmr']['empirical_fdr']
        dml_checks.append(dict(dataset='real',dml_identical=digest(before/'real/epykit/dml.tsv')==digest(out/'real/epykit/dml.tsv')))
        assert all(c['dml_identical'] for c in dml_checks),'Unexpected CpG output changes; inspect before interpreting region results'
        write_json(out/'regression_checks.json',dml_checks)
        original=pd.read_csv(before/'simulated/simulation_metrics.csv');original['source_run']='genome_baseline';original['reused']=True
        updated=pd.read_csv(out/'simulated/simulation_metrics.csv');updated['source_run']=out.name;updated['reused']=False
        comparison=pd.concat([original,updated],ignore_index=True)
        comparison.to_csv(out/'simulation_comparison.csv',index=False)
        old=original[original.tool=='epykit'].set_index('dataset');new=updated.set_index('dataset')
        changes=[]
        for dataset in old.index:
            for field in old.select_dtypes(include='number').columns.intersection(new.select_dtypes(include='number').columns):
                a,b=old.at[dataset,field],new.at[dataset,field]
                changes.append(dict(dataset=dataset,metric=field,before=a,after=b,delta=b-a))
        pd.DataFrame(changes).to_csv(out/'epykit_changes.csv',index=False)
        real_comp=out/'real_comparison';real_comp.mkdir()
        sets={'epykit_before':before/'real/epykit/dmr.tsv','epykit_after':out/'real/epykit/dmr.tsv','DSS':before/'real/DSS/dmr.tsv'}
        real_rows=[]
        for name,path in sets.items():
            dest=real_comp/name;dest.mkdir();shutil.copy2(path,dest/'dmr.tsv')
            calls=pd.read_csv(path,sep='\t');is_new=name=='epykit_after'
            manifests=json.loads(((out if is_new else before)/'real/run_manifest.json').read_text())
            tool='DSS' if name=='DSS' else 'epykit';monitor=manifests[tool]['monitoring']
            real_rows.append(dict(tool=name,dmr_count=len(calls),wall_seconds=monitor['wall_seconds'],peak_rss_bytes=monitor['peak_process_tree_rss_bytes'],reused=not is_new))
        pd.DataFrame(real_rows).to_csv(out/'real_comparison.csv',index=False)
        subprocess.run(['Rscript',str(ROOT/'analysis/annotate_real_dmrs.R'),str(real_comp),str(real_comp/'concordance'),','.join(sets)],check=True)
        real_paper_comparison(sets,out)
        assert source_inventory(frozen/'epykit')==source_hashes
        record['status']='complete';write_json(out/'run.json',record)
        print('COMPLETE:',out,flush=True)
        print(updated[['dataset','tool','dmr_called_regions','dmr_matched_regions','dmr_region_precision','dmr_region_recall','dmr_region_f1']].to_string(index=False))
        print(pd.DataFrame(real_rows).to_string(index=False))
    except BaseException as exc:
        record.update(status='failed',error=f'{type(exc).__name__}: {exc}');write_json(out/'run.json',record)
        raise

if __name__=='__main__':main()
