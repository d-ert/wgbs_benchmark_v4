"""Evaluate frozen updates on the existing full benchmark or real cohort.

No package tests or newly generated simulations. Counts, truth and matching
rules are the recorded workflow inputs. Fresh output directories are required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import pandas as pd
import polars as pl


def write_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=str)+'\n')
    temporary.replace(path)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('dataset',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--engine',choices=['legacy','optimized_bb','score_shared','score_group'],required=True)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--paired',action='store_true')
    args=parser.parse_args()
    out=args.output
    out.mkdir(parents=True,exist_ok=True)
    if (out/'store').exists():
        raise ValueError('Preserve existing benchmark output; use a fresh destination')
    os.environ.setdefault('NUMBA_CACHE_DIR',str(out/'numba_cache'))
    import epykit as ep
    if not Path(ep.__file__).resolve().is_relative_to(args.source.resolve()):
        raise RuntimeError(f'Incorrect source imported: {ep.__file__}')
    source_config=json.loads(args.config.read_text())
    dmc=dict(source_config['dmc'])
    dmr=dict(source_config['dmr'])
    dmc.update(materialize=False,tsv=False)
    dmr.update(tsv=False)
    if args.engine=='optimized_bb':
        dmc.update(test='beta_binomial',dispersion='eb',reference='adaptive',
                   smoothing=False,beta_binomial_ci=False)
    if args.engine.startswith('score_'):
        dmc.update(test='bb_score',dispersion='eb',smoothing=False,
                   score_dispersion=args.engine.removeprefix('score_'),score_ci=False,
                   score_prior_max_groups=4096,score_prior_seed=20261009,
                   score_reference='normal',score_block_size=32768)
        dmc.pop('beta_binomial_ci',None)
    if args.paired:
        if args.engine=='optimized_bb':
            raise ValueError('The BB likelihood engine does not support donor designs')
        dmc.update(formula='~ donor + treatment',contrast='treatment')
        if args.engine=='legacy':
            dmc.update(test='glm',materialize=True)
    meta=json.loads((args.dataset/'manifest.json').read_text())
    sheet=pd.read_csv(args.dataset/'samples.tsv',sep='\t')
    if set(sheet.group)!={'treatment','control'}:
        raise ValueError('Require the original two groups')
    if args.paired:
        if 'donor' not in sheet or not sheet.groupby('donor').group.apply(lambda x:set(x)=={'treatment','control'} and len(x)==2).all():
            raise ValueError('Require complete donor pairs from the original sample sheet')
    strand_collapsed=bool(meta.get('strand_collapsed')) or meta.get('study')=='GSE64177'
    if strand_collapsed:
        import epykit.convert as conversion
        original=conversion.convert_sample
        def convert_collapsed(*a,**kw):
            kw['merge_strands']=False
            return original(*a,**kw)
        conversion.convert_sample=convert_collapsed
    config=dict(engine=args.engine,paired=args.paired,dmc=dmc,dmr=dmr,
        original_config=str(args.config),original_config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
        dataset=str(args.dataset),dataset_manifest_sha256=hashlib.sha256((args.dataset/'manifest.json').read_bytes()).hexdigest(),
        samplesheet_sha256=hashlib.sha256((args.dataset/'samples.tsv').read_bytes()).hexdigest(),
        actual_epykit_module=str(ep.__file__),assembly=meta.get('assembly','hg38'),
        input_strand_collapsed=strand_collapsed,coverage_contract='original per-sample counts; intersect present sites; no high-depth trimming',
        region_coordinates='1-based inclusive; native exclusive end minus one',
        comparisons='score engines use raw counts even when the historical paper profile smooths counts; all other candidate geometry/cutoffs retained')
    write_json(out/'runner_config.json',config)
    phases={}
    begin=time.monotonic()
    sample_csv=out/'samples.csv'
    sheet.to_csv(sample_csv,index=False)
    print('READING original counts',flush=True)
    md=ep.read_bismark(str(sample_csv),treatment_group='treatment',control_group='control',
        assembly=meta.get('assembly','hg38'),store_dir=str(out/'store'))
    ep.pp.set_unite_type(md,type='intersect')
    loaded=time.monotonic()
    phases['read_filter_s']=loaded-begin
    write_json(out/'progress.json',dict(stage='dmc',phases=phases,started_unix=time.time()))
    print('FITTING',args.engine,'paired',args.paired,flush=True)
    ep.tl.dmc(md,**dmc)
    if args.engine.startswith('score_') and md.uns['dmc']['test_used']!='bb_score':
        raise RuntimeError('Requested count engine dispatched to another model')
    write_json(out/'dmc_metadata.json',md.uns['dmc'])
    columns=['chrom','pos','pvalue','qvalue','meth_diff']
    handle=md.dmc_store
    rows,first=0,True
    with (out/'dml.tsv').open('w') as stream:
        for chrom,frame in handle.iter_chroms(columns=columns):
            frame.write_csv(stream,separator='\t',include_header=first)
            rows+=len(frame)
            first=False
            print('EXPORTED DMC',chrom,len(frame),flush=True)
    if first:
        pl.DataFrame(schema={c:pl.String if c=='chrom' else pl.Float64 for c in columns}).write_csv(out/'dml.tsv',separator='\t')
    fitted=time.monotonic()
    phases['dml_s']=fitted-loaded
    write_json(out/'progress.json',dict(stage='dmr',phases=phases,dmc_rows=rows))
    print('CALLING chain_merge',flush=True)
    ep.tl.dmr(md,**dmr)
    regions=md.uns['dmr']
    if isinstance(regions,dict):
        regions=regions.get('frame',next(iter(regions.values())))
    regions=regions.with_columns((pl.col('end')-1).alias('end'))
    regions.write_csv(out/'dmr_candidates.tsv',separator='\t')
    qcol='combined_qvalue' if 'combined_qvalue' in regions.columns else 'qvalue'
    selected=regions.filter(pl.col(qcol).is_not_null()&(pl.col(qcol)<=.05))
    selected.write_csv(out/'dmr.tsv',separator='\t')
    phases['dmr_s']=time.monotonic()-fitted
    write_json(out/'timing.json',phases)
    write_json(out/'native_metadata.json',dict(dmc=md.uns['dmc'],dmr=md.uns['dmr_params'],epykit_version=ep.__version__))
    write_json(out/'progress.json',dict(stage='complete',phases=phases,dmc_rows=rows,
        candidates_exported=len(regions),significant_regions=len(selected)))
    print('COMPLETE',rows,'CpGs',len(selected),'regions',flush=True)


if __name__=='__main__':
    main()
