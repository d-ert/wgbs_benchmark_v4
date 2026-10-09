"""Evaluate a default-path EB change on frozen chr22 counts, without permutations."""
from pathlib import Path
import json,sys
import polars as pl
import pandas as pd
import epykit as ep
from epykit.dmr import call_dmr_chain_merge,apply_region_qfilter
ROOT=Path('/scratch/wgbs_benchmark_v4');sys.path.insert(0,str(ROOT/'src'))
from wgbs_v3.score import score
OUT=Path('/scratch/wgbs_benchmark_v4/development/epykit_calibration/results/wgbs_default_eb_check');OUT.mkdir(exist_ok=False)
rows=[]
for mode in ['signal','null']:
    prior=ROOT/'results/genome_baseline/simulated'/f'genome_autosomes_{mode}'/'epykit'
    path=OUT/mode;path.mkdir()
    sheet=pd.read_csv(ROOT/'data'/f'genome_autosomes_{mode}'/'samples.tsv',sep='\t')
    md=ep.MethylData(obs=pl.from_pandas(sheet.assign(treatment=(sheet.group=='treatment').astype(int))),
       store=str(prior/'store/.cache/raw'),assembly='hg38',analysis_root=str(path.resolve()))
    ep.pp.set_unite_type(md,type='intersect')
    ep.tl.dmc(md,test='lr',dispersion='eb',reference='adaptive',smoothing=False,chromosomes=['chr22'],tsv=False)
    p=md.dmc['pvalue'].to_numpy()
    cm=dict(alpha=.05,min_abs_meth_diff=.1,min_cpgs=5,minlen_bp=50,dis_merge_bp=500,pct_sig=.5,use_q_for_sig=False)
    frame=apply_region_qfilter(call_dmr_chain_merge(md.dmc,**cm),.05)
    frame.with_columns((pl.col('end')-1).alias('end')).write_csv(path/'dmr.tsv',separator='\t')
    md.dmc.select('chrom','pos','pvalue','qvalue','meth_diff').write_csv(path/'dml.tsv',separator='\t')
    view=Path('/scratch/wgbs_benchmark_v4/development/epykit_calibration/results/wgbs_calibration/chr22_pilot')/mode/'evaluation'
    metrics=score(view,path/'dml.tsv',path/'dmr.tsv',path/'score.json')
    rows.append(dict(dataset=mode,p_below_05=float((p<.05).mean()),p_below_001=float((p<.001).mean()),**metrics['dmr'],dml_tp=metrics['dml']['tp'],dml_fp=metrics['dml']['fp']))
pd.DataFrame(rows).to_csv(OUT/'comparison.csv',index=False)
print(pd.DataFrame(rows)[['dataset','called_regions','matched_regions','region_precision','region_recall','region_f1','p_below_05','p_below_001','dml_tp','dml_fp']].to_string(index=False))
