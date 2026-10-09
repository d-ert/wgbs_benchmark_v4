"""Quantify region detection lost before vs after the q-value gate."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor310'))
import polars as pl
import pandas as pd
import networkx as nx
catalog={}
for p in (ROOT/'data/genome_autosomes_signal/truth_parts').glob('*.parquet'):
 e=pl.read_parquet(ROOT/'data/genome_autosomes_signal/eligibility_parts'/p.name)
 unfiltered=pl.read_parquet(p,columns=['chrom','pos','region_id'])
 assert np.array_equal(unfiltered['pos'].to_numpy(), e['pos'].to_numpy())
 t=unfiltered.filter(e['eligible'])
 chrom=str(t['chrom'][0]);pos=t['pos'].to_numpy();ids=t['region_id'].to_numpy()
 catalog[chrom]=(pos,ids,np.bincount(ids.clip(min=0)))
rows=[]
for run,tool in [('genome_baseline','original'),('default_search_by_20261007','region_fix'),('beta_binomial_20261007_full','BB_F')]:
 root=ROOT/'results'/run/'simulated/genome_autosomes_signal/epykit/store/.cache/dmc'
 caches=list(root.rglob('.dmr_chain_merge*.parquet'));assert len(caches)==1
 candidates=pl.read_parquet(caches[0]).to_pandas()
 if tool=='original':original=candidates
 if tool=='region_fix':
  keys=['chrom','start','end','n_cpgs','n_significant','mean_meth_diff','combined_pvalue','dmr_type']
  assert original[keys].equals(candidates[keys])
 for stage in ['before_q','after_q']:
  calls=candidates if stage=='before_q' else candidates[candidates.combined_qvalue<=.05]
  for threshold in [.5,.8]:
   graph=nx.Graph();supported=set();any_truth=set()
   for i,r in calls.iterrows():
    pos,ids,sizes=catalog[str(r.chrom)]
    lo=np.searchsorted(pos,r.start);hi=np.searchsorted(pos,r.end-1,side='right')
    span=ids[lo:hi];tids,counts=np.unique(span[span>0],return_counts=True)
    any_truth.update(int(x) for x in tids)
    for tid,inter in zip(tids,counts):
     if inter/len(span)>=threshold and inter/sizes[tid]>=threshold:
      graph.add_edge(('t',int(tid)),('p',int(i)),weight=float(inter/(len(span)+sizes[tid]-inter)))
      supported.add(int(tid))
   matches=nx.max_weight_matching(graph,maxcardinality=True,weight='weight')
   rows.append(dict(tool=tool,stage=stage,threshold=threshold,calls=len(calls),matched=len(matches),truth_regions_with_any_overlap=len(any_truth),truth_regions_with_qualifying_overlap=len(supported),recall=len(matches)/2893))
pd.DataFrame(rows).to_csv(OUT/'candidate_ceiling.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
