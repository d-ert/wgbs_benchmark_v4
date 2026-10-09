"""Exploratory thresholds on existing epykit calls; never changes benchmark outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import networkx as nx

ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
RUN=ROOT/'results/genome_baseline'; dataset=ROOT/'data/genome_autosomes_signal'
calls=pd.read_csv(RUN/'simulated/genome_autosomes_signal/epykit/dmr.tsv',sep='\t')
meta=json.loads((dataset/'manifest.json').read_text())
edges=[]; positive_calls=set()
for chrom, group in calls.groupby('chrom'):
    mask=pq.read_table(dataset/'eligibility_parts'/f'{chrom}.parquet').to_pandas()
    eligible=mask.loc[mask.eligible,'pos'].to_numpy()
    truth=pq.read_table(dataset/'truth_parts'/f'{chrom}.parquet',columns=['pos','region_id'],filters=[('region_id','>',0)]).to_pandas()
    truth=truth[truth.pos.isin(eligible)]
    pos=truth.pos.to_numpy(); ids=truth.region_id.to_numpy(); sizes=truth.region_id.value_counts()
    for j,row in group.iterrows():
        n=np.searchsorted(eligible,row.end,side='right')-np.searchsorted(eligible,row.start)
        if not n:continue
        lo=np.searchsorted(pos,row.start);hi=np.searchsorted(pos,row.end,side='right')
        tids,counts=np.unique(ids[lo:hi],return_counts=True)
        if len(tids):positive_calls.add(j)
        for tid,inter in zip(tids,counts):
            if inter/n>=.5 and inter/sizes[tid]>=.5:
                edges.append((('t',int(tid)),('p',int(j)),float(inter/(n+sizes[tid]-inter))))
null=pd.read_csv(RUN/'simulated/genome_autosomes_null/epykit/dmr.tsv',sep='\t')
rows=[]
for cutoff in [.05,.01,.001,.0001,.00001]:
    selected=set(calls.index[calls.combined_qvalue<=cutoff])
    graph=nx.Graph();graph.add_weighted_edges_from((t,p,weight) for t,p,weight in edges if p[1] in selected)
    matched=len(nx.max_weight_matching(graph,maxcardinality=True))
    n=len(selected)
    rows.append(dict(q_cutoff=cutoff,signal_calls=n,matched_regions=matched,precision=matched/n,
        recall=matched/meta['n_regions'],f1=2*matched/(n+meta['n_regions']),
        signal_calls_without_positive_eligible_cpg=len(selected-positive_calls),
        null_calls=int((null.combined_qvalue<=cutoff).sum())))
assert rows[0]['matched_regions']==1202 and rows[0]['signal_calls_without_positive_eligible_cpg']==2137
pd.DataFrame(rows).to_csv(OUT/'epykit_exploratory_threshold_tradeoff.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
