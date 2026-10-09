import numpy as np
import pandas as pd
import networkx as nx

REGION_COLUMNS=['chrom','start','end','pvalue','qvalue','direction','n_cpgs','mean_meth_diff']

def clump_sites(sites,alpha=1e-5,gap=100,minlen=50,min_cpgs=3,pct_sig=.5):
    """Raw-p CpG aggregation for methylKit; no invented region p/q-value."""
    rows=[]
    for chrom,d in sites.groupby('chrom',sort=False):
        d=d.sort_values('pos');pos=d.pos.to_numpy();p=d.pvalue.to_numpy()
        significant=np.isfinite(p)&(p<=alpha)
        seed=pos[significant]
        for group in np.split(seed,np.flatnonzero(np.diff(seed)>gap)+1):
            if not len(group):continue
            start,end=int(group[0]),int(group[-1])
            lo,hi=np.searchsorted(pos,start),np.searchsorted(pos,end,side='right')
            if end-start+1<=minlen or hi-lo<=min_cpgs or significant[lo:hi].mean()<pct_sig:continue
            effect=float(d.meth_diff.iloc[lo:hi].mean())
            rows.append(dict(chrom=chrom,start=start,end=end,pvalue=np.nan,qvalue=np.nan,
                direction=int(np.sign(effect)),n_cpgs=hi-lo,mean_meth_diff=effect))
    return pd.DataFrame(rows,columns=REGION_COLUMNS)

def clump_export(path,out):
    buffer=[];current=None;frames=[]
    for chunk in pd.read_csv(path,sep='\t',chunksize=500000):
        for chrom,d in chunk.groupby('chrom',sort=False):
            if current is not None and chrom!=current:
                frames.append(clump_sites(pd.concat(buffer,ignore_index=True)));buffer=[]
            current=chrom;buffer.append(d)
    if buffer:frames.append(clump_sites(pd.concat(buffer,ignore_index=True)))
    result=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=REGION_COLUMNS)
    result.to_csv(out,sep='\t',index=False)
    return len(result)

def paper_metrics(calls,paper):
    pairs=[]
    if 'mean_meth_diff' in calls:
        direction=np.sign(calls.mean_meth_diff.to_numpy())
    else:direction=calls.direction.to_numpy()
    for i,c in calls.reset_index(drop=True).iterrows():
        for j,p in paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)].iterrows():
            overlap=min(c.end,p.end)-max(c.start,p.start)+1
            same=direction[i]==(1 if p.direction=='hyper' else -1)
            pairs.append(dict(call_index=i,paper_index=int(j),same_direction=same,
                fraction_call=overlap/(c.end-c.start+1),fraction_paper=overlap/(p.end-p.start+1),overlap_bp=overlap))
    pairs=pd.DataFrame(pairs,columns=['call_index','paper_index','same_direction','fraction_call','fraction_paper','overlap_bp'])
    results=[]
    for cutoff in [0.,.5,.8]:
        subset=pairs[(pairs.fraction_call>=cutoff)&(pairs.fraction_paper>=cutoff)]
        same=subset[subset.same_direction.astype(bool)]
        graph=nx.Graph()
        graph.add_edges_from((('c',int(x.call_index)),('p',int(x.paper_index))) for x in same.itertuples())
        matched=len(nx.max_weight_matching(graph,maxcardinality=True))
        covered=same.paper_index.nunique();found=same.call_index.nunique()
        results.append(dict(reciprocal_threshold=cutoff,calls=len(calls),paper_total=len(paper),
            calls_overlapping_paper=found,paper_regions_covered=covered,one_to_one_matches=matched,
            coverage_precision=found/len(calls) if len(calls) else None,coverage_recall=covered/len(paper),
            one_to_one_precision=matched/len(calls) if len(calls) else None,one_to_one_recall=matched/len(paper),
            discordant_pairs=int((~subset.same_direction.astype(bool)).sum())))
    return results,pairs
