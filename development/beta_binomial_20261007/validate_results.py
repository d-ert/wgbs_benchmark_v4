"""Independent checks of exported BB results, without calling the scorer.

Recompute CpG membership from truth/eligibility, BH from its rank formula,
AP from precision at tied thresholds, and region match cardinality with
SciPy's bipartite matching rather than the scorer's NetworkX matching.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT/'results/beta_binomial_20261007_full'
sys.path[:0] = [str(ROOT/'vendor310')]
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching


def validate_dataset(name):
    dataset = ROOT/'data'/name
    target = RUN/'simulated'/name/'epykit'
    recorded = json.loads((target/'score.json').read_text())
    regions = pd.read_csv(target/'dmr.tsv',sep='\t')
    truth,all_region_ids = {},set()
    eligible_count,positive_count,full_positive = 0,0,0
    edges = {cutoff:[] for cutoff in [.5,.8]}
    for part in sorted((dataset/'truth_parts').glob('*.parquet')):
        chrom = part.stem
        t = pq.read_table(part,columns=['pos','is_dml','region_id'])
        e = pq.read_table(dataset/'eligibility_parts'/part.name,columns=['pos','eligible'])
        positions = t['pos'].to_numpy()
        np.testing.assert_array_equal(positions,e['pos'].to_numpy())
        assert np.all(np.diff(positions)>0)
        keep = e['eligible'].to_numpy().astype(bool)
        site_positive = t['is_dml'].to_numpy().astype(bool)
        ids = t['region_id'].to_numpy()
        all_region_ids.update(int(v) for v in np.unique(ids) if v>0)
        full_positive += int(site_positive.sum())
        pos,y,rid = positions[keep],site_positive[keep],ids[keep]
        eligible_count += len(pos)
        positive_count += int(y.sum())
        truth[chrom] = (pos,y,rid)
        sizes = {int(k):int(v) for k,v in zip(*np.unique(rid[rid>0],return_counts=True))}
        for index,region in regions[regions.chrom==chrom].iterrows():
            left,right = np.searchsorted(pos,region.start),np.searchsorted(pos,region.end,side='right')
            called_size = right-left
            if called_size==0:continue
            overlaps = rid[left:right]
            for tid,overlap in zip(*np.unique(overlaps[overlaps>0],return_counts=True)):
                for cutoff in edges:
                    if overlap/called_size>=cutoff and overlap/sizes[int(tid)]>=cutoff:
                        edges[cutoff].append((int(index),int(tid)))
    assert set(regions.chrom)<=set(truth)
    p_chunks,q_chunks,y_chunks = [],[],[]
    seen = {chrom:0 for chrom in truth}
    last = {chrom:0 for chrom in truth}
    tp,fp = 0,0
    for chunk in pd.read_csv(target/'dml.tsv',sep='\t',usecols=['chrom','pos','pvalue','qvalue'],chunksize=500000):
        labels = np.empty(len(chunk),dtype=bool)
        for chrom,rows in chunk.groupby('chrom',sort=False):
            assert chrom in truth
            pos,y,_ = truth[chrom]
            observed = rows.pos.to_numpy()
            assert observed[0]>last[chrom] and np.all(np.diff(observed)>0)
            idx = np.searchsorted(pos,observed)
            assert np.all(idx<len(pos))
            np.testing.assert_array_equal(pos[idx],observed)
            labels[rows.index.to_numpy()-chunk.index[0]] = y[idx]
            seen[chrom] += len(rows)
            last[chrom] = int(observed[-1])
        p,q = chunk.pvalue.to_numpy(),chunk.qvalue.to_numpy()
        found = np.nan_to_num(q,nan=1.)<=.05
        tp += int((found&labels).sum());fp += int((found&~labels).sum())
        p_chunks.append(p.copy());q_chunks.append(q.copy());y_chunks.append(labels)
    assert all(seen[chrom]==len(truth[chrom][0]) for chrom in truth), 'Incomplete eligible CpG export'
    p,q,y = np.concatenate(p_chunks),np.concatenate(q_chunks),np.concatenate(y_chunks)
    counts = dict(tp=tp,fp=fp,fn=positive_count-tp,tn=eligible_count-positive_count-fp,
                  n_eligible=eligible_count,n_truth=full_positive,
                  n_scored=int((np.nan_to_num(p,nan=1.)<1).sum()))
    for key,value in counts.items():assert value==recorded['dml'][key],(name,key,value,recorded['dml'][key])
    assert np.all((p[np.isfinite(p)]>=0)&(p[np.isfinite(p)]<=1))
    finite = np.flatnonzero(np.isfinite(p))
    order = np.argsort(p[finite])
    adjusted = len(finite)*p[finite[order]]/np.arange(1,len(finite)+1)
    adjusted = np.minimum(1,np.minimum.accumulate(adjusted[::-1])[::-1])
    expected_q = np.full(len(p),np.nan)
    expected_q[finite[order]] = adjusted
    np.testing.assert_allclose(q,expected_q,atol=2e-12,rtol=2e-12,equal_nan=True)
    safe = np.nan_to_num(p,nan=1.)
    rank = np.argsort(safe,kind='stable')
    sorted_p = safe[rank]
    ap = None
    if sorted_p[0]!=sorted_p[-1] and y.any() and (~y).any():
        endpoints = np.r_[np.flatnonzero(sorted_p[:-1]!=sorted_p[1:]),len(p)-1]
        positives = np.cumsum(y[rank],dtype=np.int64)[endpoints]
        ap = float(np.sum((positives/(endpoints+1))*np.diff(np.r_[0,positives])/y.sum()))
        np.testing.assert_allclose(ap,recorded['dml']['average_precision'],atol=1e-12,rtol=1e-12)
    else:assert recorded['dml']['average_precision'] is None
    truth_ids = sorted(all_region_ids)
    mapping = {tid:i for i,tid in enumerate(truth_ids)}
    region_checks = {}
    for cutoff,pairs in edges.items():
        matched = 0
        if len(regions) and truth_ids and pairs:
            a,b = zip(*pairs)
            matrix = csr_matrix((np.ones(len(pairs),dtype=np.int8),(a,[mapping[tid] for tid in b])),
                                shape=(len(regions),len(truth_ids)))
            matched = int((maximum_bipartite_matching(matrix,perm_type='column')>=0).sum())
        key = 'dmr' if cutoff==.5 else 'dmr_80'
        assert matched==recorded[key]['matched_regions']
        assert len(truth_ids)==recorded[key]['true_regions']
        assert len(regions)==recorded[key]['called_regions']
        region_checks[key] = dict(called=len(regions),true=len(truth_ids),matched=matched)
    return dict(dataset=name,cpg=counts,average_precision=ap,regions=region_checks,
                bh_max_absolute_difference=float(np.nanmax(np.abs(q-expected_q))),
                scorer_functions_called=False)


def main():
    assert json.loads((RUN/'run.json').read_text())['status']=='complete'
    results = [validate_dataset(name) for name in ['genome_autosomes_signal','genome_autosomes_null']]
    value = dict(status='verified',datasets=results,
                 validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (RUN/'comparison/independent_validation.json').write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(value,indent=2),flush=True)


if __name__=='__main__':main()
