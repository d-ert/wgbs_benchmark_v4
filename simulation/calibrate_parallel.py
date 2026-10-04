"""Calibrate disjoint chromosome groups, then combine exact sufficient statistics."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import numpy as np
from .calibration import calibrate
from wgbs_v3.provenance import write_json


def _worker(args):
    matrix,out,samples,source,chroms,fasta=args
    if (out/'manifest.json').exists() and (out/'sufficient_statistics.npz').exists():
        return json.loads((out/'manifest.json').read_text())
    return calibrate(matrix,out,samples,source,chunk_rows=100000,
                     chromosomes=chroms,reference_fasta=fasta)


def calibrate_genome(matrix,out,samples,source,chroms,fasta,workers=4):
    groups=[chroms[i::workers] for i in range(workers)]
    out.mkdir(parents=True,exist_ok=True)
    jobs=[(matrix,out/f'part_{i+1}',samples,source,group,fasta) for i,group in enumerate(groups) if group]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        manifests=list(pool.map(_worker,jobs))
    stats=[np.load(job[1]/'sufficient_statistics.npz') for job in jobs]
    hist=sum(s['hist'] for s in stats)
    depth=sum(s['depth_sum'] for s in stats)
    spatial_sum=sum(s['spatial_sum'] for s in stats)
    spatial_count=sum(s['spatial_count'] for s in stats)
    counts={key:sum(m['counts'][key] for m in manifests) for key in manifests[0]['counts']}
    edges=np.linspace(-4,np.log10(4),401)
    tau=[]
    for row in hist:
        k=int(np.searchsorted(np.cumsum(row),(row.sum()+1)/2)) if row.sum() else None
        tau.append(float(10**((edges[k]+edges[k+1])/2)) if k is not None else None)
    observed=np.flatnonzero([x is not None for x in tau]);empty=[]
    if not len(observed):raise ValueError('No eligible CpGs')
    for i,value in enumerate(tau):
        if value is None:
            empty.append(i+1);tau[i]=tau[observed[np.argmin(np.abs(observed-i))]]
    gap=np.array(manifests[0]['spatial_correlation']['gap_edges_bp'])
    correlations=np.divide(spatial_sum,spatial_count,out=np.zeros_like(spatial_sum),where=spatial_count>0)
    mid=np.sqrt(np.maximum(gap[:-1],1)*gap[1:])
    use=(spatial_count>=1000)&(correlations>.01)&(correlations<.99)
    length=None
    if use.sum()>=2:
        slope=np.polyfit(mid[use],np.log(correlations[use]),1,w=np.sqrt(spatial_count[use]))[0]
        if slope<0:length=float(np.clip(-1/slope,25,2000))
    result=dict(manifests[0])
    result.update(scope=dict(chromosomes=chroms,max_sites=None,kind='all_autosomes'),
        counts=counts,sample_mean_depth=(depth/counts['all']).tolist(),tau_by_mean_bin=tau,
        empty_mean_bins=empty,parquet_files=[f for m in manifests for f in m['parquet_files']],
        spatial_correlation=dict(gap_edges_bp=gap.tolist(),mean_residual_correlation=correlations.tolist(),
                                 pairs=spatial_count.tolist(),fitted_length_bp=length,validated=False),
        calibration_execution='four disjoint chromosome groups; exact sufficient statistics combined')
    write_json(out/'manifest.json',result)
    return result
