"""Development check: fixed geometry windows before data-dependent chaining."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
ROOT=Path('/scratch/wgbs_benchmark_v4');sys.path.insert(0,str(ROOT/'src'))
from wgbs_v3.score import score
OUT=Path('/scratch/wgbs_benchmark_v4/development/epykit_calibration/results/position_screen_check');OUT.mkdir(exist_ok=False)
rows=[]
for mode in ['signal','null']:
 p=ROOT/'results/genome_baseline/simulated'/f'genome_autosomes_{mode}'/'epykit'
 f=pq.ParquetFile(p/'store/.cache/dmc/lr/chrom=chr22.parquet').read(columns=['pos','pvalue','meth_diff']).to_pandas()
 pos=f.pos.to_numpy();pv=f.pvalue.to_numpy();diff=f.meth_diff.to_numpy();n=len(pos)
 starts=np.arange(n);ends=np.maximum(starts+4,np.searchsorted(pos,pos+49))
 clipped=np.minimum(ends,n-1)
 gaps=np.r_[0,np.cumsum(np.diff(pos)>500)]
 geometry=(ends<n)&(gaps[clipped]==gaps[starts])
 valid=np.isfinite(pv)&(pv>0)&(pv<=1)&np.isfinite(diff)
 z=np.zeros(n);z[valid]=np.sign(diff[valid])*norm.isf(np.clip(pv[valid]/2,1e-300,.5))
 zs=np.r_[0,np.cumsum(z)];vs=np.r_[0,np.cumsum(valid)]
 lo=starts[geometry];hi=ends[geometry];nv=vs[hi+1]-vs[lo]
 screen_p=np.ones(len(lo));good=nv>=5
 screen_p[good]=2*norm.sf(np.abs((zs[hi[good]+1]-zs[lo[good]])/np.sqrt(nv[good])))
 q=multipletests(screen_p,method='fdr_bh')[1]
 # Original full-genome q values retained for the regional correction; the
 # independent window family here is chr22 only, matching the pilot scope.
 calls=pd.read_csv(Path('/scratch/wgbs_benchmark_v4/development/epykit_calibration/results/wgbs_calibration/chr22_pilot')/mode/'asymptotic/dmr.tsv',sep='\t')
 sq=[]
 for c in calls.itertuples():
  i=np.searchsorted(pos[lo],c.start);j=np.searchsorted(pos[lo],c.end,side='right')
  chosen=q[i:j][pos[hi[i:j]]<=c.end]
  sq.append(float(chosen.min()) if len(chosen) else 1.)
 calls['position_screen_qvalue']=sq
 calls['combined_qvalue']=np.maximum(calls.combined_qvalue,calls.position_screen_qvalue)
 calls=calls[calls.combined_qvalue<.05]
 dest=OUT/mode;dest.mkdir();calls.to_csv(dest/'dmr.tsv',sep='\t',index=False)
 view=Path('/scratch/wgbs_benchmark_v4/development/epykit_calibration/results/wgbs_calibration/chr22_pilot')/mode/'evaluation'
 m=score(view,None,dest/'dmr.tsv',dest/'score.json')['dmr']
 rows.append(dict(dataset=mode,window_tests=len(q),**m))
pd.DataFrame(rows).to_csv(OUT/'comparison.csv',index=False)
print(pd.DataFrame(rows)[['dataset','window_tests','called_regions','matched_regions','region_precision','region_recall','region_f1']].to_string(index=False))
