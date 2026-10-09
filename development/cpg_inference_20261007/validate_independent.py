"""Held-out count validation of the default EB uncertainty correction."""
from pathlib import Path
import json
import logging
import subprocess
import sys
import types

ROOT=Path('/scratch/wgbs_benchmark_v4')
REPO=Path('/scratch/epykit-calibration-fresh_20261007')
OUT=Path(__file__).parent
sys.path[:0]=[str(ROOT/'vendor310'),str(REPO/'src'),str(OUT)]
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import average_precision_score
from diagnose import generate,accum
import epykit.dmc as revised
logging.getLogger('epykit').setLevel(logging.ERROR)


def main():
    module=types.ModuleType('epykit._pre_cpg_fix');module.__package__='epykit'
    sys.modules[module.__name__]=module
    original=subprocess.check_output(['git','show','f73e66d:src/epykit/dmc.py'],cwd=REPO,text=True)
    exec(compile(original,module.__name__,'exec'),module.__dict__)
    rows=[];strata=[]
    scenarios=[(3,20,False),(5,20,False),(5,20,True),(5,60,True),(10,20,True)]
    seeds=[8641021,8642023,8643019]
    for seed in seeds:
        for nrep,depth,imbalance in scenarios:
            for signal in [False,True]:
                groups,truth,mean,rho=generate(seed+int(signal),100000,nrep,depth,imbalance,signal)
                args=(*accum(groups[0]),*accum(groups[1]))
                boundary=np.any([np.all(g[0]==0,1)|np.all(g[0]==g[1],1) for g in groups],0)
                for name,engine in [('region_fix_only',module),('local_group_df',revised)]:
                    result=engine._score_finalize(*args,dispersion='eb',reference='adaptive')
                    p=np.nan_to_num(result[0],nan=1.);q=multipletests(p,method='fdr_bh')[1]
                    found=q<=.05;tp=int((found&truth).sum());fp=int((found&~truth).sum())
                    rows.append(dict(seed=seed,nrep=nrep,depth=depth,imbalance=imbalance,signal=signal,
                        version=name,sites=len(truth),true_sites=int(truth.sum()),tp=tp,fp=fp,
                        fdp=fp/(tp+fp) if tp+fp else 0.,recall=tp/truth.sum() if truth.any() else None,
                        raw_null_fpr05=float((p[~truth]<.05).mean()),raw_null_fpr001=float((p[~truth]<.001).mean()),
                        average_precision=average_precision_score(truth,-p) if truth.any() else None,
                        boundary_fp=int((found&~truth&boundary).sum()),median_df=float(np.median(result[-1]))))
                    if not signal:
                        for mu in np.unique(mean):
                            for r in np.unique(rho):
                                select=(mean==mu)&(rho==r)
                                strata.append(dict(seed=seed,nrep=nrep,depth=depth,imbalance=imbalance,
                                    version=name,mean=mu,rho=r,sites=int(select.sum()),
                                    fpr05=float((p[select]<.05).mean()),fpr001=float((p[select]<.001).mean()),
                                    bh_false_calls=int(found[select].sum())))
    table=pd.DataFrame(rows);table.to_csv(OUT/'heldout_count_validation.csv',index=False)
    pd.DataFrame(strata).to_csv(OUT/'heldout_null_strata.csv',index=False)
    summary=table.groupby(['version','nrep','depth','imbalance','signal']).agg(
        runs=('seed','size'),tp_mean=('tp','mean'),fp_mean=('fp','mean'),
        fdp_mean=('fdp','mean'),recall_mean=('recall','mean'),
        raw_null_fpr05=('raw_null_fpr05','mean'),raw_null_fpr001=('raw_null_fpr001','mean'),
        average_precision=('average_precision','mean')).reset_index()
    summary.to_csv(OUT/'heldout_count_summary.csv',index=False)
    context=dict(seeds=seeds,n_sites_per_run=100000,scenario_count=len(scenarios),
        comparison='Same generated counts per version; default lr/eb/adaptive; BH q <= .05',
        model='Independent beta-binomial counts with heterogeneous means/rho, Poisson depths; signal prevalence1%, delta .2',
        scope='Held-out seeds for bounded uncertainty patch; not genome-wide validation or spatially correlated validation',
        source_before='f73e66d',source_after=str(REPO),inference_uses='observed counts only')
    (OUT/'heldout_validation_context.json').write_text(json.dumps(context,indent=2)+'\n')
    print(summary.to_string(index=False))


if __name__=='__main__':main()
