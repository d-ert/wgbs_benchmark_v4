"""Independent count-level calibration probes for default epykit CpG inference.

Ground truth is used only for evaluation. Candidate inference receives counts.
This probe does not read the original benchmark truth or yesterday's work.
"""
from pathlib import Path
import json
import logging
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
SOURCE = Path('/scratch/epykit-calibration-fresh_20261007/src')
OUT = Path(__file__).parent
sys.path[:0] = [str(ROOT/'vendor310'), str(SOURCE)]
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import average_precision_score
from epykit.dmc import _score_finalize
logging.getLogger('epykit').setLevel(logging.ERROR)


def generate(seed, n_sites, n_rep, depth, unbalanced, signal):
    rng = np.random.default_rng(seed)
    means = rng.choice([.01,.05,.15,.35,.5,.65,.85,.95,.99], n_sites)
    rho = rng.choice([.001,.01,.05,.15,.35], n_sites)
    truth = rng.random(n_sites) < .01 if signal else np.zeros(n_sites, dtype=bool)
    case_mean = np.clip(means + np.where(means<=.5,.2,-.2)*truth, .001,.999)
    factors = np.geomspace(.2,2.,n_rep) if unbalanced else np.ones(n_rep)
    groups = []
    for mu, scale in [(case_mean,factors), (means,factors[::-1])]:
        n = np.maximum(1, rng.poisson(depth*scale[None,:], (n_sites,n_rep))).astype(float)
        concentration = (1-rho)/rho
        probability = rng.beta((mu*concentration)[:,None], ((1-mu)*concentration)[:,None], size=n.shape)
        m = rng.binomial(n.astype(int), probability).astype(float)
        groups.append((m,n))
    return groups, truth, means, rho


def accum(group):
    m,n = group
    return n.sum(1),m.sum(1),(m*m/n).sum(1),np.full(len(n),n.shape[1],dtype=np.int32)


def main():
    rows=[];strata=[]
    for n_rep,depth,unbalanced in [(3,20,False),(5,20,False),(5,20,True),(5,60,True),(10,20,True)]:
        for signal in [False,True]:
            groups,truth,mean,rho=generate(731003+n_rep+depth+int(signal),100000,n_rep,depth,unbalanced,signal)
            args=(*accum(groups[0]),*accum(groups[1]))
            for mode in ['eb','site','bb_eb']:
                output=_score_finalize(*args,dispersion=mode,reference='adaptive',
                    sn2_case=(groups[0][1]**2).sum(1),sn2_ctrl=(groups[1][1]**2).sum(1))
                p=output[0]; valid=np.isfinite(p);safe=np.where(valid,p,1.)
                q=multipletests(safe,method='fdr_bh')[1];found=q<=.05
                tp=int((found&truth).sum());fp=int((found&~truth).sum())
                rows.append(dict(n_rep=n_rep,depth=depth,unbalanced=unbalanced,signal=signal,method=mode,
                    raw_fpr05=float((safe[~truth]<.05).mean()),raw_fpr001=float((safe[~truth]<.001).mean()),
                    tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0.,recall=tp/truth.sum() if truth.sum() else None,
                    average_precision=average_precision_score(truth,-safe) if truth.any() else None,
                    median_phi=float(np.median(output[5]))))
                if not signal and mode=='eb':
                    for m in np.unique(mean):
                        for r in np.unique(rho):
                            select=(mean==m)&(rho==r)
                            strata.append(dict(n_rep=n_rep,depth=depth,unbalanced=unbalanced,mean=m,rho=r,
                                sites=int(select.sum()),raw_fpr05=float((safe[select]<.05).mean()),
                                raw_fpr001=float((safe[select]<.001).mean()),bh_calls=int(found[select].sum())))
    pd.DataFrame(rows).to_csv(OUT/'baseline_probe.csv',index=False)
    pd.DataFrame(strata).to_csv(OUT/'baseline_null_strata.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    (OUT/'probe_context.json').write_text(json.dumps(dict(seed_base=731003,n_sites=100000,
        models='Independent beta-binomial counts; heterogeneous means and rho; Poisson depths; 1% signal with delta .2',
        source=str(SOURCE),purpose='Development diagnosis only; separate evaluation seeds required'),indent=2)+'\n')


if __name__=='__main__':main()
