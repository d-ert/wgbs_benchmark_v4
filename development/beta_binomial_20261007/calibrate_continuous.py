"""Independent development simulations: inference gets counts, never truth.

This development set can guide the candidate choice. Subsequent validation
must use untouched seeds. It does not read genome-baseline or paper truth.
"""
from pathlib import Path
import json
import sys
import time

ROOT = Path('/scratch/wgbs_benchmark_v4')
SOURCE = Path('/scratch/epykit-calibration-fresh_20261007/src')
OUT = Path(__file__).parent / 'continuous_development'
OUT.mkdir(exist_ok=True)
sys.path[:0] = [str(ROOT/'vendor310'),str(SOURCE)]
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score
from statsmodels.stats.multitest import multipletests
from epykit._beta_binomial import estimate_prior, test_lognormal
from epykit.dmc import _score_finalize


def generate(seed,sites,nrep,depth,imbalance,signal,association=False):
    rng = np.random.default_rng(seed)
    mean = rng.choice([.001,.01,.05,.15,.35,.5,.65,.85,.95,.99,.999],sites)
    rho = rng.choice([.001,.01,.05,.15,.35],sites)
    if association:
        rho = np.where((mean<.05)|(mean>.95),.35,.005)
    truth = rng.random(sites)<.01 if signal else np.zeros(sites,dtype=bool)
    case_mean = np.clip(mean+np.where(mean<=.5,.2,-.2)*truth,.001,.999)
    factor = np.geomspace(.2,2,nrep) if imbalance else np.ones(nrep)
    methylated,coverage = [],[]
    for mu,scale in [(case_mean,factor),(mean,factor[::-1])]:
        n = np.maximum(1,rng.poisson(depth*scale[None,:],(sites,nrep)))
        k = (1-rho)/rho
        probability = rng.beta((mu*k)[:,None],((1-mu)*k)[:,None],size=n.shape)
        methylated.append(rng.binomial(n,probability))
        coverage.append(n)
    return np.column_stack(methylated),np.column_stack(coverage),truth,mean,rho


def accum(m,n):
    return n.sum(1).astype(float),m.sum(1).astype(float),(m.astype(float)**2/n).sum(1),np.full(len(m),m.shape[1])


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value))


def main():
    sites = int(sys.argv[1]) if len(sys.argv)>1 else 50000
    scenarios = [(3,20,False,False),(5,20,False,False),(5,20,True,False),
                 (5,60,True,False),(10,20,True,False),(5,20,True,True)]
    records,strata = [],[]
    for index,(nrep,depth,imbalance,association) in enumerate(scenarios):
        for signal in [False,True]:
            seed = 10021001+index*101+int(signal)
            label = f'n{nrep}_d{depth}_imb{imbalance}_assoc{association}_signal{signal}'
            m,n,truth,mean,rho = generate(seed,sites,nrep,depth,imbalance,signal,association)
            start = time.monotonic()
            prior = estimate_prior(m,n,nrep)
            output = test_lognormal(m,n,nrep,prior=prior,reference='F')
            seconds = time.monotonic()-start
            (OUT/f'prior_{label}.json').write_text(json.dumps(output['prior'],default=serial,indent=2)+'\n')
            old = _score_finalize(*accum(m[:,:nrep],n[:,:nrep]),*accum(m[:,nrep:],n[:,nrep:]),
                dispersion='eb',reference='adaptive',sn2_case=(n[:,:nrep].astype(float)**2).sum(1),
                sn2_ctrl=(n[:,nrep:].astype(float)**2).sum(1))[0]
            methods = [('bb_lognormal_F',output['pvalue']),('bb_lognormal_chi2',stats.chi2.sf(output['lr'],1)),('legacy_lr',old)]
            for method,p in methods:
                safe = np.where(np.isfinite(p),p,1.)
                q = multipletests(safe,method='fdr_bh')[1]
                found = q<=.05
                tp,fp = int((found&truth).sum()),int((found&~truth).sum())
                record = dict(label=label,seed=seed,sites=sites,nrep=nrep,depth=depth,imbalance=imbalance,
                    association=association,signal=signal,method=method,tp=tp,fp=fp,
                    fdp=fp/(tp+fp) if tp+fp else 0.,recall=tp/truth.sum() if truth.any() else None,
                    raw_null_fpr05=float((safe[~truth]<=.05).mean()),
                    raw_null_fpr001=float((safe[~truth]<=.001).mean()),
                    raw_null_fpr1e5=float((safe[~truth]<=1e-5).mean()),
                    average_precision=average_precision_score(truth,-safe) if truth.any() else None,
                    seconds=seconds if method.startswith('bb') else None,
                    prior_converged=output['prior']['converged'],prior_loglik=output['prior']['loglik'])
                records.append(record)
                print(json.dumps(record),flush=True)
                for mu in np.unique(mean):
                    for r in np.unique(rho):
                        mask = (mean==mu)&(rho==r)&~truth
                        if mask.any():
                            strata.append(dict(label=label,method=method,mean=mu,rho=r,sites=int(mask.sum()),
                                fpr05=float((safe[mask]<=.05).mean()),fpr001=float((safe[mask]<=.001).mean()),
                                fp=int((found&mask).sum())))
            pd.DataFrame(records).to_csv(OUT/'development_metrics.csv',index=False)
            pd.DataFrame(strata).to_csv(OUT/'development_null_strata.csv',index=False)
    (OUT/'development_context.json').write_text(json.dumps(dict(
        purpose='Development only; hold-out validation must use different seeds',sites=sites,
        seed_base=10021001,source=str(SOURCE),counts='Independent beta-binomial; Poisson depths',
        references=['F(1,nvalid-2)','chi2(1) diagnostic','frozen legacy lr'],
        region_correction='f73e66d unchanged',truth_used_by_inference=False),indent=2)+'\n')


if __name__=='__main__':main()
