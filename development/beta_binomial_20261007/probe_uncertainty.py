"""Throwaway comparison of dispersion-uncertainty approximations.

Not a released inference policy; calibration must decide whether any is useful.
"""
from pathlib import Path
import json
import sys
import time
ROOT=Path('/scratch/wgbs_benchmark_v4')
SOURCE=Path('/scratch/epykit-calibration-fresh_20261007/src')
sys.path[:0]=[str(ROOT/'vendor310'),str(SOURCE),str(Path(__file__).parent)]
import numpy as np
import pandas as pd
from numba import njit
from scipy import special,stats
from sklearn.metrics import average_precision_score
from statsmodels.stats.multitest import multipletests
from epykit import _beta_binomial as bb
from calibrate import generate
OUT=Path(__file__).parent

@njit(cache=True)
def profiles(m,n,ncase,grid):
    result=np.empty((len(m),len(grid),2))
    for i in range(len(m)):
        for b,rho in enumerate(grid):
            a,la,ia,_=bb._fit_mean(m[i],n[i],rho,0,ncase)
            c,lc,ic,_=bb._fit_mean(m[i],n[i],rho,ncase,m.shape[1])
            _,l0,_,_=bb._fit_mean(m[i],n[i],rho,0,m.shape[1])
            penalty=.5*((np.log(ia) if 0<a<1 else 0)+(np.log(ic) if 0<c<1 else 0))
            result[i,b]=2*max(0.,la+lc-l0),la+lc-penalty
    return result

rows=[]
for nrep,depth,imbalance,association in [(5,20,True,False),(5,60,True,False),(10,20,True,False),(5,20,True,True)]:
    for signal in [False,True]:
        label=f'n{nrep}_d{depth}_imb{imbalance}_assoc{association}_signal{signal}'
        prior=json.loads((OUT/f'prior_{label}.json').read_text())
        for key in ['rho_grid','rho_weights','mean_grid','mean_weights']:prior[key]=np.array(prior[key])
        m,n,truth,mean,rho=generate(9741001+nrep+depth+int(signal)+int(association)*101,30000,nrep,depth,imbalance,signal,association)
        started=time.monotonic()
        fit=profiles(m,n,nrep,prior['rho_grid'])
        print(label,'profile seconds',time.monotonic()-started,flush=True)
        conditional=stats.chi2.sf(fit[:,:,0],1)
        posterior=special.softmax(fit[:,:,1]+np.log(prior['rho_weights'])[None,:],axis=1)
        cr_p=(conditional*posterior).sum(1)
        # Exact full-model posterior rho given the discrete empirical mean prior.
        marginal=[]
        for a,z in [(m[:,:nrep],n[:,:nrep]),(m[:,nrep:],n[:,nrep:])]:
            L=bb._mixture_loglik(a,z,prior['mean_grid'],prior['rho_grid'])
            marginal.append(special.logsumexp(L+np.log(prior['mean_weights'])[None,:,None],axis=1))
        posterior_exact=special.softmax(marginal[0]+marginal[1]+np.log(prior['rho_weights'])[None,:],axis=1)
        exact_p=(conditional*posterior_exact).sum(1)
        for method,p in [('CR_posterior_tail',cr_p),('count_posterior_tail',exact_p)]:
            found=multipletests(p,method='fdr_bh')[1]<=.05
            tp,fp=int((found&truth).sum()),int((found&~truth).sum())
            row=dict(label=label,method=method,tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0,
                recall=tp/truth.sum() if truth.any() else None,
                fpr05=float((p[~truth]<=.05).mean()),fpr001=float((p[~truth]<=.001).mean()),
                fpr1e5=float((p[~truth]<=1e-5).mean()),
                AP=average_precision_score(truth,-p) if truth.any() else None,seconds=time.monotonic()-started)
            rows.append(row);print(json.dumps(row),flush=True)
        pd.DataFrame(rows).to_csv(OUT/'uncertainty_probe.csv',index=False)
