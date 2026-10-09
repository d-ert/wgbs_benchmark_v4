"""Throwaway count-mixture model comparison with explicit Bayesian-FDR semantics."""
from pathlib import Path
import json
import sys
import time
ROOT=Path('/scratch/wgbs_benchmark_v4');SOURCE=Path('/scratch/epykit-calibration-fresh_20261007/src')
OUT=Path(__file__).parent
sys.path[:0]=[str(ROOT/'vendor310'),str(SOURCE),str(OUT)]
import numpy as np
from scipy import special,optimize
from sklearn.metrics import average_precision_score
from epykit import _beta_binomial as bb
from calibrate import generate

def main():
    rows=[]
    for nrep,depth,association in [(3,20,False),(5,20,False),(5,60,False),(10,20,False),(5,20,True)]:
        imbalance=nrep!=3
        for signal in [False,True]:
            label=f'n{nrep}_d{depth}_imb{imbalance}_assoc{association}_signal{signal}'
            prior=json.loads((OUT/f'prior_{label}.json').read_text())
            for k in ['rho_grid','rho_weights','mean_grid','mean_weights']:prior[k]=np.array(prior[k])
            m,n,truth,mean,rho=generate(9861001+nrep+depth+int(signal)+101*int(association),50000,nrep,depth,imbalance,signal,association)
            start=time.monotonic()
            L1=bb._mixture_loglik(m[:,:nrep],n[:,:nrep],prior['mean_grid'],prior['rho_grid'])
            L0=bb._mixture_loglik(m[:,nrep:],n[:,nrep:],prior['mean_grid'],prior['rho_grid'])
            logmu=np.log(prior['mean_weights'])[None,:,None]
            logrho=np.log(prior['rho_weights'])[None,:]
            null=special.logsumexp(special.logsumexp(L1+L0+logmu,axis=1)+logrho,axis=1)
            alt=special.logsumexp(special.logsumexp(L1+logmu,axis=1)+special.logsumexp(L0+logmu,axis=1)+logrho,axis=1)
            logbf=alt-null
            def objective(pi):
                return -np.logaddexp(np.log(pi),np.log1p(-pi)+logbf).sum()
            fit=optimize.minimize_scalar(objective,bounds=(.5,1-1e-10),method='bounded')
            pi=fit.x
            if objective(1-1e-12)<fit.fun:pi=1-1e-12
            postnull=special.expit(np.log(pi)-np.log1p(-pi)-logbf)
            order=np.argsort(postnull)
            q=np.empty(len(m));q[order]=np.cumsum(postnull[order])/np.arange(1,len(m)+1)
            found=q<=.05
            tp,fp=int((found&truth).sum()),int((found&~truth).sum())
            row=dict(label=label,pi0=pi,tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0.,
                     recall=tp/truth.sum() if truth.any() else None,
                     AP=average_precision_score(truth,logbf) if truth.any() else None,
                     null_markov_tail05=float((logbf[~truth]>=np.log(20)).mean()),
                     null_markov_tail001=float((logbf[~truth]>=np.log(1000)).mean()),
                     seconds=time.monotonic()-start)
            rows.append(row);print(json.dumps(row),flush=True)
            (OUT/'bayes_factor_probe.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
