"""Throwaway probe: deconvolve count-residual noise before EB moderation."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from diagnose import generate,accum
import numpy as np
import pandas as pd
from scipy.stats import f
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import average_precision_score


def infer(groups):
    a,b=accum(groups[0]),accum(groups[1]);A,M,S,N=a;B,K,T,L=b
    A2=(groups[0][1]**2).sum(1);B2=(groups[1][1]**2).sum(1)
    mu1=M/A;mu2=K/B;mu=(M+K)/(A+B);d1=mu1*(1-mu1);d2=mu2*(1-mu2)
    c1=np.maximum(0,np.divide(S-M*M/A,d1,out=np.zeros(len(A)),where=d1>0))
    c2=np.maximum(0,np.divide(T-K*K/B,d2,out=np.zeros(len(B)),where=d2>0))
    info1=d1>0;info2=d2>0;nu=(N-1)*info1+(L-1)*info2
    ea=A-N-(A2/A-1);eb=B-L-(B2/B-1)
    coeff=(ea*info1+eb*info2)/np.maximum(nu,1)
    coeff=np.where(nu>0,coeff,(ea+eb)/(N+L-2))
    raw=(c1+c2)/np.maximum(nu,1)
    r=(raw-1)/np.maximum(coeff,1e-8)
    second=(raw*raw/(1+2/np.maximum(nu,1))-2*raw+1)/np.maximum(coeff*coeff,1e-8)
    leverage=((A2/A-1)/A+(B2/B-1)/B)/(1/A+1/B)
    labels=np.digitize(np.minimum(mu,1-mu),[.02,.05,.1,.2,.35])
    scale=np.zeros(len(A));dof=np.zeros(len(A))
    for lab in np.unique(labels):
        sel=labels==lab;fit=sel&(nu>0)&(coeff>0)
        if fit.sum()<100:fit=(nu>0)&(coeff>0)
        center=np.clip(np.mean(r[fit]),0,.95)
        between=max(float(np.mean(second[fit])-center**2),0.)
        prior_mean=1+coeff[sel]*center;prior_var=coeff[sel]**2*between
        shape=2+np.divide(prior_mean**2,prior_var,out=np.full_like(prior_mean,5e5),where=prior_var>1e-12)
        shape=np.minimum(shape,5e5);prior_scale=prior_mean*(shape-1)/shape
        posterior=(nu[sel]*np.maximum(raw[sel],0)+2*shape*prior_scale)/(nu[sel]+2*shape)
        rho=np.clip((posterior-1)/np.maximum(coeff[sel],1e-8),0,.95)
        scale[sel]=1+leverage[sel]*rho;dof[sel]=nu[sel]+2*shape
    loglik=lambda x,n,pr:x*np.log(np.clip(pr,1e-9,1-1e-9))+(n-x)*np.log(np.clip(1-pr,1e-9,1-1e-9))
    stat=2*(loglik(M,A,mu1)+loglik(K,B,mu2)-loglik(M,A,mu)-loglik(K,B,mu))
    return f.sf(stat/scale,1,dof)


def main():
    rows=[]
    for nrep,depth,imbalance in [(3,20,False),(5,20,False),(5,20,True),(5,60,True),(10,20,True)]:
        for signal in [False,True]:
            groups,truth,mean,rho=generate(731003+nrep+depth+int(signal),100000,nrep,depth,imbalance,signal)
            p=infer(groups);called=multipletests(p,method='fdr_bh')[1]<=.05
            tp=int(sum(called&truth));fp=int(sum(called&~truth))
            rows.append(dict(nrep=nrep,depth=depth,imbalance=imbalance,signal=signal,
                fpr05=float(np.mean(p[~truth]<.05)),fpr001=float(np.mean(p[~truth]<.001)),
                tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0,
                ap=average_precision_score(truth,-p) if truth.any() else None))
    pd.DataFrame(rows).to_csv(Path(__file__).parent/'count_prior_probe.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=='__main__':main()
