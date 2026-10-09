"""Throwaway count-variance mixture probe; not a production inference engine."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from diagnose import generate,accum
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import average_precision_score


def infer(groups):
    a,b=accum(groups[0]),accum(groups[1]);A,M,S,N=a;B,K,T,L=b
    A2=(groups[0][1]**2).sum(1);B2=(groups[1][1]**2).sum(1)
    mu1=M/A;mu2=K/B;mu=(M+K)/(A+B);d1=mu1*(1-mu1);d2=mu2*(1-mu2)
    c1=np.maximum(0,np.divide(S-M*M/A,d1,out=np.zeros(len(A)),where=d1>0))
    c2=np.maximum(0,np.divide(T-K*K/B,d2,out=np.zeros(len(B)),where=d2>0))
    informative1=d1>0;informative2=d2>0
    df=(N-1)*informative1+(L-1)*informative2
    exposure=(A-N-(A2/A-1))*informative1+(B-L-(B2/B-1))*informative2
    leverage=((A2/A-1)/A+(B2/B-1)/B)/(1/A+1/B)
    rho=np.r_[0,np.geomspace(.0005,.8,15)]
    labels=np.digitize(np.minimum(mu,1-mu),[.02,.05,.1,.2,.35])
    p=np.zeros(len(A))
    loglik=lambda x,n,pr: x*np.log(np.clip(pr,1e-9,1-1e-9))+(n-x)*np.log(np.clip(1-pr,1e-9,1-1e-9))
    stat=2*(loglik(M,A,mu1)+loglik(K,B,mu2)-loglik(M,A,mu)-loglik(K,B,mu))
    for label in np.unique(labels):
        sel=labels==label;deg=df[sel];Q=(c1+c2)[sel]
        scale=1+np.maximum(exposure[sel]/np.maximum(deg,1),0)[:,None]*rho
        ll=-(deg/2)[:,None]*np.log(scale)-Q[:,None]/(2*scale)
        weights=np.ones(len(rho))/len(rho)
        informative=deg>0
        for it in range(60):
            lp=ll+np.log(weights)
            posterior=np.exp(lp-logsumexp(lp,axis=1)[:,None])
            new=posterior[informative].mean(0)
            if np.max(np.abs(new-weights))<1e-6:break
            weights=new
        pp=chi2.sf(stat[sel,None]/(1+leverage[sel,None]*rho),1)
        p[sel]=(posterior*pp).sum(1)
    return p


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
    pd.DataFrame(rows).to_csv(Path(__file__).parent/'mixture_probe.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=='__main__':main()
