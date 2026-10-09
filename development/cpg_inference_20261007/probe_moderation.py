"""Throwaway analytic probe of scaled-F variance moderation; not package code."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from diagnose import generate
import numpy as np
import pandas as pd
from scipy.special import digamma,polygamma
from scipy.optimize import brentq
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import average_precision_score


def moderate(x,df):
    x=np.maximum(x,1e-5*np.median(x[x>0]))
    e=np.log(x)+np.log(df/2)-digamma(df/2)
    noise=float(np.var(e,ddof=1)-polygamma(1,df/2))
    if noise<=0:return np.full(len(x),x.mean()),np.full(len(x),1e6)
    df0=2*brentq(lambda z:float(polygamma(1,z))-noise,1e-8,1e8)
    prior=np.exp(e.mean()+digamma(df0/2)-np.log(df0/2))
    return (df*x+df0*prior)/(df+df0),np.full(len(x),df+df0)


def main():
    rows=[]
    for nrep,depth,imbalance in [(3,20,False),(5,20,False),(5,20,True),(5,60,True),(10,20,True)]:
        for signal in [False,True]:
            groups,truth,mean,rho=generate(731003+nrep+depth+int(signal),100000,nrep,depth,imbalance,signal)
            y1=groups[0][0]/groups[0][1];y2=groups[1][0]/groups[1][1]
            mu1=y1.mean(1);mu2=y2.mean(1);mu=(mu1+mu2)/2
            normalized=np.maximum(mu*(1-mu),1e-6)
            sample=((y1-mu1[:,None])**2).sum(1)+((y2-mu2[:,None])**2).sum(1)
            sample/=2*nrep-2
            sampling=(1/groups[0][1]).mean(1)+(1/groups[1][1]).mean(1);sampling/=2
            for method in ['logF','floored_logF','strata_floored_logF']:
                values=sample/normalized
                if method!='logF':values=np.maximum(values,sampling)
                labels=np.digitize(np.minimum(mu,1-mu),[.02,.05,.1,.2,.35]) if method.startswith('strata') else np.zeros(len(mu),int)
                variance=np.zeros(len(mu));dof=np.zeros(len(mu))
                for lab in np.unique(labels):
                    select=labels==lab;variance[select],dof[select]=moderate(values[select],2*nrep-2)
                variance=np.maximum(variance,sampling)*normalized
                z=(mu1-mu2)/np.sqrt(variance*2/nrep)
                p=2*t.sf(np.abs(z),dof);called=multipletests(p,method='fdr_bh')[1]<=.05
                tp=int(sum(called&truth));fp=int(sum(called&~truth))
                rows.append(dict(nrep=nrep,depth=depth,imbalance=imbalance,signal=signal,method=method,
                    fpr05=float(np.mean(p[~truth]<.05)),fpr001=float(np.mean(p[~truth]<.001)),
                    tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0,
                    ap=average_precision_score(truth,-p) if truth.any() else None))
    pd.DataFrame(rows).to_csv(Path(__file__).parent/'moderation_probe.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=='__main__':main()
