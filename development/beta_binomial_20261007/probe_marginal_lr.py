"""Throwaway genuine integrated-dispersion likelihood fitting experiment."""
from pathlib import Path
import json
import sys
import time
ROOT=Path('/scratch/wgbs_benchmark_v4')
SOURCE=Path('/scratch/epykit-calibration-fresh_20261007/src')
OUT=Path(__file__).parent
sys.path[:0]=[str(ROOT/'vendor310'),str(SOURCE),str(OUT)]
import numpy as np
from scipy import special,stats,optimize
from numba import njit
from sklearn.metrics import average_precision_score
from statsmodels.stats.multitest import multipletests
from epykit import _beta_binomial as bb
from calibrate import generate

@njit(cache=True)
def density(m,n,ncase,a,c,grid,weights):
    ll=np.empty(len(grid));ga=np.empty(len(grid));gc=np.empty(len(grid))
    ia=np.empty(len(grid));ic=np.empty(len(grid))
    for k,r in enumerate(grid):
        ll[k]=bb._group_loglik(m,n,a,r,0,ncase)+bb._group_loglik(m,n,c,r,ncase,len(m))+np.log(weights[k])
        ga[k],ia[k]=bb._score_information(m,n,a,r,0,ncase)
        gc[k],ic[k]=bb._score_information(m,n,c,r,ncase,len(m))
    shift=np.max(ll);post=np.exp(ll-shift);mass=np.sum(post);post/=mass
    g0=np.sum(post*ga);g1=np.sum(post*gc)
    h0=np.sum(post*(ga*ga-ia))-g0*g0
    h1=np.sum(post*(gc*gc-ic))-g1*g1
    h01=np.sum(post*ga*gc)-g0*g1
    return shift+np.log(mass),g0,g1,h0,h1,h01,np.sum(post*ia),np.sum(post*ic)

@njit(cache=True)
def fit(m,n,ncase,grid,weights,null):
    ma,na=np.sum(m[:ncase]),np.sum(n[:ncase]);mc,nc=np.sum(m[ncase:]),np.sum(n[ncase:])
    best=-np.inf;best_a=0.;best_c=0.
    for start in range(2):
        if start==0:
            a=ma/na;c=mc/nc
        else:
            a=np.mean(m[:ncase]/n[:ncase]);c=np.mean(m[ncase:]/n[ncase:])
        if null:a=c=(ma+mc)/(na+nc) if start==0 else np.mean(m/n)
        fixed_a=ma==0 or ma==na;fixed_c=mc==0 or mc==nc
        if null:fixed_a=fixed_c=ma+mc==0 or ma+mc==na+nc
        for iteration in range(40):
            ll,g0,g1,h0,h1,h01,i0,i1=density(m,n,ncase,a,c,grid,weights)
            if null:
                grad=g0+g1;hess=h0+h1+2*h01
                step=-grad/hess if hess<0 else grad/(i0+i1)
                da=dc=0. if fixed_a else step
            elif fixed_a and fixed_c:da=dc=0.
            elif fixed_a:da=0.;dc=-g1/h1 if h1<0 else g1/i1
            elif fixed_c:dc=0.;da=-g0/h0 if h0<0 else g0/i0
            else:
                det=h0*h1-h01*h01
                if h0<0 and h1<0 and det>0:
                    da=(-h1*g0+h01*g1)/det;dc=(h01*g0-h0*g1)/det
                else:da=g0/i0;dc=g1/i1
            if max(abs(da),abs(dc))<1e-9:break
            scale=1.
            accepted=False
            for _ in range(30):
                aa=a if fixed_a else min(1-1e-12,max(1e-12,a+scale*da))
                cc=c if fixed_c else min(1-1e-12,max(1e-12,c+scale*dc))
                nextll=density(m,n,ncase,aa,cc,grid,weights)[0]
                if nextll>=ll-1e-10:
                    a,c=aa,cc;accepted=True;break
                scale/=2
            if not accepted:break
        ll=density(m,n,ncase,a,c,grid,weights)[0]
        if ll>best:best,best_a,best_c=ll,a,c
    return best,best_a,best_c

@njit(cache=True)
def run(m,n,ncase,grid,weights):
    output=np.empty((len(m),3))
    for i in range(len(m)):
        full,a,c=fit(m[i],n[i],ncase,grid,weights,False)
        null,_,_=fit(m[i],n[i],ncase,grid,weights,True)
        output[i]=max(0.,2*(full-null)),a,c
    return output

def main():
    rows=[]
    for nrep,depth in [(5,20),(5,60),(10,20)]:
        for signal in [False,True]:
            label=f'n{nrep}_d{depth}_imbTrue_assocFalse_signal{signal}'
            prior=json.loads((OUT/f'prior_{label}.json').read_text())
            grid,weights=np.array(prior['rho_grid']),np.array(prior['rho_weights'])
            m,n,truth,mean,rho=generate(9791041+nrep+depth+int(signal),30000,nrep,depth,True,signal)
            started=time.monotonic();result=run(m,n,nrep,grid,weights)
            p=stats.chi2.sf(result[:,0],1)
            found=multipletests(p,method='fdr_bh')[1]<=.05
            tp,fp=int((found&truth).sum()),int((found&~truth).sum())
            row=dict(label=label,tp=tp,fp=fp,fdp=fp/(tp+fp) if tp+fp else 0.,
                recall=tp/truth.sum() if truth.any() else None,
                fpr05=float((p[~truth]<=.05).mean()),fpr001=float((p[~truth]<=.001).mean()),
                fpr1e5=float((p[~truth]<=1e-5).mean()),
                AP=average_precision_score(truth,-p) if truth.any() else None,seconds=time.monotonic()-started)
            rows.append(row);print(json.dumps(row),flush=True)
            (OUT/'marginal_lr_probe.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
