"""Untouched-seed evaluation of the frozen experimental BB-F candidate.

No genome or paper truth is read. Alternate prior starts are sensitivity
analyses, never used to choose the reported candidate or reference.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path('/scratch/wgbs_benchmark_v4')
SOURCE = Path('/scratch/epykit-calibration-fresh_20261007/src')
OUT = Path(__file__).parent / 'holdout_v3'
sys.path[:0] = [str(ROOT/'vendor310'), str(SOURCE)]
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score
from statsmodels.stats.multitest import multipletests
from epykit import _beta_binomial as bb
from calibrate_continuous import generate, serial


def evaluate(p, truth, eligible):
    p = np.where(eligible & np.isfinite(p),p,1.)
    q = multipletests(p[eligible],method='fdr_bh')[1]
    found = np.zeros(len(p),dtype=bool)
    found[eligible] = q<=.05
    tp,fp = int((found&truth).sum()),int((found&~truth).sum())
    null = eligible & ~truth
    result = dict(tp=tp,fp=fp,any_rejection=bool(found.any()),
                  fdp=fp/(tp+fp) if tp+fp else 0.,
                  recall=tp/(truth&eligible).sum() if (truth&eligible).any() else None,
                  tested_sites=int(eligible.sum()),positive_sites=int((truth&eligible).sum()),
                  average_precision=float(average_precision_score(truth[eligible],-p[eligible]))
                    if (truth&eligible).any() else None)
    for cutoff in [.05,.001,1e-5]:
        count = int((p[null]<=cutoff).sum())
        n = int(null.sum())
        lo,hi = stats.binomtest(count,n).proportion_ci(.95)
        result[f'null_tail_{cutoff:g}'] = count/n
        result[f'null_tail_{cutoff:g}_lo'] = float(lo)
        result[f'null_tail_{cutoff:g}_hi'] = float(hi)
    return result,found


def main():
    if OUT.exists():
        raise ValueError(f'Preserving existing holdout evaluation: {OUT}')
    OUT.mkdir()
    source_hash = hashlib.sha256((SOURCE/'epykit/_beta_binomial.py').read_bytes()).hexdigest()
    record = dict(status='running',seed_base=71082001,algorithm=bb.ALGORITHM_REVISION,
                  source_sha256=source_hash,truth_used_by_inference=False,
                  purpose='Held-out calibration and prior-start sensitivity; no retuning')
    (OUT/'context.json').write_text(json.dumps(record,indent=2)+'\n')
    scenarios = []
    # Three independent whole-null experiments at each replicate count.
    for nrep in [3,5,10]:
        for repeat in range(3):
            scenarios.append((f'null_n{nrep}_repeat{repeat}',25000,nrep,20,True,False,False,'none'))
    for label,nrep,depth,association,stress in [
        ('balanced',5,20,False,'none'),('high_depth',5,60,False,'none'),
        ('many_replicates',10,20,False,'none'),('mean_rho_association',5,20,True,'none'),
        ('unequal_group_rho',5,20,False,'unequal_rho'),
        ('missing_coverage',5,20,False,'missing'),
        ('strong_effect',5,60,False,'strong_effect')]:
        for signal in [False,True]:
            scenarios.append((f'{label}_signal{signal}',75000,nrep,depth,label!='balanced',signal,association,stress))
    rows,sensitivity = [],[]
    for index,(label,sites,nrep,depth,imbalance,signal,association,stress) in enumerate(scenarios):
        seed = 71082001+index*151
        m,n,truth,mean,rho = generate(seed,sites,nrep,depth,imbalance,signal,association)
        rng = np.random.default_rng(seed+8000000)
        if stress in {'unequal_rho','strong_effect'}:
            case_mean = np.clip(mean+np.where(mean<=.5,.6 if stress=='strong_effect' else .2,
                                             -.6 if stress=='strong_effect' else -.2)*truth,.001,.999)
            case_rho = np.minimum(.8,rho*3) if stress=='unequal_rho' else rho
            k = (1-case_rho)/case_rho
            latent = rng.beta((case_mean*k)[:,None],((1-case_mean)*k)[:,None],size=(sites,nrep))
            m[:,:nrep] = rng.binomial(n[:,:nrep],latent)
        if stress=='missing':
            missing = rng.random(n.shape)<.15
            m[missing] = 0
            n[missing] = 0
        eligible = ((n[:,:nrep]>0).sum(1)>=2)&((n[:,nrep:]>0).sum(1)>=2)
        start = time.monotonic()
        prior = bb.estimate_prior(m,n,nrep)
        result = bb.test_lognormal(m,n,nrep,prior=prior,reference='F')
        values,found = evaluate(result['pvalue'],truth,eligible)
        row = dict(label=label,seed=seed,sites=sites,nrep=nrep,depth=depth,signal=signal,
                   stress=stress,prior_converged=bool(prior['converged']),
                   log_rho_mean=prior['log_rho_mean'],log_rho_sd=prior['log_rho_sd'],
                   seconds=time.monotonic()-start,**values)
        rows.append(row)
        print(json.dumps(row),flush=True)
        (OUT/f'prior_{label}.json').write_text(json.dumps(prior,default=serial,indent=2)+'\n')
        # Both original deterministic starts: quantify, never substitute.
        for j,weights in enumerate(prior['start_rho_weights']):
            alternative = bb.project_log_prior(dict(prior,rho_weights=weights))
            alt = (result if np.array_equal(weights,prior['rho_weights']) else
                   bb.test_lognormal(m,n,nrep,prior=alternative,reference='F'))
            metrics,alt_found = evaluate(alt['pvalue'],truth,eligible)
            valid = eligible & np.isfinite(result['pvalue']) & np.isfinite(alt['pvalue'])
            difference = np.abs(np.log10(np.maximum(result['pvalue'][valid],1e-300))
                                -np.log10(np.maximum(alt['pvalue'][valid],1e-300)))
            sensitivity.append(dict(label=label,start=j,loglik=prior['start_logliks'][j],
                log_rho_mean=alternative['log_rho_mean'],log_rho_sd=alternative['log_rho_sd'],
                changed_calls=int((found!=alt_found).sum()),
                abs_log10_p_difference_p95=float(np.quantile(difference,.95)),**metrics))
        pd.DataFrame(rows).to_csv(OUT/'metrics.csv',index=False)
        pd.DataFrame(sensitivity).to_csv(OUT/'prior_start_sensitivity.csv',index=False)
    null_trials = [r for r in rows if r['label'].startswith('null_n')]
    k = sum(r['any_rejection'] for r in null_trials)
    ci = stats.binomtest(k,len(null_trials)).proportion_ci(.95)
    record.update(status='complete',sites_total=sum(r['sites'] for r in rows),
                  repeated_null_trials=len(null_trials),repeated_null_trials_with_rejection=k,
                  null_any_rejection_95ci=[float(ci.low),float(ci.high)])
    assert hashlib.sha256((SOURCE/'epykit/_beta_binomial.py').read_bytes()).hexdigest()==source_hash
    (OUT/'context.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    main()
