"""Development null-size and fixed-design count-bootstrap diagnostic.

Training and test counts use different seeds. Model parameters generate counts
only; the inference path receives observed counts/design and a learned prior.
This is not final holdout validation. Bootstrap p-values are pointwise model-
conditional diagnostics, not calibrated genome-wide analytical replacements.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control

from epykit import _beta_binomial as bb
from epykit._count_score import test_counts
from epykit._dispersion_prior import fit_dispersion_prior


def generate(rng,means,rhos,n):
    means,rhos = np.broadcast_arrays(np.asarray(means),np.asarray(rhos))
    means,rhos = np.broadcast_to(means,n.shape),np.broadcast_to(rhos,n.shape)
    positive = rhos > 0
    k = np.divide(1,rhos,out=np.ones(n.shape),where=positive)-1
    k[~positive] = 1.
    latent = rng.beta(np.clip(means*k,1e-10,None),np.clip((1-means)*k,1e-10,None))
    return rng.binomial(n,np.where(positive,latent,means))


def run():
    records,bootstrap = [],[]
    start = time.perf_counter()
    source_sha256 = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [Path(__file__)]+[Path(__file__).parent/'epykit/src/epykit'/name
                          for name in ('_count_score.py','_dispersion_prior.py','_beta_binomial.py')]}
    for replicates in (3,5,10):
        x = np.column_stack([np.ones(2*replicates),
                             np.r_[np.ones(replicates),np.zeros(replicates)]])
        training_rng = np.random.default_rng(202610092+replicates)
        train_n = training_rng.integers(5,101,size=(1200,2*replicates))
        train_mean = training_rng.choice([.02,.1,.25,.5],size=(len(train_n),1))
        train_rho = np.where(train_mean<=.1,.025,.12)
        train_m = generate(training_rng,train_mean,train_rho,train_n)
        training_start = time.perf_counter()
        prior = fit_dispersion_prior(train_m,train_n,replicates,max_training=1600,seed=991)
        training_seconds = time.perf_counter()-training_start
        print(json.dumps(dict(stage='prior ready',replicates=replicates,
                              elapsed_seconds=time.perf_counter()-start,
                              training_seconds=training_seconds)),flush=True)
        prior_summary = dict(bin_rho=(np.array(prior['rho_weights_by_bin'])@
                                      prior['rho_grid']).tolist(),
                             training_digest=prior['training_digest'],
                             converged=prior['bin_converged'])
        rng = np.random.default_rng(202610102+replicates)
        for mean in (.02,.1,.5):
            for rho in (0.,.025,.12):
                for pattern in ('equal','unequal'):
                    depths = (np.full(2*replicates,20) if pattern=='equal' else
                              np.tile(np.linspace(5,100,replicates).astype(int),2))
                    n = np.broadcast_to(depths,(5000,len(depths))).copy()
                    m = generate(rng,mean,rho,n)
                    for mode in ('shared','group'):
                        fit_start = time.perf_counter()
                        fit = test_counts(m,n,x,[0.,1.],prior=prior,n_case=replicates,
                                          dispersion_mode=mode,intervals=False)
                        fit_seconds = time.perf_counter()-fit_start
                        p = np.nan_to_num(fit['pvalue'],nan=1.)
                        records.append(dict(replicates=replicates,mean=mean,rho=rho,
                            depths=pattern,mode=mode,sites=len(m),prior=prior_summary,
                            fit_seconds=fit_seconds,training_seconds=training_seconds,
                            estimable=int(fit['estimable'].sum()),median_rho=float(np.median(fit['rho'])),
                            iteration_converged=float(fit['dispersion_converged'].mean()),
                            minimum_p=float(p.min()),
                            rejection_rates={str(a):float((p<=a).mean()) for a in (.05,.01,.001)},
                            bh_calls_05=int((false_discovery_control(p)<=.05).sum())))
        if replicates == 5:
            # Three observed fixtures with fixed original depths, one of them
            # containing a moderate contrast. Never hold coverage at its mean.
            n = np.tile([5,10,20,50,100,5,10,20,50,100],(3,1))
            means = np.array([[.1]*10,[.5]*10,[.6]*5+[.35]*5])
            m = generate(rng,means,.12,n)
            original = test_counts(m,n,x,[0.,1.],prior=prior,n_case=5,intervals=False)
            bb_prior = bb.project_log_prior(dict(rho_grid=np.array(prior['rho_grid']),
                                       rho_weights=np.array(prior['pooled_rho_weights'])))
            bb_fit = bb.test_lognormal(m,n,5,prior=bb_prior,intervals=False)
            for i in range(len(m)):
                size = 4000
                bn = np.broadcast_to(n[i],(size,n.shape[1])).copy()
                bm = generate(rng,original['mu_null'][i],original['rho'][i],bn)
                simulated = test_counts(bm,bn,x,[0.,1.],prior=prior,n_case=5,intervals=False)
                statistic = np.nan_to_num(simulated['statistic'],nan=0.)
                exceed = int((statistic>=original['statistic'][i]).sum())
                p = (exceed+1)/(size+1)
                bootstrap.append(dict(observed_m=m[i].tolist(),fixed_n=n[i].tolist(),
                    observed_statistic=float(original['statistic'][i]),
                    normal_p=float(original['pvalue'][i]),bb_f_p=float(bb_fit['pvalue'][i]),
                    null_count_bootstrap_p=p,monte_carlo_se=float(np.sqrt(p*(1-p)/(size+1))),
                    resamples=size,tail_exceedances=exceed,
                    generated_null_rho=float(original['rho'][i]),
                    limits='Conditional plug-in null parameters; each replicate refits rho '
                           'against the frozen observed-count prior. Not a posterior-predictive '
                           'integration over parameter uncertainty or genome-wide tail guarantee.'))
    artifact = dict(kind='estimated-rho development diagnostic',
        source_sha256=source_sha256,
        seeds='202610092+replicates (training); 202610102+replicates (test)',
        elapsed_seconds=time.perf_counter()-start,records=records,bootstrap=bootstrap,
        recommendation='Experimental only pending independent held-out calibration.',
        limits='Independent sites, one realization per scenario; BH calls do not estimate FDR.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(artifact,indent=2)+'\n')
    print(json.dumps(dict(scenarios=len(records),elapsed_seconds=artifact['elapsed_seconds'],
        worst_size_05=max(records,key=lambda r:r['rejection_rates']['0.05']),
        worst_size_001=max(records,key=lambda r:r['rejection_rates']['0.001']),bootstrap=bootstrap)))


if __name__=='__main__':
    run()
