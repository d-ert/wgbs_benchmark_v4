"""Bounded warm-kernel diagnostics; timings are not genome benchmark reruns."""
from pathlib import Path
import sys, os, importlib.util, json, time
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor310'))
os.environ['NUMBA_CACHE_DIR']='/tmp/epykit-report-numba'
import numpy as np
source=ROOT/'results/beta_binomial_20261007_full/source_snapshot/epykit/_beta_binomial.py'
spec=importlib.util.spec_from_file_location('report_bb',source)
bb=importlib.util.module_from_spec(spec);spec.loader.exec_module(bb)
rng=np.random.default_rng(202610091)
prior=dict(log_rho_mean=np.log(.1),log_rho_sd=1.)
rows=[]
for depth in [10,20,80]:
    n=np.full((5000,10),depth,dtype=np.int64)
    latent=rng.beta(4.5,4.5,size=n.shape)
    m=rng.binomial(n,latent).astype(np.int64)
    rho=np.full(len(m),.1)
    # Compile/warm both kernels before measurement.
    bb.fit_fixed_dispersion(m[:4],n[:4],5,rho[:4])
    bb.test_lognormal(m[:4],n[:4],5,prior=prior)
    for method,fn in [
        ('fixed_rho',lambda:bb.fit_fixed_dispersion(m,n,5,rho)),
        ('continuous_MAP_rho',lambda:bb.test_lognormal(m,n,5,prior=prior))]:
        times=[]
        for _ in range(3):
            t=time.perf_counter();result=fn();times.append(time.perf_counter()-t)
        rows.append(dict(depth=depth,method=method,sites=len(m),median_seconds=float(np.median(times)),repeats=times))
    # Count objective evaluations using the Python outer loop calling the
    # same compiled inner objective. This is a call count, not a speed test.
    original=bb._log_dispersion_objective
    counter=[0]
    def counted(*args):
        counter[0]+=1
        return original(*args)
    bb._log_dispersion_objective=counted
    bb._lognormal_batch.py_func(m[:100],n[:100],5,prior['log_rho_mean'],prior['log_rho_sd'])
    bb._log_dispersion_objective=original
    rows[-1]['mean_objective_evaluations_per_site']=counter[0]/100
    # A small confidence-interval subset isolates an option skipped in the
    # experimental full benchmark, so no timing is conflated with that run.
    bb.test_lognormal(m[:2],n[:2],5,prior=prior,intervals=True)
    t=time.perf_counter();bb.test_lognormal(m[:100],n[:100],5,prior=prior,intervals=True)
    rows.append(dict(depth=depth,method='MAP_plus_conditional_CI',sites=100,median_seconds=time.perf_counter()-t,repeats=None))
(OUT/'bb_kernel_profile.json').write_text(json.dumps(dict(seed=202610091,prior=prior,notes='5000 synthetic independent sites, 5+5 replicates, warmed compiled kernels, three timing repeats. Objective counts on 100 sites use Python outer loop with compiled objective. No performance extrapolation or inferential equivalence claimed.',rows=rows),indent=2)+'\n')
print(json.dumps(rows,indent=2))
