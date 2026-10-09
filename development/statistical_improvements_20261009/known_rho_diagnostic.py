"""Development-only oracle diagnostic: rho is deliberately supplied as known.

This diagnoses the analytical reference separately from dispersion estimation.
It is not the held-out evaluation and must never train the production prior.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control

from epykit._count_score import test_score


def run():
    rng = np.random.default_rng(202610091)
    size = 5000
    records = []
    start = time.perf_counter()
    for replicates in (3, 5, 10):
        x = np.column_stack([np.ones(2*replicates),
                             np.r_[np.ones(replicates), np.zeros(replicates)]])
        for mean in (.01, .1, .5):
            for rho in (0., .05, .2):
                for pattern in ('equal', 'unequal'):
                    depths = (np.full(2*replicates, 20) if pattern == 'equal' else
                              np.tile(np.linspace(5, 100, replicates).astype(int), 2))
                    n = np.broadcast_to(depths, (size, len(depths))).copy()
                    if rho:
                        k = 1/rho-1
                        latent = rng.beta(mean*k, (1-mean)*k, size=n.shape)
                    else:
                        latent = mean
                    m = rng.binomial(n, latent)
                    result = test_score(m, n, x, np.array([0., 1.]), rho, intervals=False)
                    p = np.nan_to_num(result['pvalue'], nan=1.)
                    records.append(dict(replicates=replicates, mean=mean, rho=rho,
                        depths=pattern, sites=size, estimable=int(result['estimable'].sum()),
                        minimum_p=float(p.min()),
                        rejection_rates={str(a): float((p <= a).mean())
                                         for a in (.05, .01, .001, .0001)},
                        bh_calls_05=int((false_discovery_control(p) <= .05).sum())))
    artifact = dict(kind='known-rho analytical reference development diagnostic',
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),Path(__file__).parent/'epykit/src/epykit/_count_score.py')},
        seed=202610091, sites_per_scenario=size, elapsed_seconds=time.perf_counter()-start,
        limits='Independent CpGs; fixed known rho; one realization per scenario. '
               'BH calls on one null realization do not estimate FDR.', records=records)
    Path(__file__).with_suffix('.json').write_text(json.dumps(artifact, indent=2)+'\n')
    print(json.dumps(dict(scenarios=len(records), elapsed_seconds=artifact['elapsed_seconds'],
        worst_size_05=max(records, key=lambda r:r['rejection_rates']['0.05']),
        worst_size_001=max(records, key=lambda r:r['rejection_rates']['0.001']))))


if __name__ == '__main__':
    run()
