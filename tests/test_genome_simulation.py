import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from simulation.simulate import Scenario, simulate, coverage_parameters, _counts
from simulation.validate_streaming import validate_dataset


def calibration_fixture(root):
    root.mkdir()
    files = []
    for chrom in ('chr1', 'chr2'):
        path = root / f'{chrom}_00001.parquet'
        pq.write_table(pa.table({'pos': np.arange(1, 12001, 10), 'mu': np.full(1200, .4),
            'coverage_mean': np.linspace(5, 50, 1200), 'zero_fraction': np.linspace(0, .4, 1200)}), path)
        files.append(str(path))
    (root / 'manifest.json').write_text(json.dumps({'release_ready': False,
        'parquet_files': files, 'tau_by_mean_bin': [.1]*99, 'sample_mean_depth': [10, 20, 40],
        'spatial_correlation': {'fitted_length_bp': 100, 'validated': False}}))
    return root


def read_counts(root, sample='control_1'):
    return pd.read_csv(root/'counts'/f'{sample}.cov.gz', sep='\t', header=None,
                       names=['chrom','pos','end','percent','m','u']).set_index(['chrom','pos'])


def test_weighted_dropout_normalization():
    frame = pd.DataFrame({'coverage_mean': [5., 100.], 'zero_fraction': [.1, .8]})
    depth, dropout = coverage_parameters(frame)
    assert np.isclose(np.mean(depth * (1-dropout)), 30)


def test_reproducible_simulation_and_nested_thinning(tmp_path):
    cal = calibration_fixture(tmp_path/'cal')
    paths = [tmp_path/x for x in ('a','b','thin','signal')]
    for path, coverage, mode in zip(paths, (30,30,10,30), ('null','null','null','dmr')):
        simulate(cal,path,Scenario(771,coverage=coverage,mode=mode,regions=4),chromosomes=['chr1','chr2'])
        assert validate_dataset(path, ['chr1','chr2'])['passed']
    a,b,thin,signal = paths
    assert (a/'counts/control_1.cov.gz').read_bytes() == (b/'counts/control_1.cov.gz').read_bytes()
    high = read_counts(a)
    low = read_counts(thin).reindex(high.index,fill_value=0)
    assert (low.m <= high.m).all() and (low.u <= high.u).all()
    # Coverage draws cannot change when effects are changed.
    sig = read_counts(signal,'treatment_1')
    null = read_counts(a,'treatment_1')
    pd.testing.assert_series_equal(sig.m+sig.u,null.m+null.u)
    factors=json.loads((a/'manifest.json').read_text())['library_factors']
    assert len(factors)==10 and np.isclose(np.mean(list(factors.values())),1)


def test_beta_binomial_sampling_moments():
    size=100000
    mu,tau,n=.3,.2,25
    x=_counts(np.full(size,mu),np.full(size,tau),np.full(size,n),np.random.default_rng(991),
              False,np.arange(size),None,count_rng=np.random.default_rng(992))
    rho=tau/(1+tau)
    assert abs(x.mean()-n*mu)<.08
    assert abs(x.var()-n*mu*(1-mu)*(1+(n-1)*rho))<.4


def test_validator_rejects_corrupt_counts(tmp_path):
    import gzip
    import pytest
    cal=calibration_fixture(tmp_path/'cal')
    out=tmp_path/'null'
    simulate(cal,out,Scenario(771,mode='null'),chromosomes=['chr1','chr2'])
    target=out/'counts/control_1.cov.gz'
    with gzip.open(target,'rt') as f:
        rows=f.readlines()
    fields=rows[0].rstrip().split('\t')
    fields[4]='-1'
    rows[0]='\t'.join(fields)+'\n'
    with gzip.open(target,'wt') as f:f.writelines(rows)
    with pytest.raises(AssertionError):
        validate_dataset(out,['chr1','chr2'])


def test_parallel_calibration_combines_statistics_not_averages(tmp_path, monkeypatch):
    from simulation import calibrate_parallel as c
    def fake_worker(args):
        out=args[1];out.mkdir(parents=True)
        first=args[4]==['chr1']
        count=100 if first else 300
        hist=np.zeros((99,400),dtype=np.int64)
        hist[:,10 if first else 300]=1 if first else 3
        np.savez(out/'sufficient_statistics.npz',hist=hist,
                 depth_sum=np.array([200.] if first else [3000.]),
                 spatial_sum=np.zeros(7),spatial_count=np.zeros(7,dtype=np.int64))
        return dict(counts={'all':count,'valid':count},release_ready=True,
                    spatial_correlation={'gap_edges_bp':[0,25,50,100,200,500,1000,2000]},
                    parquet_files=[str(out/'part.parquet')])
    class Pool:
        def __init__(self,**kwargs):pass
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def map(self,fn,jobs):return list(map(fake_worker,jobs))
    monkeypatch.setattr(c,'ProcessPoolExecutor',Pool)
    result=c.calibrate_genome(None,tmp_path/'cal',[],None,['chr1','chr2'],None,workers=2)
    assert result['counts']['all']==400
    assert result['sample_mean_depth']==[8.]
    edges=np.linspace(-4,np.log10(4),401)
    assert np.allclose(result['tau_by_mean_bin'],10**((edges[300]+edges[301])/2))
    assert result['scope']['chromosomes']==['chr1','chr2']
