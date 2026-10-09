"""Independent oracles for removal of count-only likelihood constants."""
import numpy as np
import pytest
from scipy.special import gammaln
from scipy.stats import betabinom, binom

from epykit import _beta_binomial as bb


@pytest.mark.parametrize('rho', [0., 1e-8, .01, .4])
def test_relative_group_likelihood_matches_normalized_oracle(rho):
    assert hasattr(bb, '_relative_group_loglik')
    n = np.array([0, 5, 20, 300], dtype=np.int64)
    m = np.array([0, 1, 7, 271], dtype=np.int64)
    mu = .31
    if rho == 0:
        normalized = binom.logpmf(m, n, mu).sum()
    else:
        k = (1-rho)/rho
        normalized = betabinom.logpmf(m, n, mu*k, (1-mu)*k).sum()
    constants = (gammaln(n+1)-gammaln(m+1)-gammaln(n-m+1)).sum()
    actual = bb._relative_group_loglik(m, n, mu, rho, 0, len(m))
    np.testing.assert_allclose(actual, normalized-constants, atol=3e-6, rtol=1e-8)


@pytest.mark.parametrize('mu,m', [(0., [0,0]), (1., [5,300])])
def test_relative_group_endpoints_remain_exact(mu, m):
    assert hasattr(bb, '_relative_group_loglik')
    n = np.array([5,300], dtype=np.int64)
    assert bb._relative_group_loglik(np.array(m),n,mu,.1,0,2) == 0.


def test_public_loglik_and_intervals_keep_normalized_contract():
    m = np.array([[1,3,7,18,17,1,2,4,6,9]])
    n = np.full_like(m,20)
    prior = {'log_rho_mean':np.log(.06), 'log_rho_sd':1.2}
    result = bb.test_lognormal(m,n,5,prior=prior,intervals=True)
    rho = result['rho'][0]
    k = (1-rho)/rho
    a,b = result['mu_case'][0],result['mu_control'][0]
    expected = (betabinom.logpmf(m[0,:5],n[0,:5],a*k,(1-a)*k).sum()
                +betabinom.logpmf(m[0,5:],n[0,5:],b*k,(1-b)*k).sum())
    np.testing.assert_allclose(result['loglik_full'],[expected],atol=1e-9)
    assert result['ci_lo'][0] <= a-b <= result['ci_hi'][0]
    assert result['ci_hi'][0]-result['ci_lo'][0] < 2


def test_optimizer_diagnostics_cover_rho_search():
    result = bb.test_lognormal(np.array([[1,3,7,18,17,1,2,4,6,9]]),
        np.full((1,10),20),5,prior={'log_rho_mean':np.log(.06),'log_rho_sd':1.2})
    for key in ['rho_objective_evaluations','rho_search_converged',
                'rho_search_boundary','rho_objective_gap','rho_fallback_used']:
        assert key in result, key
        assert np.shape(result[key]) == (1,)
    assert result['rho_objective_evaluations'][0] >= 17
    assert result['rho_search_converged'][0]
    assert result['rho_objective_gap'][0] >= 0


def test_endpoint_dispersion_search_has_explicit_prior_only_diagnostics():
    result = bb.test_lognormal(np.zeros((1,4),dtype=int),np.full((1,4),20),2,
        prior={'log_rho_mean':np.log(.06),'log_rho_sd':1.2})
    assert 'rho_objective_evaluations' in result
    assert result['rho_objective_evaluations'][0] <= 3
    assert result['rho_search_converged'][0]
    assert not result['rho_fallback_used'][0]


def test_optimizer_diagnostics_reach_public_dmc_store(synth_md_filtered):
    import epykit as ep
    md = synth_md_filtered
    ep.tl.dmc(md, test='beta_binomial', chromosomes=['chr1'],
              beta_binomial_ci=False, tsv=False)
    assert {'bb_rho_evaluations','bb_rho_search_converged',
            'bb_rho_boundary','bb_rho_fallback_used','bb_rho_bracket_width'}.issubset(md.dmc.columns)
    assert md.dmc['bb_rho_evaluations'].min() >= 1
    assert md.dmc['bb_rho_search_converged'].all()
    assert 'local bracket' in md.uns['dmc']['beta_binomial_optimization']
