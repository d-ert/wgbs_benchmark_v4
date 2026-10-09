"""Regression checks for dropped extreme evidence and incomplete null scans."""
import numpy as np
import pytest

from epykit import dmr
from types import SimpleNamespace


def test_zero_tail_is_retained():
    p = dmr._stouffer_combine_signed(np.array([0., .5]), np.ones(2))
    assert np.isfinite(p) and p < 1e-100
    assert np.isfinite(dmr._stouffer_combine_signed(np.zeros(2), np.ones(2)))


def test_invalid_probabilities_and_effects_are_excluded():
    assert np.isnan(dmr._stouffer_combine_signed(
        np.array([np.nan, -.1, 1.1, np.inf]), np.ones(4)))
    assert dmr._stouffer_combine_signed(
        np.array([.5, 1e-20]), np.array([1., np.inf])) == pytest.approx(.5)


def test_invalid_weights_do_not_add_evidence():
    assert dmr._stouffer_combine_signed(
        np.array([.5, 1e-20]), np.ones(2), np.array([1., np.nan])) == pytest.approx(.5)
    assert np.isnan(dmr._stouffer_combine_signed(
        np.array([.5, .5]), np.ones(2), np.array([-1., np.inf])))


@pytest.mark.parametrize('results,n_perm', [([(False, None)], 1), ([], 1)])
def test_region_count_ratio_refuses_failed_or_missing_scan(results, n_perm):
    with pytest.raises(RuntimeError, match='Incomplete'):
        dmr._aggregate_region_perm_results(observed_pvalues=np.array([.01]),
            results=results, n_perm=n_perm, fdr_method='region')


def test_clean_empty_and_self_scans_keep_their_distinct_contracts():
    p, q, set_fdr = dmr._aggregate_region_perm_results(
        observed_pvalues=np.array([.01, .02]),
        results=[(False, np.array([.005])), (False, np.array([])),
                 (True, np.array([1e-9]))], n_perm=3, fdr_method='region')
    np.testing.assert_allclose(q, [.25, .25])
    assert set_fdr == pytest.approx(.25)


def test_permutation_diagnostics_report_used_scans_and_uncertainty():
    assert hasattr(dmr, '_region_permutation_diagnostics')
    result = dmr._region_permutation_diagnostics(
        np.array([.01, .02, .03, .04]),
        [(False, np.array([])), (False, np.full(4, .1)),
         (False, np.full(8, .1)), (True, np.array([1e-9]))], 'region')
    assert result['empirical_n_perm_requested'] == 4
    assert result['empirical_n_perm_used'] == 3
    assert result['empirical_null_regions'] == 12
    assert result['empirical_fdr_set_mc_se'] == pytest.approx(1/np.sqrt(3))
    assert result['empirical_inference'] == 'count-ratio estimate; not a finite-sample FDR guarantee'


def test_invalid_null_scores_cannot_disappear_from_region_pool():
    with pytest.raises(ValueError, match='Invalid'):
        dmr._aggregate_region_perm_results(observed_pvalues=np.array([.01]),
            results=[(False, np.array([np.nan]))], n_perm=1, fdr_method='region')


def test_changed_regional_score_invalidates_both_caches(monkeypatch):
    store = SimpleNamespace(manifest={'input_sig': 'same-count-data'})
    sliding = lambda: dmr._dmr_sliding_cache_key(store, 500, 5, 3, .05, .1, 'pvalue')
    chain = lambda: dmr._dmr_chain_merge_cache_key(store, .05, .1, 500, 5, .5, 50, 'pvalue')
    old = sliding(), chain()
    monkeypatch.setattr(dmr, 'REGIONAL_SCORE_REVISION', 'different-score', raising=False)
    assert (sliding(), chain()) != old
