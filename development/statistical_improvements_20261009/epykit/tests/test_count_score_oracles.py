"""Known-dispersion quasi-score checks derived independently of the kernel."""

import importlib
import importlib.util

import numpy as np
import pytest
from scipy.special import expit
from scipy.stats import norm


def core():
    assert importlib.util.find_spec("epykit._count_score") is not None, (
        "Count-score kernel is missing"
    )
    return importlib.import_module("epykit._count_score")


def design(n_case=3):
    return np.column_stack([np.ones(2 * n_case), np.r_[np.ones(n_case), np.zeros(n_case)]])


@pytest.mark.parametrize("rho", [0.0, 0.1, 0.4])
def test_efficient_score_matches_hand_computed_binomial_inflation(rho):
    m = np.array([[16, 15, 17, 4, 5, 3]])
    n = np.full_like(m, 20)
    r = core().test_score(m, n, design(), np.array([0.0, 1.0]), rho)
    inflation = 1 + 19 * rho
    # Null mu=.5; treatment score=18/inflation. Removing intercept
    # information leaves 7.5/inflation, giving statistic 43.2/inflation.
    np.testing.assert_allclose(r["score"], [18 / inflation], rtol=1e-8)
    np.testing.assert_allclose(r["variance"], [7.5 / inflation], rtol=1e-8)
    np.testing.assert_allclose(r["statistic"], [43.2 / inflation], rtol=1e-8)
    np.testing.assert_allclose(r["pvalue"], [2 * norm.sf(np.sqrt(43.2 / inflation))], rtol=1e-8)
    np.testing.assert_allclose(r["effect"], [0.6], atol=1e-7)
    assert r["ci_lo"][0] < 0.6 < r["ci_hi"][0]


def test_quasi_mean_uses_saturating_read_weights():
    n = np.array([[5, 5, 5, 5, 100, 5, 5, 5, 5, 100]])
    m = np.array([[1, 1, 1, 1, 90, 1, 1, 1, 1, 10]])
    rho = 0.1
    r = core().fit_quasi_means(m, n, design(5), rho)
    weights = n[0] / (1 + (n[0] - 1) * rho)
    np.testing.assert_allclose(r["effective_weights"][0], weights)
    assert weights[4] / weights[0] == pytest.approx((100 / 10.9) / (5 / 1.4))
    a = (4 / 1.4 + 90 / 10.9) / (20 / 1.4 + 100 / 10.9)
    b = (4 / 1.4 + 10 / 10.9) / (20 / 1.4 + 100 / 10.9)
    np.testing.assert_allclose(r["mu"][0], [a] * 5 + [b] * 5, atol=1e-8)
    assert abs(a - m[0, :5].sum() / n[0, :5].sum()) > 0.25


def test_donor_adjusted_mean_matches_weighted_statsmodels():
    import statsmodels.api as sm

    x = np.column_stack(
        [np.ones(8), np.r_[np.ones(4), np.zeros(4)], np.tile(np.eye(4)[:, 1:], (2, 1))]
    )
    n = np.array([[10, 20, 30, 40, 15, 25, 35, 45]])
    m = np.array([[8, 12, 21, 29, 4, 9, 14, 22]])
    rho = 0.07
    oracle = sm.GLM(
        np.column_stack([m[0], n[0] - m[0]]),
        x,
        family=sm.families.Binomial(),
        var_weights=1 / (1 + (n[0] - 1) * rho),
    ).fit()
    r = core().fit_quasi_means(m, n, x, rho)
    np.testing.assert_allclose(r["coefficients"][0], oracle.params, atol=1e-7)
    np.testing.assert_allclose(r["mu"][0], oracle.fittedvalues, atol=1e-8)


def test_nuisance_adjusted_score_matches_independent_null_glm_projection():
    import statsmodels.api as sm

    group = np.r_[np.ones(4), np.zeros(4)]
    batch = np.array([0.0, 1.0, 2.0, 1.0, 1.0, 0.0, 2.0, 0.0])
    nuisance = np.column_stack([np.ones(8), batch])
    x = np.column_stack([nuisance, group])
    n = np.array([[10, 20, 30, 40, 15, 25, 35, 45]])
    m = np.array([[8, 12, 21, 29, 4, 9, 14, 22]])
    rho = np.array([0.05, 0.08, 0.06, 0.12, 0.04, 0.1, 0.09, 0.07])
    d = 1 + (n[0] - 1) * rho
    null = sm.GLM(
        np.column_stack([m[0], n[0] - m[0]]),
        nuisance,
        family=sm.families.Binomial(),
        var_weights=1 / d,
    ).fit()
    residual = (m[0] - n[0] * null.fittedvalues) / d
    w = n[0] * null.fittedvalues * (1 - null.fittedvalues) / d
    projection = np.linalg.solve(nuisance.T @ (w[:, None] * nuisance), nuisance.T @ (w * group))
    adjusted = group - nuisance @ projection
    expected_score = adjusted @ residual
    expected_variance = np.sum(w * adjusted**2)
    r = core().test_score(m, n, x, np.array([0.0, 0.0, 1.0]), rho)
    np.testing.assert_allclose(r["score"], expected_score, atol=1e-8)
    np.testing.assert_allclose(r["variance"], expected_variance, atol=1e-8)
    np.testing.assert_allclose(
        r["pvalue"], 2 * norm.sf(abs(expected_score) / np.sqrt(expected_variance)), rtol=1e-7
    )


def test_group_swap_preserves_p_and_reverses_response_effect():
    m = np.array([[16, 12, 17, 4, 5, 7]])
    n = np.array([[20, 18, 22, 19, 24, 20]])
    a = core().test_score(m, n, design(), np.array([0.0, 1.0]), 0.08)
    order = [3, 4, 5, 0, 1, 2]
    b = core().test_score(m[:, order], n[:, order], design(), np.array([0.0, 1.0]), 0.08)
    np.testing.assert_allclose(a["pvalue"], b["pvalue"], rtol=1e-7, atol=1e-12)
    np.testing.assert_allclose(a["effect"], -b["effect"], atol=1e-8)


def test_signed_and_scaled_treatment_contrast_preserves_observed_counterfactuals():
    m, n = np.array([[16, 12, 17, 4, 5, 7]]), np.full((1, 6), 20)
    positive = core().test_score(m, n, design(), [0.0, 1.0], 0.08)
    negative = core().test_score(m, n, design(), [0.0, -1.0], 0.08)
    scaled = core().test_score(m, n, design(), [0.0, 2.0], 0.08)
    np.testing.assert_allclose(positive["pvalue"], negative["pvalue"], rtol=1e-10)
    np.testing.assert_allclose(positive["pvalue"], scaled["pvalue"], rtol=1e-10)
    np.testing.assert_allclose(positive["effect"], -negative["effect"], atol=1e-8)
    np.testing.assert_allclose(positive["effect"], scaled["effect"], atol=1e-8)


def test_local_missing_coverage_masks_unestimable_contrast():
    n = np.array([[20] * 6, [0, 0, 0, 20, 20, 20]])
    m = np.array([[16, 12, 17, 4, 5, 7], [0, 0, 0, 4, 5, 7]])
    r = core().test_score(m, n, design(), np.array([0.0, 1.0]), 0.08)
    assert np.isfinite(r["pvalue"][0])
    assert np.isnan(r["pvalue"][1])
    assert not r["estimable"][1]


def test_zero_coverage_observation_equals_removing_that_row():
    x = design()
    m = np.array([[0, 12, 17, 4, 5, 7]])
    n = np.array([[0, 18, 22, 19, 24, 20]])
    a = core().test_score(m, n, x, np.array([0.0, 1.0]), 0.08)
    b = core().test_score(m[:, 1:], n[:, 1:], x[1:], np.array([0.0, 1.0]), 0.08)
    np.testing.assert_allclose(a["pvalue"], b["pvalue"], atol=1e-12)
    np.testing.assert_allclose(a["effect"], b["effect"], atol=1e-8)


def test_endpoint_groups_do_not_get_artificially_narrow_intervals():
    m = np.array([[0] * 6, [20] * 6, [20, 20, 20, 0, 0, 0]])
    n = np.full_like(m, 20)
    r = core().test_score(m, n, design(), np.array([0.0, 1.0]), 0.1)
    np.testing.assert_allclose(r["pvalue"][:2], [1.0, 1.0])
    assert r["pvalue"][2] < 0.01
    assert not r["interval_available"].any()
    assert np.isnan(r["ci_lo"]).all() and np.isnan(r["ci_hi"]).all()


@pytest.mark.parametrize(
    "m,n,rho",
    [([[2]], [[1]], 0.1), ([[0.5]], [[10]], 0.1), ([[0]], [[0]], 1.0), ([[0]], [[0]], -0.1)],
)
def test_invalid_count_or_dispersion_is_refused(m, n, rho):
    with pytest.raises(ValueError):
        core().fit_quasi_means(np.array(m), np.array(n), np.ones((1, 1)), rho)


def test_aliased_design_and_zero_contrast_are_refused():
    cs = core()
    m, n = np.array([[1, 2, 3, 1, 2, 3]]), np.full((1, 6), 10)
    with pytest.raises(ValueError, match="rank"):
        cs.fit_quasi_means(m, n, np.ones((6, 2)), 0.1)
    with pytest.raises(ValueError, match="contrast"):
        cs.test_score(m, n, design(), np.zeros(2), 0.1)


def test_high_depth_information_saturates_with_positive_rho():
    n = np.full((1, 6), 1000000, dtype=np.int64)
    m = np.array([[600000, 700000, 650000, 300000, 400000, 350000]])
    fit = core().fit_quasi_means(m, n, design(), 0.1)
    assert fit["effective_weights"].max() <= 10.0
    r = core().test_score(m, n, design(), np.array([0.0, 1.0]), 0.1)
    assert np.isfinite(r["pvalue"][0]) and r["pvalue"][0] > 0.001
