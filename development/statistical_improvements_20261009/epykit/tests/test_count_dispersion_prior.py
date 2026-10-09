"""Observed-count prior and estimated-rho contracts, with independent oracles."""

import importlib
import importlib.util

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import betabinom, binom


def prior_module():
    assert importlib.util.find_spec("epykit._dispersion_prior") is not None, (
        "Count prior is missing"
    )
    return importlib.import_module("epykit._dispersion_prior")


def counts(seed=710, sites=120):
    rng = np.random.default_rng(seed)
    n = rng.integers(10, 81, size=(sites, 10))
    latent = rng.beta(2.0, 5.0, size=n.shape)
    return rng.binomial(n, latent), n


def simple_prior():
    return dict(
        rho_grid=[0.0, 0.05, 0.2],
        pooled_rho_weights=[0.2, 0.5, 0.3],
        rho_weights_by_bin=[[0.2, 0.5, 0.3]] * 4,
        mean_bin_edges=[0.0, 0.05, 0.15, 0.3, 0.5],
        algorithm_revision="independent-fixture",
    )


def design():
    return np.column_stack([np.ones(10), np.r_[np.ones(5), np.zeros(5)]])


def test_bounded_prior_preserves_binomial_component_and_is_serializable():
    import json

    m, n = counts()
    prior = prior_module().fit_dispersion_prior(m, n, 5, max_training=80, seed=22)
    assert prior["rho_grid"][0] == 0.0
    assert prior["n_training_groups"] <= 80
    assert np.isfinite(prior["rho_weights_by_bin"]).all()
    np.testing.assert_allclose(np.sum(prior["rho_weights_by_bin"], axis=1), 1.0)
    assert prior["algorithm_revision"] and prior["training_digest"]
    assert "truth" not in prior
    json.dumps(prior)


def test_training_is_group_swap_row_order_and_seed_reproducible():
    m, n = counts()
    pm = prior_module()
    a = pm.fit_dispersion_prior(m, n, 5, max_training=80, seed=22)
    order = np.random.default_rng(91).permutation(len(m))
    b = pm.fit_dispersion_prior(
        m[order][:, np.r_[5:10, 0:5]], n[order][:, np.r_[5:10, 0:5]], 5, max_training=80, seed=22
    )
    assert a["training_digest"] == b["training_digest"]
    np.testing.assert_allclose(a["rho_weights_by_bin"], b["rho_weights_by_bin"], atol=1e-12)
    c = pm.fit_dispersion_prior(m, n, 5, max_training=80, seed=23)
    assert a["training_digest"] != c["training_digest"]


def test_unequal_group_sizes_preserve_count_prior_on_group_swap():
    m, n = counts(sites=80)
    m, n = m[:, :8], n[:, :8]
    pm = prior_module()
    a = pm.fit_dispersion_prior(m, n, 3, max_training=80, seed=22)
    order = np.r_[3:8, 0:3]
    b = pm.fit_dispersion_prior(m[:, order], n[:, order], 5, max_training=80, seed=22)
    assert a["training_digest"] == b["training_digest"]
    np.testing.assert_allclose(a["rho_weights_by_bin"], b["rho_weights_by_bin"], atol=1e-12)


def test_missing_coverage_contributes_no_likelihood_or_training_observation():
    m, n = counts(sites=80)
    m[:, 0], n[:, 0] = 0, 0
    pm = prior_module()
    a = pm.posterior_dispersion(m, n, np.full(m.shape, 0.3), simple_prior())
    b = pm.posterior_dispersion(m[:, 1:], n[:, 1:], np.full(m[:, 1:].shape, 0.3), simple_prior())
    np.testing.assert_allclose(a["weights"], b["weights"], atol=1e-12)
    # An entirely absent site's group must not enter the training bound.
    c = pm.fit_dispersion_prior(
        np.vstack([m, np.zeros((1, 10), int)]),
        np.vstack([n, np.zeros((1, 10), int)]),
        5,
        max_training=1000,
    )
    d = pm.fit_dispersion_prior(m, n, 5, max_training=1000)
    assert c["training_digest"] == d["training_digest"]


def test_conditional_posterior_matches_normalized_scipy_count_likelihood():
    pm = prior_module()
    m = np.array([[1, 8, 0], [3, 6, 12]])
    n = np.array([[10, 20, 0], [8, 20, 30]])
    mu = np.array([[0.1, 0.4, np.nan], [0.3, 0.3, 0.5]])
    p = simple_prior()
    actual = pm.posterior_dispersion(m, n, mu, p)
    likelihood = []
    for rho in p["rho_grid"]:
        if rho == 0:
            logp = binom.logpmf(m, n, mu)
        else:
            k = 1 / rho - 1
            logp = betabinom.logpmf(m, n, mu * k, (1 - mu) * k)
        likelihood.append(np.where(n > 0, logp, 0.0).sum(axis=1))
    logw = np.array(likelihood).T + np.log(p["pooled_rho_weights"])
    expected = np.exp(logw - logsumexp(logw, axis=1)[:, None])
    np.testing.assert_allclose(actual["weights"], expected, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(actual["rho_mean"], expected @ p["rho_grid"])
    np.testing.assert_allclose(
        actual["rho_variance"],
        expected @ (np.array(p["rho_grid"]) ** 2) - (expected @ p["rho_grid"]) ** 2,
    )


def test_uninformative_endpoint_posterior_retains_prior_and_no_coverage_status():
    m, n = np.zeros((2, 6), int), np.array([[20] * 6, [0] * 6])
    r = prior_module().posterior_dispersion(m, n, np.zeros((2, 6)), simple_prior())
    np.testing.assert_allclose(r["weights"], [simple_prior()["pooled_rho_weights"]] * 2)
    assert r["valid_counts"].tolist() == [6, 0]


def test_near_zero_rho_posterior_matches_exact_two_read_probabilities():
    pm = prior_module()
    m, n = np.array([[1, 0, 2]]), np.full((1, 3), 2)
    mu = np.array([[0.4, 0.6, 0.3]])
    p = simple_prior()
    p["rho_grid"] = [0.0, 1e-12, 0.2]
    products = []
    for rho in p["rho_grid"]:
        products.append(
            2 * 0.4 * 0.6 * (1 - rho) * (0.4 * (1 - 0.6 * (1 - rho))) * (0.3 * (0.3 + 0.7 * rho))
        )
    expected = np.array(products) * p["pooled_rho_weights"]
    expected /= expected.sum()
    actual = pm.posterior_dispersion(m, n, mu, p)
    np.testing.assert_allclose(actual["weights"][0], expected, atol=1e-12, rtol=1e-12)


def test_mean_dependent_training_recovers_dispersion_ordering_better_than_pooled():
    rng = np.random.default_rng(122)
    size = 500
    means = np.r_[np.full(size, 0.1), np.full(size, 0.4)]
    rho = np.r_[np.full(size, 0.008), np.full(size, 0.18)]
    k = 1 / rho - 1
    n = np.full((2 * size, 10), 60)
    m = rng.binomial(n, rng.beta((means * k)[:, None], ((1 - means) * k)[:, None], size=n.shape))
    p = prior_module().fit_dispersion_prior(m, n, 5, max_training=1600, seed=29)
    trend = np.array(p["rho_weights_by_bin"]) @ p["rho_grid"]
    pooled = np.array(p["pooled_rho_weights"]) @ p["rho_grid"]
    assert trend[1] < trend[3]
    assert abs(trend[1] - 0.008) + abs(trend[3] - 0.18) < abs(pooled - 0.008) + abs(pooled - 0.18)
    assert not p["bin_fallback"][1] and not p["bin_fallback"][3]


def test_general_design_prior_records_conditional_per_sample_mean_approximation():
    m, n = counts(sites=70)
    x = np.column_stack([design(), np.tile(np.arange(5), 2)])
    p = prior_module().fit_dispersion_prior(m, n, 5, design=x, max_training=100, seed=23)
    assert p["training_method"] == "conditional fitted per-sample means"
    assert "uncertainty" in p["limitations"]
    assert p["design_columns"] == 3


@pytest.mark.parametrize("mode", ["shared", "group"])
def test_estimated_dispersion_score_refits_and_preserves_swap(mode):
    prior_module()
    cs = importlib.import_module("epykit._count_score")
    assert hasattr(cs, "test_counts"), "Estimated-dispersion score wrapper is missing"
    m, n = counts(sites=20)
    a = cs.test_counts(
        m, n, design(), [0.0, 1.0], prior=simple_prior(), n_case=5, dispersion_mode=mode
    )
    order = np.r_[5:10, 0:5]
    b = cs.test_counts(
        m[:, order],
        n[:, order],
        design(),
        [0.0, 1.0],
        prior=simple_prior(),
        n_case=5,
        dispersion_mode=mode,
    )
    np.testing.assert_allclose(a["pvalue"], b["pvalue"], atol=1e-10, rtol=1e-7)
    np.testing.assert_allclose(a["effect"], -b["effect"], atol=1e-8)
    assert np.all(a["rho"] >= 0) and np.all(a["rho"] < 1)
    assert np.all(a["rho_uncertainty"] >= 0)
    assert np.all(a["dispersion_iterations"] <= 4)
    assert "conditional" in a["reference"] and a["experimental"]


def test_conditional_prior_mean_fits_are_bounded_before_training_selection(monkeypatch):
    pm = prior_module()
    cs = importlib.import_module("epykit._count_score")
    m, n = counts(sites=300)
    x = np.column_stack([design(), np.tile(np.arange(5), 2)])
    observed = []
    original = cs.fit_quasi_means

    def capture(m, n, *args, **kwargs):
        observed.append(len(m))
        return original(m, n, *args, **kwargs)

    monkeypatch.setattr(cs, "fit_quasi_means", capture)
    pm.fit_dispersion_prior(m, n, 5, design=x, max_training=40, seed=12)
    assert observed and max(observed) <= 40


def test_group_specific_fit_uses_explicit_shared_fallback_for_low_support():
    prior_module()
    cs = importlib.import_module("epykit._count_score")
    assert hasattr(cs, "test_counts"), "Estimated-dispersion score wrapper is missing"
    m, n = counts(sites=3)
    m[1, :4], n[1, :4] = 0, 0
    r = cs.test_counts(
        m, n, design(), [0.0, 1.0], prior=simple_prior(), n_case=5, dispersion_mode="group"
    )
    assert r["dispersion_fallback"][1]
    assert r["rho_case"][1] == pytest.approx(r["rho_control"][1])


@pytest.mark.parametrize("kwargs", [{"max_training": 0}, {"n_case": 0}, {"n_case": 10}])
def test_invalid_training_configuration_is_refused(kwargs):
    m, n = counts(sites=3)
    options = dict(n_case=5, max_training=80)
    options.update(kwargs)
    with pytest.raises(ValueError):
        prior_module().fit_dispersion_prior(m, n, **options)
