"""Experimental quasi-score inference with depth-aware biological variance.

These are working-model score/IRLS calculations, not BB likelihood fits.
Normal tails and interior delta-method intervals require calibration with
small biological sample sizes. Coverage is conditioned on, and N=0 is missing.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg, special, stats

from ._beta_binomial import _check_rho, _counts

ALGORITHM_REVISION = "bb-quasi-score-v2-conditional-eb-experimental"


def _inputs(m, n, design, rho):
    m, n = _counts(m, n, 2)
    x = np.asarray(design, dtype=float)
    if x.ndim != 2 or x.shape[0] != m.shape[1] or not np.isfinite(x).all():
        raise ValueError("Design must be a finite sample-by-coefficient matrix")
    if x.shape[1] and np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("Design is rank deficient")
    r = _check_rho(rho)
    if r.ndim == 1 and r.shape == (len(m),):
        r = r[:, None]
    try:
        r = np.broadcast_to(r, m.shape)
    except ValueError as exc:
        raise ValueError("rho must broadcast by site/sample") from exc
    return m.astype(float), n.astype(float), x, r


def _rank_mask(observed, x):
    if not x.shape[1]:
        return observed.any(axis=1)
    if observed.all():
        return np.ones(len(observed), dtype=bool)
    patterns, inverse = np.unique(observed, axis=0, return_inverse=True)
    valid = np.array(
        [
            mask.sum() >= x.shape[1] and np.linalg.matrix_rank(x[mask]) == x.shape[1]
            for mask in patterns
        ]
    )
    return valid[inverse]


def _solve_systems(information, rhs, valid):
    """Batched solve, isolating numerical singularities instead of adding ridge."""
    result = np.full_like(rhs, np.nan, dtype=float)
    usable = valid.copy()
    indices = np.flatnonzero(valid)
    if not len(indices):
        return result, usable
    try:
        result[indices] = np.linalg.solve(information[indices], rhs[indices])
    except np.linalg.LinAlgError:
        for i in indices:
            try:
                result[i] = np.linalg.solve(information[i], rhs[i])
            except np.linalg.LinAlgError:
                usable[i] = False
    usable &= np.isfinite(result).all(axis=(1, 2))
    return result, usable


def _information(x, weights):
    return np.einsum("sp,bs,sq->bpq", x, weights, x, optimize=True)


def _objective(m, n, x, beta, inflation):
    eta = beta @ x.T
    return ((m * eta - n * np.logaddexp(0.0, eta)) / inflation).sum(axis=1)


def _fit(m, n, x, rho, max_iter, tol):
    sites, samples = m.shape
    p = x.shape[1]
    observed = n > 0
    inflation = 1 + (n - 1) * rho
    effective = np.divide(n, inflation)
    valid = _rank_mask(observed, x)
    beta = np.zeros((sites, p))
    converged = np.zeros(sites, dtype=bool)
    iterations = np.zeros(sites, dtype=np.int32)
    if not p:
        return dict(
            coefficients=beta,
            mu=np.full((sites, samples), 0.5),
            information=np.empty((sites, 0, 0)),
            effective_weights=effective,
            valid_counts=observed.sum(axis=1),
            converged=valid,
            estimable=valid,
            boundary=np.zeros(sites, dtype=bool),
            iterations=iterations,
        )

    gamma = special.logit((m + 0.5) / (n + 1))
    initial = _information(x, effective)
    initial_rhs = ((effective * gamma) @ x)[..., None]
    solved, valid = _solve_systems(initial, initial_rhs, valid)
    beta[valid] = solved[valid, :, 0]
    # A constant in the design span permits common exact endpoint means.
    constant_fit = x @ np.linalg.lstsq(x, np.ones(samples), rcond=None)[0]
    has_constant = np.allclose(constant_fit, 1.0, atol=1e-10, rtol=0.0)
    total_m, total_n = m.sum(axis=1), n.sum(axis=1)
    endpoint = valid & has_constant & ((total_m == 0) | (total_m == total_n))
    converged[endpoint] = True
    active = valid & ~endpoint

    for _ in range(max_iter):
        indices = np.flatnonzero(active)
        if not len(indices):
            break
        mu = np.clip(special.expit(beta @ x.T), 1e-12, 1 - 1e-12)
        info = _information(x, effective * mu * (1 - mu))
        gradient = ((m - n * mu) / inflation) @ x
        delta, solvable = _solve_systems(info, gradient[..., None], active)
        active &= solvable
        indices = np.flatnonzero(active)
        if not len(indices):
            break
        step = delta[indices, :, 0]
        max_step = np.max(np.abs(step), axis=1)
        done = (max_step <= tol) | (np.max(np.abs(gradient[indices]), axis=1) <= tol)
        converged[indices[done]] = True
        active[indices[done]] = False
        indices = indices[~done]
        step = step[~done]
        if not len(indices):
            continue
        step *= np.minimum(1.0, 5 / np.maximum(np.max(np.abs(step), axis=1), 1e-20))[:, None]
        old_value = _objective(m[indices], n[indices], x, beta[indices], inflation[indices])
        scale = np.ones(len(indices))
        proposal = beta[indices] + step
        for _ in range(20):
            value = _objective(m[indices], n[indices], x, proposal, inflation[indices])
            bad = value < old_value - 1e-10
            if not bad.any():
                break
            scale[bad] *= 0.5
            proposal = beta[indices] + scale[:, None] * step
        beta[indices] = proposal
        iterations[indices] += 1

    mu = special.expit(beta @ x.T)
    mu[endpoint] = (total_m[endpoint] > 0)[:, None].astype(float)
    info = _information(x, effective * mu * (1 - mu))
    boundary = endpoint | (((mu < 1e-6) | (mu > 1 - 1e-6)) & observed).any(axis=1)
    mu[~valid] = np.nan
    beta[~valid | endpoint] = np.nan
    return dict(
        coefficients=beta,
        mu=mu,
        information=info,
        effective_weights=effective,
        valid_counts=observed.sum(axis=1),
        converged=converged,
        estimable=valid,
        boundary=boundary,
        iterations=iterations,
    )


def fit_quasi_means(m, n, design, rho, *, max_iter=40, tol=1e-8):
    """Fit mean regression under fixed biological rho using quasi-score IRLS.

    `effective_weights` are N/[1+(N-1)rho], separate from the information
    weights that also include mu(1-mu). Boundary coefficients are unavailable
    for identical endpoint observations even though their point means exist.
    """
    if max_iter < 1 or not np.isfinite(tol) or tol <= 0:
        raise ValueError("Require positive iteration budget and tolerance")
    return _fit(*_inputs(m, n, design, rho), max_iter, tol)


def test_score(m, n, design, contrast, rho, *, intervals=True):
    """Efficient one-dimensional null score with nuisance terms projected out.

    The contrast is a coefficient vector c testing c beta=0. Response effects
    compare standardized counterfactuals at the observed minimum/maximum of
    Xv, where v=c/(c'c), holding the orthogonal nuisance design fixed. A binary
    treatment coefficient compares its observed 0/1 levels, with sign reversed
    for a negative contrast and unchanged by positive rescaling. General
    contrasts use the observed range, which may not represent actual subjects.
    Intervals are interior conditional
    delta-method approximations, not an inversion of the null score test.
    """
    mf, nf, x, r = _inputs(m, n, design, rho)
    c = np.asarray(contrast, dtype=float)
    if c.shape != (x.shape[1],) or not np.isfinite(c).all() or not np.any(c):
        raise ValueError("Require a finite nonzero one-dimensional contrast")
    basis = linalg.null_space(c[None, :])
    nuisance = x @ basis
    direction = x @ (c / (c @ c))
    null = _fit(mf, nf, nuisance, r, 40, 1e-8)
    full = _fit(mf, nf, x, r, 40, 1e-8)
    valid = full["estimable"] & null["estimable"] & (null["valid_counts"] > x.shape[1])
    endpoint = (
        valid
        & null["boundary"]
        & full["boundary"]
        & ((mf.sum(axis=1) == 0) | (mf.sum(axis=1) == nf.sum(axis=1)))
    )
    regular = valid & null["converged"] & ~endpoint
    inflation = 1 + (nf - 1) * r
    residual = (mf - nf * null["mu"]) / inflation
    weights = full["effective_weights"] * null["mu"] * (1 - null["mu"])
    score = residual @ direction
    variance = weights @ (direction * direction)
    if nuisance.shape[1]:
        cross = (weights * direction) @ nuisance
        solved, regular = _solve_systems(null["information"], cross[..., None], regular)
        score -= np.einsum("bp,bp->b", solved[:, :, 0], residual @ nuisance)
        variance -= np.einsum("bp,bp->b", cross, solved[:, :, 0])
    regular &= (variance > 0) & np.isfinite(variance) & np.isfinite(score)
    statistic = np.full(len(mf), np.nan)
    pvalue = np.full(len(mf), np.nan)
    statistic[regular] = score[regular] ** 2 / variance[regular]
    pvalue[regular] = np.exp(np.log(2.0) + stats.norm.logsf(np.sqrt(statistic[regular])))
    pvalue[endpoint], statistic[endpoint] = 1.0, 0.0
    score[~regular & ~endpoint], variance[~regular & ~endpoint] = np.nan, np.nan
    score[endpoint], variance[endpoint] = 0.0, 0.0

    low, high = float(direction.min()), float(direction.max())
    x0 = x - np.outer(direction - low, c)
    x1 = x0 + (high - low) * c
    mean0 = special.expit(full["coefficients"] @ x0.T)
    mean1 = special.expit(full["coefficients"] @ x1.T)
    effect = (mean1 - mean0).mean(axis=1)
    effect[endpoint] = 0.0
    lo, hi = np.full(len(mf), np.nan), np.full(len(mf), np.nan)
    interval_ok = valid & full["converged"] & ~full["boundary"]
    if intervals:
        identity = np.broadcast_to(np.eye(x.shape[1]), full["information"].shape)
        covariance, interval_ok = _solve_systems(full["information"], identity, interval_ok)
        gradient = ((mean1 * (1 - mean1)) @ x1 - (mean0 * (1 - mean0)) @ x0) / x.shape[0]
        effect_var = np.einsum("bp,bpq,bq->b", gradient, covariance, gradient)
        interval_ok &= (effect_var > 0) & np.isfinite(effect_var)
        half = stats.norm.ppf(0.975) * np.sqrt(effect_var[interval_ok])
        lo[interval_ok] = np.maximum(-1.0, effect[interval_ok] - half)
        hi[interval_ok] = np.minimum(1.0, effect[interval_ok] + half)
    else:
        interval_ok[:] = False
    return dict(
        score=score,
        variance=variance,
        statistic=statistic,
        pvalue=pvalue,
        effect=effect,
        ci_lo=lo,
        ci_hi=hi,
        interval_available=interval_ok,
        estimable=regular | endpoint,
        converged=null["converged"],
        full_converged=full["converged"],
        mu_null=null["mu"],
        mu_full=full["mu"],
        coefficients=full["coefficients"],
        valid_counts=null["valid_counts"],
        effective_weights=full["effective_weights"],
        effect_estimand="standardized response contrast over observed design-projection range",
        reference="normal approximation conditional on supplied rho",
    )


def test_counts(m, n, design, contrast, *, prior, n_case, dispersion_mode="shared", intervals=True):
    """Experimental plug-in EB biological-dispersion score for count regressions.

    At most four mean/posterior updates use a fixed learned prior. A shared rho
    is the default; separate group dispersions require >=3 covered replicates
    per group, otherwise use the explicit shared fallback. Posterior standard
    deviations diagnose uncertainty but do not adjust normal reference tails.
    The fitted-mean, dispersion and prior uncertainty are not integrated into
    the reference law. Development null bootstrap is required before use as a
    recommended small-sample inference engine.
    """
    from ._dispersion_prior import posterior_dispersion

    m, n = _counts(m, n, 2)
    if not isinstance(n_case, (int, np.integer)) or not 0 < n_case < m.shape[1]:
        raise ValueError("n_case must split two nonempty groups")
    if dispersion_mode not in ("shared", "group"):
        raise ValueError("dispersion_mode must be shared or group")
    x = np.asarray(design, dtype=float)
    grid = _check_rho(prior["rho_grid"])
    weights = np.asarray(prior["pooled_rho_weights"], dtype=float)
    if (
        weights.shape != grid.shape
        or not np.isfinite(weights).all()
        or np.any(weights < 0)
        or weights.sum() <= 0
    ):
        raise ValueError("Invalid pooled dispersion prior")
    initial = float((weights / weights.sum()) @ grid)
    rho = np.full(m.shape, initial)
    iterations = np.zeros(len(m), dtype=np.int32)
    stable = np.zeros(len(m), dtype=bool)
    group_support = ((n[:, :n_case] > 0).sum(axis=1) >= 3) & ((n[:, n_case:] > 0).sum(axis=1) >= 3)
    fallback = ~group_support if dispersion_mode == "group" else np.zeros(len(m), bool)
    posterior_sd = np.full(m.shape, np.nan)
    last_mean = np.full(m.shape, np.nan)
    for iteration in range(4):
        fit = fit_quasi_means(m, n, x, rho)
        usable = fit["estimable"] & fit["converged"]
        fallback |= ~usable
        pm, pn = np.where(usable[:, None], m, 0), np.where(usable[:, None], n, 0)
        mu = np.where(usable[:, None], fit["mu"], 0.5)
        shared = posterior_dispersion(pm, pn, mu, prior)
        updated = np.broadcast_to(shared["rho_mean"][:, None], m.shape).copy()
        uncertainty = np.broadcast_to(np.sqrt(shared["rho_variance"])[:, None], m.shape).copy()
        if dispersion_mode == "group":
            groups = posterior_dispersion(
                pm,
                pn,
                mu,
                prior,
                group_labels=np.r_[np.zeros(n_case, int), np.ones(m.shape[1] - n_case, int)],
            )
            for label, (lo, hi) in enumerate(((0, n_case), (n_case, m.shape[1]))):
                rows = group_support & usable
                updated[rows, lo:hi] = groups["rho_mean"][rows, label, None]
                uncertainty[rows, lo:hi] = np.sqrt(groups["rho_variance"][rows, label, None])
        change = np.max(np.abs(updated - rho), axis=1)
        mean_change = np.max(np.abs(mu - last_mean), axis=1)
        now_stable = (change <= 1e-4) & (mean_change <= 1e-5) & usable
        # Retain completed rows exactly, so iteration budget is explicit.
        active = ~stable
        rho[active], posterior_sd[active] = updated[active], uncertainty[active]
        iterations[active] = iteration + 1
        last_mean = mu
        stable |= now_stable
        if stable.all():
            break
    result = test_score(m, n, x, contrast, rho, intervals=intervals)
    result.update(
        rho=rho.mean(axis=1),
        rho_case=rho[:, :n_case].mean(axis=1),
        rho_control=rho[:, n_case:].mean(axis=1),
        rho_uncertainty=posterior_sd.mean(axis=1),
        rho_case_uncertainty=posterior_sd[:, :n_case].mean(axis=1),
        rho_control_uncertainty=posterior_sd[:, n_case:].mean(axis=1),
        dispersion_iterations=iterations,
        dispersion_converged=stable,
        dispersion_fallback=fallback,
        dispersion_mode=dispersion_mode,
        experimental=True,
        prior_revision=prior.get("algorithm_revision", "unknown"),
        algorithm_revision=ALGORITHM_REVISION,
        reference="normal approximation conditional on plug-in EB rho; "
        "fitted-mean, dispersion and prior uncertainty not integrated",
    )
    return result
