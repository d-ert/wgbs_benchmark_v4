"""Bounded, observed-count EB priors for biological beta-binomial rho.

The empirical mean strata and conditional regression variant are working
approximations. Neither removes fitted-mean or learned-prior uncertainty.
The zero-rho atom is retained; sequencing depth enters the likelihood.
"""

from __future__ import annotations

import hashlib

import numpy as np
from numba import njit
from scipy import special

from ._beta_binomial import _check_rho, _counts, _fit_mixture, _mixture_loglik, _relative_logpmf

ALGORITHM_REVISION = "mean-stratified-count-mixture-v2-bounded-experimental"
MEAN_BIN_EDGES = np.array([0.0, 0.05, 0.15, 0.30, 0.50])
MEAN_GRID = (1 - np.cos(np.linspace(0, np.pi, 41))) / 2
RHO_GRID = np.r_[0.0, np.geomspace(1e-4, 0.8, 25)]


@njit(cache=True)
def _small_rho_relative_loglik(m, n, mu, rho):
    result = np.zeros(len(m))
    for i in range(len(m)):
        for j in range(m.shape[1]):
            if n[i, j] > 0:
                result[i] += _relative_logpmf(m[i, j], n[i, j], mu[i, j], rho)
    return result


def _bin_index(means, edges):
    folded = np.minimum(means, 1 - means)
    return np.clip(np.searchsorted(edges, folded, side="right") - 1, 0, len(edges) - 2)


def _conditional_loglik(m, n, mu, rho_grid):
    """Normalized count likelihood at fixed per-sample means, bounded by rows.

    Smallest positive training rho is 1e-4, keeping concentration manageable.
    Endpoint means are degenerate counts. Missing observations add zero.
    """
    covered = n > 0
    mu = np.where(covered, mu, 0.5)
    choose = special.gammaln(n + 1) - special.gammaln(m + 1) - special.gammaln(n - m + 1)
    result = np.empty((len(m), len(rho_grid)))
    binomial = choose + special.xlogy(m, mu) + special.xlog1py(n - m, -mu)
    interior = (mu > 0) & (mu < 1)
    for i, rho in enumerate(rho_grid):
        if rho == 0:
            logp = binomial
        elif rho < 1e-5:
            result[:, i] = _small_rho_relative_loglik(m, n, mu, rho) + choose.sum(axis=1)
            continue
        else:
            k = (1 - rho) / rho
            a, b = np.where(interior, mu, 0.5) * k, np.where(interior, 1 - mu, 0.5) * k
            logp = choose + special.betaln(m + a, n - m + b) - special.betaln(a, b)
            logp = np.where(interior, logp, binomial)
        result[:, i] = np.where(covered, logp, 0.0).sum(axis=1)
    return result


def _fit_rho_mixture(log_likelihood, initial=None, max_iter=400, tol=1e-5):
    shift = log_likelihood.max(axis=1)
    if not np.isfinite(shift).all():
        raise ValueError("Dispersion grid gives zero likelihood to observed counts")
    likelihood = np.exp(log_likelihood - shift[:, None])
    weights = (
        np.full(likelihood.shape[1], 1 / likelihood.shape[1])
        if initial is None
        else np.array(initial, dtype=float, copy=True)
    )
    previous, converged = -np.inf, False
    for iteration in range(max_iter):
        weighted = likelihood * weights
        marginal = weighted.sum(axis=1)
        loglik = float((np.log(marginal) + shift).sum())
        if abs(loglik - previous) <= tol * len(likelihood):
            converged = True
            break
        previous = loglik
        weights = np.maximum((weighted / marginal[:, None]).mean(axis=0), 1e-300)
        weights /= weights.sum()
    return dict(rho_weights=weights, loglik=loglik, converged=converged, iterations=iteration + 1)


def _bounded_groups(m, n, n_case, mu, max_training, seed):
    width = max(n_case, m.shape[1] - n_case)
    groups_m, groups_n, groups_mu = [], [], []
    for lo, hi in ((0, n_case), (n_case, m.shape[1])):
        gm, gn = (
            np.zeros((len(m), width), dtype=np.int64),
            np.zeros((len(m), width), dtype=np.int64),
        )
        gm[:, : hi - lo], gn[:, : hi - lo] = m[:, lo:hi], n[:, lo:hi]
        order = np.lexsort((gm, gn), axis=1)
        gm, gn = np.take_along_axis(gm, order, 1), np.take_along_axis(gn, order, 1)
        groups_m.append(gm)
        groups_n.append(gn)
        if mu is not None:
            gu = np.full((len(m), width), 0.5)
            gu[:, : hi - lo] = mu[:, lo:hi]
            groups_mu.append(np.take_along_axis(gu, order, 1))
    gm, gn = np.vstack(groups_m), np.vstack(groups_n)
    gu = np.vstack(groups_mu) if mu is not None else None
    valid = (gn > 0).sum(axis=1) >= 2
    if gu is not None:
        valid &= np.isfinite(np.where(gn > 0, gu, 0.0)).all(axis=1)
    gm, gn = gm[valid], gn[valid]
    gu = gu[valid] if gu is not None else None
    if not len(gm):
        raise ValueError(
            "Dispersion training requires two covered replicates in at least one group"
        )
    # Content-based seeded priorities preserve row, group and within-group swaps.
    prefix = str(int(seed)).encode()
    fingerprints = []
    for i in range(len(gm)):
        h = hashlib.blake2b(digest_size=16)
        h.update(prefix)
        h.update(gm[i].astype("<i8").tobytes())
        h.update(gn[i].astype("<i8").tobytes())
        if gu is not None:
            h.update(np.round(np.where(gn[i] > 0, gu[i], 0.5), 12).astype("<f8").tobytes())
        fingerprints.append(h.hexdigest())
    indices = np.argsort(fingerprints, kind="stable")[:max_training]
    digest = hashlib.sha256("".join(fingerprints[i] for i in indices).encode()).hexdigest()
    return (
        gm[indices],
        gn[indices],
        gu[indices] if gu is not None else None,
        digest,
        int(valid.sum()),
    )


def fit_dispersion_prior(m, n, n_case, *, design=None, max_training=4096, seed=20261009):
    """Learn pooled/mean-stratified mixing weights from observed counts only.

    Each training unit is a replicate group with at least two covered samples.
    Sparse strata (<64 non-endpoint groups) use the pooled prior. Sampling uses
    seeded content hashes so group and row permutations do not change it.
    """
    m, n = _counts(m, n, 2)
    if not isinstance(n_case, (int, np.integer)) or not 0 < n_case < m.shape[1]:
        raise ValueError("n_case must split two nonempty groups")
    if not isinstance(max_training, (int, np.integer)) or max_training < 1:
        raise ValueError("max_training must be positive")
    available = ((n[:, :n_case] > 0).sum(axis=1) >= 2) | ((n[:, n_case:] > 0).sum(axis=1) >= 2)
    indices = np.flatnonzero(available)
    if not len(indices):
        raise ValueError(
            "Dispersion training requires two covered replicates in at least one group"
        )
    if len(indices) > max_training:
        # Bound sites before conditional regression or group-matrix allocation.
        # Canonical group order makes priorities invariant to swapping groups.
        width = max(n_case, m.shape[1] - n_case)
        prefix = str(int(seed)).encode()
        priorities = []
        for i in indices:
            payloads = []
            for lo, hi in ((0, n_case), (n_case, m.shape[1])):
                row = np.zeros((2, width), dtype="<i8")
                row[0, : hi - lo], row[1, : hi - lo] = m[i, lo:hi], n[i, lo:hi]
                payloads.append(row.tobytes())
            priorities.append(
                hashlib.blake2b(prefix + b"".join(sorted(payloads)), digest_size=16).hexdigest()
            )
        indices = indices[np.argsort(priorities, kind="stable")[:max_training]]
    m, n = m[indices], n[indices]
    mu = None
    if design is not None:
        from ._count_score import fit_quasi_means

        x = np.asarray(design)
        if x.ndim != 2 or x.shape[0] != m.shape[1]:
            raise ValueError("Design must be sample-aligned")
        mean_block = min(
            max_training, 32768, max(1, 134217728 // (8 * (8 * x.shape[1] ** 2 + 12 * m.shape[1])))
        )
        mu = np.empty(m.shape, dtype=float)
        for lo in range(0, len(m), mean_block):
            hi = lo + mean_block
            fit = fit_quasi_means(m[lo:hi], n[lo:hi], design, 0.05)
            mu[lo:hi] = fit["mu"]
            mu[lo:hi][~fit["converged"]] = np.nan
    gm, gn, gu, digest, n_available = _bounded_groups(m, n, n_case, mu, max_training, seed)
    if gu is None:
        means = gm.sum(axis=1) / gn.sum(axis=1)
        loglik = _mixture_loglik(gm, gn, MEAN_GRID, RHO_GRID)
        initial_mean = np.full(len(MEAN_GRID), 1 / len(MEAN_GRID))

        def fit_subset(mask, initial=None):
            starts = [np.full(len(RHO_GRID), 1 / len(RHO_GRID))]
            if initial is not None:
                starts.append(initial)
            else:
                warm = np.exp(-0.5 * ((np.log(np.maximum(RHO_GRID, 1e-5)) - np.log(0.05)) / 2) ** 2)
                starts.append(warm / warm.sum())
            fits = [_fit_mixture(loglik[mask], initial_mean, s, 1e-5, 400) for s in starts]
            return max(fits, key=lambda f: f["loglik"])

        method = "integrated replicate-group mean count mixture"
    else:
        means = np.where(gn > 0, gu, 0.0).sum(axis=1) / (gn > 0).sum(axis=1)
        loglik = _conditional_loglik(gm, gn, gu, RHO_GRID)

        def fit_subset(mask, initial=None):
            return _fit_rho_mixture(loglik[mask], initial)

        method = "conditional fitted per-sample means"
    pooled = fit_subset(np.ones(len(gm), dtype=bool))
    bin_indices = _bin_index(means, MEAN_BIN_EDGES)
    informative = (gm.sum(axis=1) > 0) & (gm.sum(axis=1) < gn.sum(axis=1))
    weights, fallback, counts, convergence = [], [], [], []
    for i in range(len(MEAN_BIN_EDGES) - 1):
        mask = bin_indices == i
        enough = int((mask & informative).sum()) >= 64
        fitted = fit_subset(mask, pooled["rho_weights"]) if enough else pooled
        weights.append(fitted["rho_weights"].tolist())
        fallback.append(not enough)
        counts.append(int((mask & informative).sum()))
        convergence.append(bool(fitted["converged"]))
    return dict(
        rho_grid=RHO_GRID.tolist(),
        mean_grid=MEAN_GRID.tolist(),
        mean_bin_edges=MEAN_BIN_EDGES.tolist(),
        rho_weights_by_bin=weights,
        pooled_rho_weights=pooled["rho_weights"].tolist(),
        bin_fallback=fallback,
        bin_informative_groups=counts,
        bin_converged=convergence,
        pooled_converged=bool(pooled["converged"]),
        n_training_groups=len(gm),
        n_available_groups=n_available,
        n_training_sites=len(m),
        n_available_sites=int(available.sum()),
        training_digest=digest,
        seed=int(seed),
        max_training=int(max_training),
        training_method=method,
        design_columns=None if design is None else int(np.asarray(design).shape[1]),
        algorithm_revision=ALGORITHM_REVISION,
        limitations="Empirical observed-mean strata; approximate within-stratum mean/rho "
        "independence; fitted regression mean and learned-prior uncertainty "
        "are not fully integrated. Analytical score reference is experimental.",
    )


def posterior_dispersion(m, n, mu, prior, *, group_labels=None):
    """Conditional count-likelihood posterior rho moments, retaining rho=0.

    Group labels request separate posteriors. This does not add a reference df.
    No observations retain the appropriate prior and report valid_counts=0.
    """
    m, n = _counts(m, n, 2)
    mu = np.asarray(mu, dtype=float)
    try:
        mu = np.broadcast_to(mu, m.shape)
    except ValueError as exc:
        raise ValueError("Means must broadcast to site-by-sample counts") from exc
    covered = n > 0
    if np.any(covered & (~np.isfinite(mu) | (mu < 0) | (mu > 1))):
        raise ValueError("Covered observations require finite means in [0,1]")
    if group_labels is not None:
        labels = np.asarray(group_labels)
        if labels.shape != (m.shape[1],):
            raise ValueError("group_labels must label each sample")
        unique = np.unique(labels)
        fits = [
            posterior_dispersion(m[:, labels == v], n[:, labels == v], mu[:, labels == v], prior)
            for v in unique
        ]
        return dict(
            {
                key: np.stack([f[key] for f in fits], axis=1)
                for key in ("weights", "rho_mean", "rho_variance", "valid_counts", "bin_index")
            },
            group_labels=unique,
        )
    grid = _check_rho(prior["rho_grid"])
    if grid.ndim != 1 or not len(grid) or np.any(np.diff(grid) <= 0):
        raise ValueError("Prior rho grid must be strictly increasing")
    edges = np.asarray(prior["mean_bin_edges"], dtype=float)
    weights = np.asarray(prior["rho_weights_by_bin"], dtype=float)
    if (
        edges.ndim != 1
        or len(edges) < 2
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
        or edges[0] != 0
        or edges[-1] != 0.5
        or weights.shape != (len(edges) - 1, len(grid))
        or not np.isfinite(weights).all()
        or np.any(weights < 0)
        or np.any(weights.sum(axis=1) <= 0)
    ):
        raise ValueError("Invalid mean-stratified dispersion prior")
    counts = covered.sum(axis=1)
    means = np.divide(
        np.where(covered, mu, 0.0).sum(axis=1), counts, out=np.full(len(m), 0.5), where=counts > 0
    )
    indices = _bin_index(means, edges)
    selected = weights[indices] / weights[indices].sum(axis=1)[:, None]
    with np.errstate(divide="ignore"):
        logw = _conditional_loglik(m, n, mu, grid) + np.log(selected)
    normalizer = special.logsumexp(logw, axis=1)
    if not np.isfinite(normalizer).all():
        raise ValueError("Conditional fitted means give zero count likelihood")
    posterior = np.exp(logw - normalizer[:, None])
    mean = posterior @ grid
    variance = np.maximum(0.0, posterior @ (grid * grid) - mean * mean)
    return dict(
        weights=posterior,
        rho_mean=mean,
        rho_variance=variance,
        valid_counts=counts,
        bin_index=indices,
    )
