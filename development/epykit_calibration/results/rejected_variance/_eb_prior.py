"""Mean/depth-aware moderation of Pearson variance estimates.

The scaled-F log-moment fit subtracts residual sampling noise before estimating
prior heterogeneity. It is an approximation for count data, not a guarantee of
finite-sample calibration. No truth labels or external dispersion are used.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import digamma, polygamma


def fit_scaled_f(variance, residual_df):
    """Fit inverse-chi-square prior scale and df from noisy sample variances."""
    variance=np.asarray(variance,dtype=float)
    residual_df=np.asarray(residual_df,dtype=float)
    positive=variance[variance>0]
    floor=1e-5*float(np.median(positive)) if positive.size else 1e-5
    x=np.maximum(variance,floor)
    adjusted=np.log(x)+np.log(residual_df/2)-digamma(residual_df/2)
    mean=float(adjusted.mean())
    excess=float(np.var(adjusted,ddof=1)-np.mean(polygamma(1,residual_df/2)))
    if excess<=0:
        return float(x.mean()),float(residual_df.sum())
    shape=brentq(lambda z:float(polygamma(1,z))-excess,1e-8,1e8)
    prior_df=min(2*shape,float(residual_df.sum()))
    scale=float(np.exp(mean+digamma(prior_df/2)-np.log(prior_df/2)))
    return scale,prior_df


def moderate_pearson(pearson, residual_df, pooled_mu, depth, usable, minimum=100):
    """Moderate within mean/depth strata, retaining the sampling uncertainty.

    Fits use raw variance estimates, before their binomial lower bound. Sparse
    strata use the global fit. Sites with unestimable within-group variance use
    the fitted stratum prior with its own df, not chromosome-wide residual df.
    """
    pearson=np.asarray(pearson,dtype=float)
    df=np.maximum(np.asarray(residual_df,dtype=float),1.)
    raw=np.maximum(pearson/df,0)
    usable=np.asarray(usable,dtype=bool)&np.isfinite(raw)&np.isfinite(pooled_mu)&np.isfinite(depth)
    if usable.sum()<minimum:
        return np.maximum(raw,1),df.copy()
    global_scale,global_df=fit_scaled_f(raw[usable],df[usable])
    mu_labels=np.digitize(np.minimum(pooled_mu,1-pooled_mu),[.02,.05,.1,.2,.35])
    cuts=np.unique(np.quantile(depth[usable],[.25,.5,.75]))
    labels=mu_labels*4+np.digitize(depth,cuts)
    scale=np.full(len(raw),global_scale);prior_df=np.full(len(raw),global_df)
    for label in np.unique(labels):
        selected=labels==label;fit=selected&usable
        if fit.sum()>=minimum:
            s,d=fit_scaled_f(raw[fit],df[fit])
            scale[selected]=s;prior_df[selected]=d
    posterior=(pearson+prior_df*scale)/(df+prior_df)
    posterior=np.where(usable,posterior,scale)
    effective_df=np.where(usable,df+prior_df,prior_df)
    return np.maximum(posterior,1),np.maximum(effective_df,1)
