# Empirical Bayes audit — 7 October 2026

Audited package commit: b5f11353602ffd017b2f7192a24d6b8b9df56d40.
The running benchmark source remains frozen; this audit changes no fit.

## Legacy lr engine, dispersion=eb

`dmc.py:_score_finalize` estimates per-site Pearson dispersion, clips it at
one, fits inverse-gamma shape from the observed mean and variance, and
computes a weighted average with the chromosome-pooled dispersion.

The prior fit does not separate sampling variation of estimated dispersions
from actual variation between latent dispersions. Clipping and selecting
only sites with non-boundary means further alter the fitted distribution.
Consequently this is data-adaptive, inverse-gamma-inspired shrinkage; it is
not a demonstrated likelihood fit of the stated hierarchical model.

There is also a mismatch between the stated posterior target and formula.
For IG(a,b), prior mean m=b/(a-1), and a scaled-chi-square sampling model
with site dispersion s and residual df d, the posterior variance mean is:

    (d*s + 2*(a-1)*m) / (d + 2*(a-1))

The implementation instead weights its chromosome anchor by 2*a. With
a=3, b=4, d=8, s=1.5, the documented-model posterior mean is 1.666667;
the implementation's expression with the prior mean as anchor is 1.714286.
The actual chromosome anchor may differ from the moment-fitted prior mean.

Important qualification: weighting by 2*a can be legitimate for the inverse
posterior mean of precision when the anchor is the prior scale b/a. It is
therefore not inherently an invalid weight. However, that is a different
target/scale from the source's posterior-mean and moment-fitting explanation.
Smyth's 2004 paper explicitly corrects the posterior-variance interpretation
in its 2009 erratum. The model, estimator target and hyperparameter fit must
be made consistent before claiming a faithful conjugate EB implementation.

Adding the learned weight to reference degrees of freedom also needs a
sampling-distribution derivation for the count/quasi-likelihood statistic;
Gaussian moderated-t theory alone does not establish it here. These code
issues do not isolate the cause of all observed false CpG calls.

## New beta_binomial engine

The prior is first fitted to actual replicate-count likelihoods using a
factorized mean/rho mixture. Log-rho moments then define an approximate
continuous normal prior. The per-site penalized count likelihood finds a
joint/profile MAP in log-rho coordinates and group means. This is an
empirical-Bayes shrinkage procedure; the numerical likelihood and objective
have independent oracle tests.

Its limitations remain material: grid approximation, moment projection,
mean/rho independence, bounded per-chromosome prior training, shared rho
between groups, and weak prior identifiability. Inference holds fitted rho
fixed and uses an approximate F reference rather than integrating posterior
dispersion or hyperparameter uncertainty. Implementation verification is
not proof of calibrated FDR or useful power in WGBS data.

The separate optional `shrink_meth_diff` helper implements normal-prior
effect shrinkage using estimated standard errors. It is not used by the
current benchmark and cannot repair dispersion/p-value calibration.

Additional source check: the `se_from='coef_se'` branch takes a GLM
standard error on the logit-coefficient scale, but continues to shrink
`meth_diff` on the methylation-proportion scale. The documentation's claim
that this branch instead shrinks a logit-scale effect does not match the
code. Those scales must be aligned before using that branch. The `ci`
branch's conversion by 2*1.96 assumes a symmetric normal/Wald interval;
it is not generally valid for t/F or profile-likelihood intervals.

Primary reference: Gordon Smyth (2004), with 2009 erratum:
https://www.math.tau.ac.il/~yekutiel/eBayes/smyth_2004.pdf
