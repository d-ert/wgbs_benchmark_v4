# Statistical inference status

The corrected chain-merge caller adjusts over its complete coordinate-defined
interval family. This addresses search multiplicity only when each fixed
interval's raw p-value is valid. Its signed Stouffer statistic still assumes
independent CpG scores; correlated/smoothed data require additional calibration.

The legacy sliding-window and rule-segmentation callers select boundaries from
the same observations that supply their p-values. Their p/q columns are
exploratory ranking outputs. Selected-span BH and per-chromosome segment BH
are not demonstrated genome-wide region FDR control. The API records this
status in `dmr_params.inference_status`.

Complete-scan max-T randomization includes successful empty scans and retains
self/mirror assignments. Its adjusted p-values control complete-null FWER under
exchangeable assignments. Strong control under partial alternatives needs
additional assumptions. Failed scans abort inference; no second BH correction
is applied to these adjusted p-values.

Count-ratio regional FDR is an estimate using eligible permutation survivor
counts. It excludes self/mirror assignments, retains successful empty scans,
and aborts on failed scans. Outputs include requested/used permutation counts,
null-region counts, and the Monte Carlo standard error of the untruncated
set-level count ratio. This is not a standard error or confidence bound for
each selected q-value. Zero sampled null counts, including a zero sample
standard error, do not prove zero error. Record the design, assignment scheme,
Monte Carlo resolution and limitations with reported estimates.

Missing count coverage is not zero methylation. Smoothing averages alter
sampling variance and induce dependence; unsmoothed binomial/BB variances
cannot automatically be reused for smoothed pseudo-counts.

The development `bb_score` backend separates count sampling depth from
biological dispersion and supports nuisance regression through quasi-score
IRLS. Its mean-stratified count-mixture prior retains rho=0. Shared and separate
group dispersion variants record posterior standard deviations, iteration
convergence and low-support fallback. These quantities do not imply reference
degrees of freedom. Normal score tails and interior response intervals condition
on plug-in dispersion and do not integrate fitted-mean, dispersion or learned-
prior uncertainty. Boundary intervals are unavailable. General contrast effects
standardize the orthogonal nuisance design and compare the observed range along
the contrast projection; these counterfactuals may not describe actual subjects.

Initial estimated-dispersion development diagnostics found anti-conservative
small-sample tails, including when a learned mean trend understates heterogeneous
site dispersion. This backend remains experimental. A fixed-design count null
bootstrap can diagnose an analytical tail conditional on estimated generating
parameters; a few such diagnostics do not establish genome-wide FDR or resolve
model misspecification. Preserve the corrected workflow as the comparison
baseline until independent evaluation supports a new recommendation.

Relevant precedents include the [DSS general-design count model](https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html)
and [edgeR quasi-likelihood inference](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).
The present working model must earn its own calibration; results from those
implementations do not validate this score approximation.
