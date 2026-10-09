# Epykit statistical improvements and comparative validation

Status: approved by the user in “i approve /goal resume”; implementation-plan review is next.

## Intent and success

The user requests implementation of the improvements identified in the statistical review, using Superpowers, followed by simulation and real-data tests to determine whether epykit improves. The deliverable is working package code, regression tests, reproducible benchmark commands, source snapshots, and an evidence-based before/after comparison. Success is trustworthy additional detection with practical runtime and memory, not a larger count at the same nominal q-value.

The investigation is documented in [the review](/scratch/wgbs_benchmark_v4/analysis/epykit_statistical_review_20261009/review.md). The corrected chain-merge workflow is the statistical comparison baseline. The original engine and BB-F remain additional comparators. New methods remain explicit experimental choices until their calibration is supported. A negative benchmark result must be reported and must not be hidden by promoting a weaker method or changing truth definitions.

## Approaches and choice

1. Numerical patches alone are inexpensive but do not address the principal calibration and region-power limitations.
2. **Recommended: a fast count-score engine plus sample-aware regional inference**, with numerical, resampling and memory fixes. This preserves Parquet I/O and offers a useful comparison before a large rewrite.
3. A joint Bayesian spatial/change-point backend is a later research alternative. It is not necessary to implement the present improvements, and its computational/prior costs are not yet justified by comparative evidence. Kinship mixed models are conditional future work, not an established explanation of the current failures.

The specification follows approach 2. It includes the full simulation/real-data evaluation, not merely a small diagnostic demonstration.

## Source isolation and compatibility

Create `development/statistical_improvements_20261009/epykit/` as an isolated checkout or source copy of `/scratch/epykit-calibration-fresh_20261007`. That source contains the complete-family correction and experimental BB model; `/scratch/epykit3` still matches the original relevant files. Keep a manifest of original hashes and a distributable patch against this recorded base.

The workspace permits writes under `/scratch/wgbs_benchmark_v4` and `/tmp`. Git metadata in the present workspace and the external worktree are read-only; use a writable local checkout or source-copy fallback rather than modify external worktrees or request broader access. Frozen result snapshots and source datasets remain intact.

Preserve Python >=3.10 compatibility, the package's Windows-supported core, the canonical DMC/DMR schemas, logging conventions, deprecated surfaces and existing defaults during experimental comparison. Add any new engine choice consistently to config, API, CLI, registry, provenance and resume fingerprints. Do not silently route a requested new count engine through the legacy binomial GLM.

## Work package A: numerical and inference contracts

- Fix zero-p-value handling in signed Stouffer: finite p=0 is underflow, not an observation to discard. Preserve negative/greater-than-one rejection and NaN handling; prefer signed-statistic/log-tail propagation where available.
- Separate normalized probability evaluation from relative log-likelihood evaluation. Omit count-only log-binomial coefficients in optimization and LR differences, while retaining normalized public `logpmf` behavior. Compare against SciPy/direct optimization at rho=0, near-zero rho, endpoint means, unequal depths and missing coverage.
- Record dispersion optimization diagnostics separately from mean-fit convergence. Retain robust fallback and cancellation-safe branches.
- Make failed permutation scans fatal for inference rather than conditionally dropping them. Retain successful empty scans. Preserve the observed/mirror assignment handling required by the selected randomization procedure; do not apply a second BH adjustment to scan-adjusted max-T p-values. Label count-ratio FDR estimates as estimates, retain their Monte Carlo uncertainty, and do not present zero estimated false counts as proof of zero error.
- Mark the legacy adaptive sliding-window/segmentation p/q outputs as exploratory and document their selection and covariance limitations. The new regional backend supplies valid-family inference; do not suggest the chain-merge correction has fixed every old backend.

## Work package B: fast count-level inference

Add an explicit experimental `bb_score` engine, with its mathematical/numerical kernel in `_count_score.py` and prior training in `_dispersion_prior.py`.

The working observation model is independent biological-sample counts conditional on observed coverage:

`E(M_s)=N_s mu_s`, `Var(M_s)=N_s mu_s(1-mu_s)[1+(N_s-1)rho]`, and `logit(mu_s)=X_s beta`.

Mean fitting uses quasi-score IRLS weights `N_s/[1+(N_s-1)rho]`, which saturate with read depth. N=0 excludes a sample's observation. Biological dispersion is rho, not a read-depth-dependent Pearson scale.

Train a bounded observed-count dispersion prior with a mean-dependent trend and retain a binomial component. Reuse and extend count-likelihood prior machinery rather than fit noisy floored site scales with an asserted inverse-gamma conjugacy. Compare shared-rho and group-specific-rho two-group variants in development; group-specific fits need enough observations and an explicit fallback. Nuisance fitting must use observed counts only and must not access simulation truth.

Compute the group/contrast efficient quasi-score after projecting out intercept and specified nuisance design terms. The initial analytical tail is explicitly an approximation, not an exact F law. Do not fabricate denominator degrees of freedom from a shrinkage weight. Export score, variance, fitted rho, valid samples, convergence, prior revision and tail method. Calibration comparisons must include count-model null bootstrap diagnostics with original coverage/design preserved and existing BB likelihood oracles. A candidate that cannot produce useful calibrated discoveries is not promoted.

Support formulas, donor effects, batches, rank checks and explicit contrasts in this engine. Reject aliased designs and contrasts without usable group replication. Estimate methylation differences on the response scale with intervals whose conditioning/approximation is documented. Preserve p-value/effect swap invariance and count validation.

## Work package C: region-first inference

Add `method="multiscale"`, implemented in focused geometry, sample-summary and inference modules. This caller reads replicate counts directly and does not require significant CpG seeds.

Define a deterministic family from coordinates and prespecified coverage eligibility before inspecting group effects. Start with CpG-count scales `(3,5,10,20,40)`, half-window offsets, a 500 bp maximum internal gap and a 2,000 bp maximum width. Count and export the actual family; densely anchoring each scale is not assumed to be inexpensive. Record all tested intervals, including nonsignificant intervals. Minimum CpG/width parameters must be explicit and fingerprinted.

Retain one contribution per biological sample or paired donor in each region. Use nuisance-adjusted contributions from the count model and an across-sample regional variance estimate, with robust moderation over comparable regions. Respect site baselines and missingness so differing coverage cannot create an apparent group contrast by averaging different CpGs. Do not model a sum of heterogeneous CpGs as an exact single ordinary BB observation.

The tested null is absence of the specified directional average group effect for the declared interval. It is distinct from every CpG being null, minimum fraction of affected CpGs, or differential variability. Export the tested interval identity, mean effect, uncertainty, score, p/q, sample support and inference method.

Initially correct the complete fixed interval family genome-wide using BY; assess BH only with a supported dependence argument/calibration, not as a power switch. Overlapping interval output retains interval-level hypotheses and q-values. Any deduplicated/merged display is a presentation view and must not inherit a claim of region-level FDR. Add independently defined parent-region omnibus testing or complete-discovery resampling before producing inferential merged-region outputs.

The spatial pooling model must account for within-sample covariance, rather than sum marginal CpG z-scores with an independence denominator. Pairing is handled in the sample model and within-pair randomization where valid. For observational covariates/unequal variances, use a documented model-based resampling scheme instead of assuming label exchangeability.

Existing smoothing paths retain their explicit limitations. New inferential pooling uses raw counts plus regional variance propagation, addressing the demonstrated averaging/floor problem without globally permitting arbitrary underdispersion. Add tests demonstrating correct covariance behavior in independent, perfectly correlated and intermediate cases.

## Work package D: memory and computational scaling

Separate bounded prior training from blockwise count fitting. Add blocks of at most 32,768 sites and halos where regional geometry requires overlap; emitted rows occur once. Keep the store as the result handle and materialize only on request. Avoid growing sample-by-whole-chromosome matrices in the new engine.

Bound threads and workers; avoid nested oversubscription. Use one thread for primary timing comparisons, with a separately reported parallel scaling experiment. Use exact genome-wide correction; introduce disk-backed sorting if the global p/q vectors dominate measured memory. Cache keys include engine/prior/tail revisions, design, eligibility, geometry and smoothing parameters.

Optimize repeated BB fitting only with likelihood/gradient oracles and boundary tests. A language rewrite is outside the selected approach. Report warm and cold timings separately for compiled kernels.

## Evaluation design

All runners write fresh result directories, source inventories, input hashes, seeds, exact settings, timings, sampled process-tree RSS and exit status. Cache validity and result completeness are verified from actual files and process handles. Do not restart a live job merely because a polling call times out. Determine whether existing full real-data jobs are actually live before sharing timing resources.

1. **Regression and numerical validation:** run relevant existing tests, add failing tests for the observed bugs, then implementation and green tests. Use independent count likelihood/score oracles, not tests that repeat the implementation. Broaden to the supported core suite after integration and record unavailable optional tests honestly.
2. **Development matrix:** use independent new seeds with null and mixed-signal data; 3, 5 and 10 samples/group; depth 5,20,80; equal/unequal libraries; mean-dependent and unequal group dispersion; boundaries near 0/1; missingness and outliers. Include paired donors and batches. Develop on these runs, then freeze code/settings before final holdout.
3. **Untouched holdout:** generate repeated null/signal region datasets with realistic spatial dependence and varied lengths/CpG density, plus misspecified backgrounds. Existing `correlated=False` results remain historical development comparisons. The actual calibration field is `spatial_correlation`: it contains measured correlations and a 2,000 bp fitted length but `validated=false`, and the simulator refuses correlated scenarios in that state. Validate a separate correlated generator against simulated count correlations and the retained empirical summaries; do not change the flag without evidence. Include short-range dependence and persistent sample/donor effects rather than assume a single decay parameter explains every measured correlation.
4. **Recorded full simulations:** run changed methods on the existing 22-autosome signal/null data with the original common coverage mask, alongside the corrected baseline. Use native-eligibility experiments separately, with both eligible and full-truth recall denominators.
5. **Real cohorts:** GSE64177 with all six donor pairs and donor-adjusted contrasts; GSE263850 with recorded group/replicate structure. Verify sample independence, pairing identifiers, assembly and coordinate conventions before use. Compare the original profiles and a harmonized new-method profile. Real paper-list overlap is descriptive concordance, not biological precision/recall or a tuning target.
6. **Comparative report:** include CpG ranking, discoveries, average FDP over simulation runs, Monte Carlo uncertainty, complete-null rejection probability, power, effect interval coverage, direction errors, 50%/80% reciprocal region matching, wholly null regions, base-level false coverage, candidate-stage losses, boundaries, runtime and memory. Compare nominal-q calibration separately from power curves at common realized error. Retain failures and incomplete outputs explicitly.

DSS and dmrseq remain existing cross-tool reference workflows. Add DSS-general/edgeR method comparators where supported by the installed R environment; distinguish unavailable comparators from zero calls. Avoid rerunning unrelated expensive callers without a concrete comparison need.

## Acceptance and completion

The goal is complete only when code, regression/numerical validation, fresh simulation evidence, both real-data runs, and the final comparative analysis are present and inspected. A smoke dataset or a background job launch is not completion.

Aim for new primary workflows within 2x corrected-baseline runtime and 1.25x memory on the common full simulation, measured under comparable execution conditions; these are engineering targets, not a justification to abandon the statistical improvements. Report tradeoffs if targets are missed.

Promote a new recommended mode only if held-out evidence supports error control and useful power. Demonstrate greater power at comparable realized error with uncertainty, rather than increased nominal-q call counts. An inconclusive result requires more evidence; a negative result requires a clear comparison and retained experimental labeling. Integration into external source trees or public release is a separate final action after the implementation and benchmarks are reviewable.

After written-spec review, write separate implementation plans for A/B and C/D/evaluation with concrete signatures and TDD tasks. The first plan must connect to the remaining work so numerical patches cannot be mistaken for completion of the broader objective.
