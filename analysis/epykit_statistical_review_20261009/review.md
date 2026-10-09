# Epykit statistical and power review — 9 October 2026

The supplied report is substantially sound. Preserve the fast store and corrected chain-merge baseline, keep BB-F experimental, and develop regional inference alongside a better CpG engine. A spatial Bayesian rewrite is a research option, but it should follow a cheaper region-first prototype. The largest additions to the report are a smoothing/variance-floor interaction, inefficient replicate weighting, unresolved selection correction in other region backends, and a reproducible zero-p-value bug.

This is a source and recorded-results review, with bounded diagnostic calculations. No callers, frozen snapshots, benchmark inputs, or existing reports were modified. Full WGBS workflows were not rerun. The supporting [script](checks.py) and [results](checks.json) preserve the new checks.

## What I verified

- The current `/scratch/epykit3/src/epykit/dmc.py` and `dmr.py` are byte-identical to the original benchmark snapshot. Corrected and BB results must not be presented as behavior already installed in that checkout.
- Sixteen signal region-comparison rows have consistent precision, recall and F1 arithmetic. This checks the arithmetic from recorded matching counts, not a fresh independent matching of every caller.
- Recounting the 22-chromosome geometry family gives **776,121,901 intervals**, harmonic factor **21.04703582**, **10,785 candidates**, and **980 retained calls**. A separate BY calculation reproduces the cached q-values.
- The GSE263850 native candidate cache contains **51 candidates before correction and 51 after correction**. Its configuration uses seed p < **1e-5**, minimum **4 CpGs**, minimum **51 bp**, 100 bp merging, 50% seed fraction, and count smoothing over 500 bp.
- The report's F versus chi-square example is correct: statistic 30 gives **0.000589388** under F(1,8), versus **4.32046e-8** under chi-square(1).
- Source inspection confirms the depth-dependent legacy scale, full-model MAP dispersion held fixed in the BB null mean fit, independence-based Stouffer denominator, and nested BB optimization.

The published precision/recall tables should remain descriptions of recorded workflows. Strict reciprocal matching failure is not the same event as a wholly null region. In particular, the stored signal summary reports that none of corrected epykit's 980 calls are entirely devoid of positive eligible CpGs; its 89 strict nonmatches cannot all be called biological false discoveries. The full local report already explains this distinction, and the pasted summary should retain it.

## Where I would revise the reasoning

**A pre-q loss is not automatically a geometry loss.** Candidate formation combines the CpG model, raw-p threshold, effect threshold, smoothing, geometry and seed fraction. The 51-candidate result excludes final BY as the immediate cause of that profile's low call count, but does not distinguish these upstream mechanisms. The real-study 1e-5 seed gate is 5,000 times smaller than the simulation's 0.05 gate. Neither is necessarily calibrated at its nominal value. Compare gates on held-out data with the complete inference procedure, not by optimizing agreement with the already-inspected paper list.

**Holding dispersion fixed is not intrinsically invalid.** A coherent moderated test can estimate nuisance dispersion once and account for uncertainty through its construction. Independent dispersion optimization under each hypothesis is not mandatory for all valid tests; with very small groups, unregularized alternative fits can themselves be unstable. The actual problem is that this fitted conditional LR has no demonstrated small-sample reference after MAP estimation and prior learning. MethylSig proposed a similar t-squared approximation, which explains the precedent but does not validate this modified EB engine. [MethylSig original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4147891/).

**Region power need not wait for genome-wide significant CpGs.** Regional testing can aggregate weak evidence across sites while retaining the biological sample as the unit of replication. Properly calibrated regional resampling can also use a ranking score that is not itself a calibrated CpG p-value. Improve CpG inference and regional discovery in parallel as development priorities.

**The full-family correction has a defensible conditional argument.** If every fixed interval has a valid p-value, replacing unselected intervals' p-values by one preserves marginal validity. BY can then address arbitrary dependence across the complete family. The large penalty is therefore a cost of this implementation's family and selection scheme, not an unavoidable price of WGBS. It does not repair invalid fixed-interval p-values. [Benjamini–Yekutieli](https://doi.org/10.1214/aos/1013699998).

## Additional statistical problems

### 1. Count smoothing can conflict with the dispersion floor

The active `smoothing=True` path averages methylated counts and coverage over neighboring CpGs, then passes those fractional averages into the legacy engine. See [averaging code](/scratch/wgbs_benchmark_v4/results/genome_baseline/source_snapshot/epykit/dmc.py:491) and [scale construction](/scratch/wgbs_benchmark_v4/results/genome_baseline/source_snapshot/epykit/dmc.py:775). This is distinct from `use_smoothed=True`, which reconstructs rounded pseudo-counts from smoothed methylation proportions in `_smoothed_store.py`.

For K independent homogeneous CpGs, each with depth N and biological dispersion rho, averaging counts gives a pseudo-count variance scale

`phi_smoothed = [1 + (N - 1) rho] / K`.

With K=10, N=20, rho=0.10, this is **0.29**, while average coverage remains 20. Flooring this at one overstates the known-model variance by **3.45 times**. A seeded 20,000-window, 5+5-replicate diagnostic gives mean raw Pearson scale **0.2912**, with **99.96%** below one. Underdispersion here is induced by averaging; it is not evidence of impossible biology.

This calculation establishes an implementation/model interaction, not the size of its effect in the real datasets. Positive within-sample spatial covariance changes the gain from averaging; heterogeneous site means and coverage complicate it further. Overlapping smoothed windows also correlate regional scores. These effects can create conservatism in one layer and anticonservatism in another.

**Proposed remedy:** treat smoothing as estimation with an explicit propagated variance, or test raw counts using a spatial mean model. A variance formula for unsmoothed integer counts should not be reused automatically for averages. Simply allowing phi < 1 everywhere is not a justified general fix. BSmooth explicitly combines smoothing with biological-replicate variability; it is useful precedent for separating those responsibilities. [BSmooth original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3491411/).

### 2. Legacy read weighting can waste biological information

Legacy group means are pooled M/sum(N). Under independent beta-binomial replicate proportions with a common mean, inverse-variance weights are proportional to

`N / [1 + (N - 1) rho]`.

They approach 1/rho as depth grows. At rho=0.10, a 100x sample carries about 9.17 effective reads and a 5x sample about 3.57; their information ratio is about **2.57**, not the **20** used by raw-depth pooling.

For replicate depths [5,5,5,5,100], the exact variance of the pooled mean is **1.82 times** that of the known-rho inverse-variance weighted mean. This is a controlled efficiency comparison, not an estimated speed or power gain for a fitted engine. The BB optimizer already uses these weights for initialization and goes on to fit a BB likelihood.

**Proposed remedy:** include depth-saturating weights in a fast quasi-score/IRLS prototype, coupled to robust biological-dispersion estimation. Correcting only the variance of the existing read-pooled estimator may repair calibration while retaining avoidable power loss. Do not treat these weights as the exact BB MLE or as known after estimating rho.

### 3. The correction does not cover all region backends

The [sliding-window caller](/scratch/wgbs_benchmark_v4/results/default_search_by_20261007/source_snapshot/epykit/dmr.py:652) selects and merges windows using the data, then applies BH only to the resulting spans. The [segmentation caller](/scratch/wgbs_benchmark_v4/results/beta_binomial_20261007_full/source_snapshot/epykit/dmr_segment.py:173) chooses segments from observed effects, combines their CpG p-values, and applies BH per chromosome. Both retain selection problems. Per-chromosome adjustment also does not generally control a combined genome-wide discovery list.

A simple counterexample: one independent uniform null p-value per chromosome, each tested at 0.05, has whole-genome FDR `1 - 0.95^22 = 0.6765`. This is an illustrative counterexample, not measured epykit FDR.

**Proposed remedy:** give every backend an explicit tested family, null hypothesis, selection treatment and correction scope. Until validated, label its p/q output as exploratory where appropriate. Making candidate boundaries more sophisticated does not supply selection-aware inference. Window-level FDR also does not automatically survive merging significant windows into regions; csaw documents this distinction and a geometry-based region aggregation approach. [csaw original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4797262/).

### 4. Zero p-values are discarded by Stouffer combination

The [validity mask](/scratch/wgbs_benchmark_v4/results/default_search_by_20261007/source_snapshot/epykit/dmr.py:249) requires p > 0 before clipping. I executed the extracted function:

| Inputs, both positive direction | Returned regional p |
| --- | --- |
| [0, 0.5] | 0.5 |
| [smallest positive normal float, 0.5] | 8.56e-161 |
| [0, 0] | NaN |

**Proposed remedy:** distinguish invalid inputs from numerical underflow. Prefer log-tail or signed-statistic propagation; at minimum, clip valid zero values before combination rather than omit them. The report's recorded CpG tail summaries contain no zero values, so this does **not** explain the reported benchmark failures.

### 5. Permutation shortcuts need stronger contracts

The corrected snapshot's max-T branch includes successful empty scans, aborts on failed scans, retains self/mirror assignments and uses a +1 Monte Carlo correction. Its alternate `region` count-ratio branch still excludes failed and self/mirror scans and can return zero estimated q-values. See [aggregation](/scratch/wgbs_benchmark_v4/results/default_search_by_20261007/source_snapshot/epykit/dmr.py:1826). A zero estimated false-count ratio is not proof of zero error; omitting computational failures can bias the null sample.

Use exact enumeration when the valid assignment space is small, or validated Monte Carlo inference with uncertainty. For two-sided statistics invariant to reversing every label, the smallest nonzero exact tail probability is ordinarily at least **2/64=0.03125** for six pairs and **2/252=0.00794** for 5+5 samples. Pooled cross-region null distributions can yield finer numerical resolution, but require comparable null statistics; the number of pooled regions is not the number of independent biological assignments. [Permutation p-values paper](https://arxiv.org/abs/1603.05766), [dmrseq methods](https://academic.oup.com/biostatistics/article/20/3/367/4899074).

Unequal group dispersion, confounded depth and observational covariates also mean naive label permutation may not test the desired equal-mean null. Match resampling to the study design. A count-model null bootstrap preserving coverage/design is an alternative with model assumptions; a studentized permutation is not automatically finite-sample exact.

## The development strategy I would favor

### A. Build a fast region-first prototype

1. Specify a coordinate-only family before looking at treatment effects: bounded physical widths and/or CpG-count scales, multiple offsets, and gaps respected. Count the actual family on the catalogue. A few scales do not guarantee a small family if every CpG anchors every scale.
2. Preserve one contribution per biological sample within each interval. Use fixed-site-weight methylation summaries or sums of suitably centered count-model scores, with depth-aware measurement uncertainty. Use consistent CpG support or model site baselines; otherwise differing coverage patterns can make groups average different biological sites.
3. Estimate variability across samples or donor differences. Moderate residual variances across comparable intervals, with robust protection for hypervariable regions. This avoids estimating a large unrestricted CpG covariance matrix from ten samples, but remains a small-sample model requiring calibration.
4. Test and correct the declared interval or region family. If the product reports merged regions, implement region-level inference for that output, rather than attaching window q-values to a selected union. One tractable design uses fixed parent regions with prespecified internal multiscale tests and valid within-parent/global correction. Boundary refinement can initially be descriptive, with the tested parent retained explicitly.
5. Support `~ donor + group`, batches and specified contrasts at the first prototype stage. Treat missing coverage as missing, not zero methylation.

With a fixed small number of scales, prefix sums of sample contributions can make interval scoring cheap; paired summaries are particularly simple. This is an architectural expectation, not a measured runtime claim. dmrseq supplies a relevant precedent for explicitly regional inference, while robust moderated variance methods supply a statistical precedent for protecting against unusually variable features. Neither transfers its guarantees automatically. [dmrseq implementation](https://github.com/kdkorthauer/dmrseq), [robust empirical Bayes paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5373812/).

The existing tile backend can supply some I/O machinery, but summing heterogeneous CpGs and treating the resulting totals as a single ordinary BB observation is not generally an exact count model. Preserve the distinction between a working mean/variance model and a full likelihood.

### B. Compare two CpG engines before choosing a new default

**Fast candidate:** depth-aware quasi-score/IRLS with biological dispersion, a learned mean-dispersion trend where supported, robust shrinkage and justified moderation. This is a development proposal, not an endorsement of the existing `bb_eb` moment option. Use a sample-by-covariate design matrix and a common implementation for paired and unpaired studies. Ordinary small-sample sandwich SEs alone are not sufficient at five samples per group.

**Likelihood reference:** retained BB count likelihood with a consistent penalized/profile or integrated nuisance treatment, including boundary cases and heteroscedastic alternatives. Independent profiling under both hypotheses is an option, not a universal requirement. Use selected difficult sites and simulation strata to compare accuracy and cost rather than forcing nested optimization across every site by default.

Add **DSS-general and edgeR's methylation workflow** as method-development comparators. DSS-general uses transformed beta-binomial regression for general designs; edgeR models methylated and unmethylated counts with a specialized design and supports established likelihood and quasi-likelihood testing. They offer useful speed/design precedents without asserting superiority on this benchmark. [DSS-general original paper](https://academic.oup.com/bioinformatics/article/32/10/1446/1743267), [edgeR methylation workflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC5747346.1/), [edgeR v4 paper](https://academic.oup.com/nar/article/53/2/gkaf018/7973897).

### C. Increase power through valid information use

- Replace the benchmark's all-samples intersection in a separately labeled native-eligibility experiment with prespecified minimum replication and usable-information requirements. Keep the original common-mask comparison for continuity. Report eligible-site and full-truth recall separately. Relaxation should not silently change the benchmark denominator.
- Evaluate independent filtering or hypothesis weighting only after the p-values are calibrated conditionally on the proposed covariates. Depth currently predicts null p-values, so applying depth-based weighting now could amplify the existing error. Genomic dependence must also be respected in cross-fitting. [Independent filtering](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865/), [IHW](https://www.nature.com/articles/nmeth.3885).
- Allow the scientific hypothesis to be `|Delta beta| > delta_min`, with an effect interval, rather than a zero-effect test followed by an observed-effect filter presented as evidence of a minimum biological effect. This asks a harder but more useful question. DSS offers posterior threshold probabilities for two-group tests; TREAT supplies the general threshold-testing principle. [DSS manual](https://www.bioconductor.org/packages/release/bioc/manuals/DSS/man/DSS.pdf), [TREAT original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654802/).
- Distinguish regional average change, any nonzero site, coherent same-direction change and differential variability. These are different null hypotheses; a single q-value should not imply all of them.

### D. Optimize the statistical workload

The report's algorithm-first advice is correct. One refinement: the log-binomial coefficients cancel from likelihood differences and parameter optimization. An internal relative-log-likelihood kernel can omit them entirely, retaining normalized probabilities for the public API. This is stronger than caching them, provided all compared branches use the same convention.

For repeated mean/rho fitting, compare joint or alternating derivative updates, warm starts and adaptive fallback. Keep the stable small-rho expansions and boundary branches: direct differences of digamma/trigamma functions can lose precision near the binomial limit. The current count summations are partly a numerical safeguard, not just accidental inefficiency. Require objective/gradient/boundary diagnostics, not only successful mean-fit flags.

Train priors on a bounded but representative sample stratified by mean, depth and missingness, with stability checks. Reuse fitted design structure, stream count blocks, and run expensive likelihood diagnostics only on a subset or formally valid screened family. Any screening used to save testing cost must preserve the original hypothesis family or provide an independent selection argument.

### E. Keep the spatial Bayesian rewrite as a later experiment

A joint change-point/count model may improve boundaries and pool moderate effects. But it adds prior sensitivity, state assumptions, approximate inference and a different meaning of posterior error control. Compare it against the cheaper regional baseline first. A published Bayesian methylome change-point model supports feasibility, not an expectation that this will be the fastest or best default for epykit. [Hirt and colleagues](https://arxiv.org/abs/2211.07311).

## Validation gates

Freeze the candidate family, prior procedure and tuning rules before final evaluation. Use repeated untouched seeds with both whole-null and mixed-null genomes. Include realistic spatial residual dependence, abrupt and gradual boundaries, varying lengths/densities, balanced and unequal depths, group-specific dispersions, paired donors, batches, outliers and structured missingness. Preserve empirical background methylation patterns as well as dependence; an independent beta-binomial simulator alone favors its own model family.

Report average FDP across runs with Monte Carlo uncertainty, power, effect-interval coverage, direction errors, strict 50%/80% matching, wholly null regions, base-level false coverage and boundary error. Distinguish rank performance at matched realized error from nominal-q calibration. Thresholds chosen using evaluation truth are diagnostic curves, not deployable thresholds.

Measure null tails by depth, mean, dispersion, sample count and missingness, including the discovery-relevant extreme tail. A good histogram near p=0.05 does not establish accurate genome-wide p-values. For extremely small tails, use mathematical/numerical validation and carefully designed rare-event methods where needed; modest simulation runs cannot empirically certify probabilities near 1e-9.

Keep paper-list recovery descriptive and annotate when the reference list comes from the same method family being compared. Add independent validation studies or orthogonal evidence when available. Pairing and covariates should be used consistently. Do not tune on the published list and then describe its overlap as external validation.

## Review coverage

**Assessment: share with caveats; revisions needed for a development specification.** Review complete for the scoped source questions and bounded checks; genome-wide reruns and empirical validation of proposed engines remain future work. Counts below describe this review's claim groups, not exhaustive package verification.

| Category | Observed defects | Assessment |
| --- | --- | --- |
| Usefulness/completeness of recommendations | 1 / 3 | Data architecture and optimization advice are useful; the statistical roadmap needs a region-first branch and the new smoothing/weighting findings. |
| Analytical clarity | 1 / 4 | Refine attribution of pre-q losses; preserve the distinction between strict matching and biological FDP already made in the full report. |
| Visual/interaction behavior | N/A | This review concerns the pasted technical report and source, not the existing web app. |

| Category | Observed defects | Assessment |
| --- | --- | --- |
| Source authority/confidence | 0 / 4 | Original, corrected and BB snapshots distinguished; current checkout comparison verified. |
| Numerical accuracy | 0 / 4 | Region arithmetic, F tails, complete-family correction and real candidate counts checked; no fresh full cross-tool matching or CpG tally. |
| Within-chart agreement | N/A | No chart audit requested or performed. |
| Source detail completeness | 0 / 3 | Code, recorded results and primary literature linked; some publisher full-text opens were unavailable, so use is limited to retrieved material. |
| Cross-artifact consistency | 0 / 2 | Key summary tables and recorded profiles agree; source deployment status explicitly retained. |
| Data-quality/eligibility controls | 0 / 2 | Intersection mask and paired-design limitations identified; raw BAM quality and all upstream conversions not audited. |
| Conclusion support | 2 / 5 | Avoid isolating candidate geometry from upstream testing; fixed nuisance estimates need justified inference but are not intrinsically invalid. Remaining major diagnoses supported. |

The caller-level problems above are findings about epykit, distinct from defects in the report itself. All proposed statistical replacements still require comparative implementation and validation; no power gain or speedup is claimed from this review.
