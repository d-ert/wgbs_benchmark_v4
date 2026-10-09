# Epykit regional inference and comparative evaluation implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans for native execution or superpowers:subagent-driven-development if the user selects delegation. Read the approved spec and part-one plan. Finishing requires inspected simulation and both real-data results, not only unit tests or launched jobs.

**Goal:** Recover weak coherent methylation changes with sample-aware multiscale testing and determine whether the implemented improvements outperform corrected epykit at defensible error rates.

**Architecture:** Define intervals from coordinates, aggregate contributions within biological samples, and fit/moderate regional contrasts with explicit global families. Preserve tested interval identities; use fixed parent hypotheses for inferential region output. A separate reproducible evaluation pipeline freezes code/tuning before untouched holdout runs.

**Tech stack:** Python >=3.10, NumPy/SciPy, Polars/PyArrow, existing immutable-truth scorer, installed DSS/edgeR/dmrseq R packages.

**Spec:** [Approved design](../specs/2026-10-09-epykit-statistical-improvements-design.md).

## Global constraints

- Apply part one's Python, filesystem, experimental-labeling, schema and provenance constraints.
- Geometry defaults: scales `(3,5,10,20,40)`, half-window offsets, maximum internal gap 500 bp and maximum width 2,000 bp.
- New blocks have at most 32,768 sites plus explicit geometry halos; emit each site/interval once.
- Retain the complete family; BY initially addresses interval dependence but never substitutes for calibrated raw tests.
- Primary timing uses one thread; engineering targets are <=2x baseline runtime and <=1.25x memory, with missed targets reported.
- Preserve original data/scoring. Paper-list concordance is descriptive and cannot train the new method.

## Review focus

- Overlapping windows do not create extra biological replicates; perfectly correlated CpGs must not gain sqrt(K) independent evidence.
- Variable coverage/site baselines must not manufacture a group contrast; exercise sparse and confounded coverage fixtures.
- Boundary/halo processing must emit each prespecified hypothesis exactly once, including small chromosomes and long gaps.
- Merged display boundaries cannot inherit window q-values; tested parent identity must remain visible.
- Simulation seeds, tuning and source hashes must be frozen before untouched testing; failed/pending runs cannot become zero discoveries.

## Paths and dependencies

P is the isolated package path defined in part one. Evaluation code lives in `development/statistical_improvements_20261009/evaluation/`; fresh outputs live under `results/statistical_improvements_20261009/`. These locations are writable. Reuse part one's `test_counts`, count-block iterator and source manifest. Existing comparison functions are `wgbs_v3.score.dml_metrics`, `dmr_metrics`, `score` and `wgbs_real.comparison.compare`; use them without redefining truth or matching rules.

### Task 1: Complete multiscale geometry and test identities

**Files:** create P's `_region_geometry.py`, `tests/test_multiscale_geometry.py`.

**Interfaces:** `multiscale_intervals(positions,*,scales=(3,5,10,20,40),max_gap=500,max_width=2000) -> structured ndarray` with half-open index and genomic boundaries, scale and stable test ID. `parent_regions(positions,*,max_cpgs=40,max_gap=500,max_width=2000)` returns disjoint coordinate-only parents with their child intervals.

- [ ] Write a brute-force geometry oracle covering empty inputs, 2/3-CpG chromosomes, duplicate/unsorted positions, exact 500 bp gaps, 2,000 bp widths and overlapping scales. Assert half-window stepping uses `max(1,k//2)` and tail inclusion is deterministic.
- [ ] Write parent assignment/halo equivalence tests: every tested child belongs to one fixed parent, no interval crosses a forbidden gap, and different processing blocks produce identical IDs.
- [ ] Run red tests, implement validated geometry and coordinate-only parent partitioning, then run oracle tests.
- [ ] Export the actual interval and parent-family counts; commit/export.

### Task 2: Biological-sample summaries and moderated region tests

**Files:** create `_region_sample_scores.py`, `_region_inference.py`, `tests/test_region_sample_inference.py`.

**Interfaces:** `region_contributions(m,n,positions,intervals,*,null_means,rho,design) -> dict` returns sample/region contributions and usable support. `fit_region_variance_prior(residual_variance,residual_df,*,robust=True) -> dict`. `test_region_contrasts(contributions,design,contrast,*,variance_prior) -> dict` returns effects, intervals, p/statistic, residual df and diagnostics.

- [ ] Write independent sample-level regression/contrast oracles for unpaired, donor-paired and batch designs; test rank deficiency and fewer covered samples than model coefficients.
- [ ] Write covariance fixtures at independent, intermediate and perfectly correlated CpGs. Duplicating identical CpG observations must not spuriously multiply independent-sample information. A homogeneous independent fixture must exhibit the expected averaging gain without flooring its variance at the unsmoothed count variance.
- [ ] Write null fixtures with site-dependent methylation and group-dependent missing coverage; assert consistent support or baseline adjustment prevents a systematic contrast. Export exclusion reasons when support is insufficient.
- [ ] Run red tests.
- [ ] Aggregate nuisance-adjusted count-model contributions with prefix sums, retaining each sample/donor. Use across-sample residual variance and robust empirical-Bayes moderation of unfloored residual variances; estimate prior strength from the sampling-variance model, with numerical/stratum diagnostics. Do not treat marginal CpG p-values as independent.
- [ ] Compute response-scale regional effects from explicit consistent-site summaries and keep them separate from standardized score contrasts. Document approximate tails and effect-interval conditioning.
- [ ] Run regression/covariance oracles plus a bounded correlated-null diagnostic; commit/export.

### Task 3: Interval and parent-level inference, API and resampling

**Files:** create `_region_multiscale.py`; modify `tl.py`, `cli/_dmr.py`, metadata/exports; create `tests/test_multiscale_engine.py` and `tests/test_parent_region_inference.py`.

**Interfaces:** `call_dmr_multiscale(md,*,formula=None,contrast=None,scales=(3,5,10,20,40),max_gap=500,max_width=2000,fdr_method="fdr_by",output="parents",...) -> pl.DataFrame`. Detailed interval tests are stored separately with IDs and complete-family q-values. Parent output uses independently declared parent nulls.

- [ ] Write an API fixture with weak concordant effects and no significant individual CpG seeds; assert intervals are tested and the parent hypothesis can be discovered when evidence is adequate.
- [ ] Write an independent parent omnibus oracle: for each fixed parent, `p_parent=min(1,number_of_children*min(child_p))`, treating unavailable children as p=1. Apply genome-wide parent BY and retain parent boundaries/IDs. Merged display views must not acquire inferential q-values.
- [ ] Write metadata/cache and two-study sample-design tests; expose unavailable statistics and experimental inference labels.
- [ ] Run red tests and implement streaming region analysis plus interval/parent correction and raw-count API/CLI dispatch.
- [ ] Add complete-pipeline paired randomization and count-model null-bootstrap diagnostics with original coverage/design preserved. Enumerate small paired assignment spaces where valid; reject failures rather than drop them. Export attainable resolution and separate FWER p-values from estimated FDR.
- [ ] Run integration/geometry/covariance/permutation tests; commit/export. Compare parent vs interval output under the unchanged one-to-one region scorer, documenting their different boundary targets.

### Task 4: Bounded memory, exact global correction and timing instrumentation

**Files:** create P's `_external_correction.py`; extend count/region store modules and tests; create evaluation `measure.py`.

**Interfaces:** `adjust_store(store,*,method="fdr_by",max_rows_in_memory=1000000) -> store` performs exact correction through bounded sorted runs and reverse cumulative minima. `run_measured(command,out,*,thread_budget=1) -> dict` records live process handle, status, timing and process-tree RSS.

- [ ] Write small-family BH/BY reference tests with ties, NaNs, p=0/1 and multiple sorted runs; match SciPy/statsmodels over the same full hypothesis count.
- [ ] Write block/halo output-equivalence tests and allocation limits growing sites/sample count; verify materialize=False avoids full result collection.
- [ ] Run red tests, implement bounded correction and timers, and verify values against in-memory correction.
- [ ] Run warm/cold kernel profiles and one-/bounded-multiworker probes without oversubscribed libraries. Save all timings; commit/export.

### Task 5: Correlated simulation generator and independent validation

**Files:** create evaluation `simulate.py`, `validate_generator.py`, `tests/test_generator_contracts.py`, `development_matrix.json`, `holdout_protocol.json`.

**Interfaces:** `generate_counts(config,out) -> manifest` with immutable latent truth, separately seeded biological/count/coverage/effect streams, count Parquet and compatible samples/truth files. `validate_generator(manifest) -> dict` checks model/count/spatial/design properties.

- [ ] Write tests for integer 0<=M<=N, missingness, donor/batch labels, ground-truth intervals and coordinate semantics; verify count sampling noise and biological variance independently.
- [ ] Write correlation tests checking independent, short-range Gaussian-copula and persistent sample/donor components. Use numerical binomial/BB marginals and large seeded validation draws to verify count correlations within Monte Carlo intervals.
- [ ] Run red tests and implement a separate generator. Preserve the old empirical manifest with `validated=false`; do not bypass its simulator guard. Compare generated distance-bin correlations against retained empirical summaries and label departures/misspecification.
- [ ] Populate the development matrix with 3/5/10 samples/group, depths 5/20/80, balanced/imbalanced libraries, shared/unequal dispersion, endpoints, missingness, outliers, paired/batch effects and region length/density scenarios. Use fresh seeds derived from a recorded development namespace.
- [ ] Run development comparisons for corrected baseline, BB-F, bb_score and multiscale. Iterate defects/tuning under TDD and retain all runs. Use scope-preserving plan adjustments if calibration exposes a better justified test; never select a tail solely because it makes more discoveries.
- [ ] Freeze source, parameters and protocol hashes. Generate final holdout seeds from a separate recorded namespace only after freeze; changing tuned code invalidates that final evaluation and requires a new untouched namespace.
- [ ] Run at least 50 repeated 100k-site complete-null and 50 mixed-signal holdouts across prespecified spatial/design strata. Record Monte Carlo uncertainty and selection of any additional runs. These do not certify probabilities near 1e-9; mathematical tail/oracle evidence remains required.

### Task 6: Full simulation and real-data comparisons

**Files:** create evaluation `run.py`, `run_epykit.py`, `run_reference.R`, `compare.py`; create immutable source/input/run inventories.

**Interfaces:** `python evaluation/run.py --stage {development,holdout,full-simulation,real,compare} --run-dir PATH` resumes only verified completed stages and polls live handles before restart. Native eligibility is a separately named experiment.

- [ ] Write runner tests for fresh destinations, failed/pending statuses, source/input hash changes, coordinate conversion and unchanged matching definitions. Use tiny fixtures only for integration, not performance claims.
- [ ] Implement runner commands selecting each implemented engine/caller explicitly. Preserve canonical outputs and failed-stage diagnostics; prohibit accidental legacy-GLM dispatch for bb_score.
- [ ] Run existing 22-autosome signal/null data with the unchanged common mask. Rerun corrected epykit alongside candidates for comparable timing; retain original/BB-F and existing DSS/dmrseq results as labeled reference workflows.
- [ ] Run GSE64177 on all 12 samples with `~ donor + group`, verifying six pairs and hg19 coordinates; compare paired baseline/new count inference and corresponding regions. Preserve the historical unpaired profile as a descriptive reference.
- [ ] Run GSE263850 on all six samples, three treatment/three control, hg38, explicit unpaired design. Keep common-parental-line/clone limitations with the comparison.
- [ ] Run installed DSS-general and edgeR reference profiles on development/holdout and both real studies as practical; record exact design and hypothesis differences. Verify whether existing expensive real callers are live before timing overlapping work.
- [ ] Inspect completed outputs, recompute candidate counts and confusion/matching metrics from immutable truth/reference, verify source/input hashes and inspect live/terminal process state. Every full new-method simulation and both real studies must finish before completion.

### Task 7: Comparative evidence, review and delivery

**Files:** evaluation `report.md`, `metrics.csv`, `validation.json`, final package patch/source inventory and `state.json`.

- [ ] Independently recompute CpG/region metrics and denominator arithmetic. Report nominal-q calibration and ranking/power at matched realized error separately, including Monte Carlo uncertainty.
- [ ] Include direction/interval coverage, 50%/80% overlap, wholly null regions, false base coverage, boundary errors, candidate-stage losses, paired-design impact, time/memory and every failed/inconclusive stratum. Paper overlap remains descriptive.
- [ ] Run the integrated supported core suite and appropriate slow numerical/statistical tests; compare against diagnosed baseline failures. Run lint/type checks supported by the actual environment and record unavailable tooling/platforms.
- [ ] Use requesting-code-review and verification-before-completion for an independent final branch review. Fix material issues, rerun affected tests/evaluations and invalidate holdout results if post-freeze inferential behavior changed.
- [ ] Audit every approved-spec requirement against inspected files and completed outputs. Deliver a runnable source copy/patch and the comparative report. Keep experimental labeling where evidence is insufficient and state negative results clearly.
- [ ] Only after the complete objective is demonstrated, mark the goal complete. Integration into external trees/public release remains a separate reviewable action.
