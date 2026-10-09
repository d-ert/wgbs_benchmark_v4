# Epykit numerical contracts and count inference implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans for native execution, or superpowers:subagent-driven-development if the user selects delegation. Follow the approved specification and task-by-task TDD. This plan is part one; part two owns regional inference and the complete simulation/real-data evaluation.

**Goal:** Implement the reviewed numerical corrections and a fast experimental depth-aware count engine with covariate support.

**Architecture:** Retain Parquet stores, existing engine schemas and corrected chain-merge. Add focused prior, quasi-score and block-loading modules; integrate them through config/API/CLI without routing formulas through legacy inference. Reference tails remain explicitly approximate pending the evaluation plan.

**Tech stack:** Python >=3.10, NumPy, SciPy, Numba, Polars, PyArrow, Patsy; existing pytest/R tools.

**Spec:** [Approved design](../specs/2026-10-09-epykit-statistical-improvements-design.md).

## Global constraints

- Writes remain under `/scratch/wgbs_benchmark_v4` or `/tmp`; external and frozen source trees remain inputs.
- Preserve Python >=3.10 compatibility, canonical result columns, deprecated surfaces and current default selections.
- Keep new count/region methods experimental until independent evaluation supports calibration and power.
- New fitting blocks contain at most 32,768 sites; missing N=0 is excluded, never interpreted as zero methylation.
- Main timing comparisons use one thread; cache fingerprints include design, eligibility, algorithm/prior/tail revisions.
- No arbitrary F-denominator degrees of freedom, no tuning on final holdout truth, and no completion after numerical patches alone.

## Review focus

- Formula dispatch must retain `bb_score`, including `materialize=False`; test metadata and actual kernel selection.
- Missing samples can destroy local design rank; mask the site and retain its hypothesis-family membership.
- Endpoint counts and rho approaching zero need stable likelihood/score behavior and honest interval availability.
- Failed/empty/self/mirror permutation scans have different inferential meanings; test each explicitly.
- Training, batching and row order must not leak truth or silently change group-swap results/cache validity.

## Paths and verification environment

`P` below means `/scratch/wgbs_benchmark_v4/development/statistical_improvements_20261009/epykit`, an isolated copy/check-out of `/scratch/epykit-calibration-fresh_20261007`. All package paths are relative to P. Package commands use its source through `PYTHONPATH`, with bytecode/Numba/matplotlib caches under `/tmp`. Save every red/green result under the development directory. A writable local git repository may carry commits; when unavailable, retain source inventories and patch files rather than attempt read-only git writes.

Baseline is recorded in `development/statistical_improvements_20261009/baseline-summary.json`: 791 passed, 28 failed, 6 skipped. Twenty-five failures require bioframe, one requires optional diptest, and two tests assert superseded reference/permutation behavior. Establish a writable isolated test environment using available dependencies or supported tooling; report unavailable extras rather than weaken unrelated assertions.

### Task 1: Isolated source and zero-tail/permutation contracts

**Files:** modify `src/epykit/dmr.py`, `src/epykit/dmr_segment.py`; create `tests/test_statistical_contracts.py`; update affected legacy permutation expectation tests and statistical limitations docs.

**Interfaces:** existing `_stouffer_combine_signed(pvals, meth_diffs, weights=None)` and `_aggregate_region_perm_results(...)` retain signatures.

- [ ] Verify P's source hashes against the recorded base before editing; create a writable local checkout/source copy through using-git-worktrees.
- [ ] Write `test_zero_tail_is_retained`: concordant `[0,.5]` gives finite p < 1e-100 and `[0,0]` does not become NaN. Write invalid-input cases for NaN, negative p, p>1 and invalid weights.
- [ ] Write `test_region_count_ratio_refuses_failed_scan`: `results=[(False,None)]` raises, while `(False, empty_array)` is a successful zero contribution. Keep observed/mirror assignments according to the selected method's explicit contract.
- [ ] Run these tests and retain expected failures before code changes.
- [ ] Implement underflow handling before clipping/filtering; refuse failed scans in both aggregation branches and expose Monte Carlo sample counts/limitations. Mark legacy selected-span q-values exploratory through metadata/documentation without removing API columns.
- [ ] Run new tests plus existing scan-denominator, region-count-ratio, segmentation and raw-p-window suites. Update the stale max-T test to require refusal, preserving independent denominator tests.
- [ ] Commit or export the verified patch and source inventory.

### Task 2: Relative BB likelihood and optimizer diagnostics

**Files:** modify `_beta_binomial.py`; create `tests/test_bb_relative_likelihood.py`; extend likelihood/inference oracle tests.

**Interfaces:** public `logpmf`, `fit_mean`, `fit_fixed_dispersion`, `test_lognormal` retain normalized likelihood contracts. Internal `_relative_group_loglik(m,n,mu,rho,start,stop) -> float` omits count-only coefficients. Public results add diagnostic fields without replacing existing ones.

- [ ] Write an independent oracle test: normalized LL minus `sum(gammaln(n+1)-gammaln(m+1)-gammaln(n-m+1))` equals relative LL for rho `0,1e-8,.01,.4`, depths `5,20,300`, missing observations and endpoint means.
- [ ] Write `test_optimizer_diagnostics_cover_rho_search`: result includes evaluations, search convergence, boundary indicator and objective gap/fallback status independently of mean-fit flags.
- [ ] Run red tests.
- [ ] Add the relative kernel and use it consistently inside objectives and LR calculations; restore constants only where normalized LL is exposed. Retain cancellation-safe small-rho branches and a dense/bracketed fallback when diagnostics fail.
- [ ] Run likelihood, profile-interval and direct-objective oracles; compare warm/cold count-kernel timing on identical generated counts. Promote the optimization only with equivalent public values.
- [ ] Commit/export the verified patch.

### Task 3: Count-score kernel with known biological dispersion

**Files:** create `src/epykit/_count_score.py`, `tests/test_count_score_oracles.py`.

**Interfaces:**

`fit_quasi_means(m,n,design,rho,*,max_iter=40,tol=1e-8) -> dict[str,np.ndarray]` returns coefficients, fitted probabilities, information, valid counts and convergence.

`test_score(m,n,design,contrast,rho,*,intervals=True) -> dict[str,np.ndarray]` returns efficient score, variance, statistic, pvalue, response-scale effect/interval and diagnostics. m/n have shape `(sites,samples)`, design `(samples,p)`, contrast `(p,)`; rho broadcasts by site/sample.

- [ ] Write known-rho oracle tests using independently constructed `X.T W X`, score vectors and nuisance projection. At rho=0 verify the binomial efficient score, not a likelihood-ratio statistic.
- [ ] Write depth-saturation and group-swap tests; at rho=.1 weights for depths 100 and 5 have ratio `(100/10.9)/(5/1.4)`.
- [ ] Write invalid count/design tests, locally missing/rank-deficient designs, identical endpoint groups and extreme-depth stability cases. Uninformative sites get no rejection; unavailable boundary intervals are explicitly flagged rather than reported artificially narrow.
- [ ] Run red tests.
- [ ] Implement stable quasi-score IRLS using `U=X.T[(m-n*mu)/(1+(n-1)*rho)]` and information `X.T diag[n*mu*(1-mu)/(1+(n-1)*rho)] X`. Project nuisance terms for the one-dimensional contrast. Use log-normal survival for the signed score approximation; record that it is approximate.
- [ ] Implement full-model response predictions under group counterfactuals with nuisance values fixed and delta-method intervals for regular interior fits. Reject unestimable contrasts and expose interval limitations near separation.
- [ ] Run all independent oracles and a known-rho small-simulation calibration diagnostic, recording limits rather than assuming asymptotic accuracy.
- [ ] Commit/export.

### Task 4: Mean-dependent count prior and estimated-dispersion inference

**Files:** create `_dispersion_prior.py`, `tests/test_count_dispersion_prior.py`; extend `_count_score.py`.

**Interfaces:**

`fit_dispersion_prior(m,n,n_case,*,design=None,max_training=4096,seed=20261009) -> dict`.

`posterior_dispersion(m,n,mu,prior,*,group_labels=None) -> dict[str,np.ndarray]`.

`test_counts(m,n,design,contrast,*,prior,n_case,dispersion_mode="shared",intervals=True) -> dict`.

- [ ] Write tests retaining a rho=0 component, finite posterior weights, no use of truth fields, group-swap invariance, missing-coverage exclusion and reproducible bounded training.
- [ ] Write a synthetic count-likelihood prior test with two mean strata and different dispersions; verify separate trends outperform a single pooled prior at recovering the known ordering, without asserting exact noisy estimates.
- [ ] Run red tests.
- [ ] Train on observed replicate counts only. Start with folded mean-bin edges `(0,.05,.15,.30,.50)`, count-mixture mean nodes `(1-cos(linspace(0,pi,41)))/2`, rho nodes `[0]+geomspace(1e-4,.8,25)`, deterministic bounded sampling and pooled-prior fallback when a bin has fewer than 64 informative replicate groups. Reuse normalized count-likelihood machinery; do not project away the binomial component.
- [ ] For general designs, use fitted per-sample means when evaluating count likelihood and record the conditional training approximation; do not silently assume all donor means equal. Estimate rho posterior moments and refit means with depth-saturating weights until stable, with a small documented iteration budget and fallback diagnostics.
- [ ] Implement shared and group-specific two-group variants; retain a shared-prior fallback when local group support is inadequate. Record both biological dispersion and posterior uncertainty; neither becomes invented reference df.
- [ ] Run null-count bootstrap diagnostics at fixed original coverage/design. Compare analytical score tails against the existing BB likelihood and bootstrap before choosing any recommended reference. Keep the engine experimental if small-sample tails fail.
- [ ] Commit/export.

### Task 5: Blockwise API/CLI integration and provenance

**Files:** create `_count_blocks.py`, `_count_score_store.py`; modify `_dmc_engines.py`, `_dmc_config.py`, `_dmc_stages.py`, `dmc.py`, `tl.py`, `cli/_dmc.py`; create `tests/test_count_score_engine.py`.

**Interfaces:** `iter_count_blocks(store,chrom,samples,canonical_positions,*,block_size=32768)` yields coordinates and m/n blocks. `run_count_score_store(md,cfg,*,design,contrast) -> DMCStore` owns bounded training then fitting. The public choice is `test="bb_score"`; model options default to shared dispersion and analytical diagnostic tails and are explicitly recorded.

Add keyword-only API/config fields `score_dispersion="shared"` (choices shared/group), `score_ci=True`, `score_reference="normal"` (the supported initial analytical approximation), `score_block_size=32768`, `score_prior_max_groups=4096`, and `score_prior_seed=20261009`. Add CLI flags `--score-dispersion`, `--no-score-ci`, `--score-reference`, and `--score-block-size`. Fingerprint all fields and resolved design/contrast; legacy `reference="adaptive"` is not silently described as the new engine's reference. General formulas use existing `_glm.build_design(..., return_design_info=True)` and `_glm.resolve_contrast`, preserving the group direction explicitly.

- [ ] Write end-to-end tests for two-group, `~ donor + group`, `~ batch + group`, explicit contrast and `materialize=False` execution. Assert the kernel is actually bb_score and metadata/key is `dmc_bb_score`, not `glm_contrast`.
- [ ] Write chunk-size invariance and peak-allocation tests with >32,768 sites, sparse union coverage, varying sample count and duplicate-coordinate validation. Assert no whole-chromosome sample matrix enters the new fitting path.
- [ ] Write cache invalidation tests for design, prior/tail revision, eligibility and inference options; cover API and CLI validation.
- [ ] Run red tests.
- [ ] Register/integrate the new engine before generic formula dispatch, using existing Patsy design construction and nuisance/contrast rank checks. Add canonical results plus `score_stat`, `score_variance`, `rho`, valid samples, convergence and prior/tail revision columns. Keep adjusted effects separate from raw descriptive means.
- [ ] Stream Parquet output and use exact global correction. Retain the full hypothesis count, assigning untestable sites p=1 for correction while preserving unavailable-statistic status in output.
- [ ] Run the new integration tests, supported DMC/CLI/config/cache tests and the full supported baseline suite. Address stale statistical tests by updating their independent mathematical reference, not by weakening assertions.
- [ ] Commit/export and write a source snapshot manifest used by part two.

## Handoff to part two

Part one completion supplies the known-rho and estimated-rho kernels, count-block iterator, frozen source revision and stable schemas. It does not demonstrate improved genome-wide power or finish the user objective. Continue with [regional inference and evaluation](2026-10-09-epykit-regions-and-evaluation.md) after its review gate.
