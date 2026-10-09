# Fresh default-path region correction — 2026-10-07

This work starts from commit `0566519`, the package source frozen with the
original genome benchmark. It does not include yesterday's permutation fix or
position-window screen. Development checkout:
`/scratch/epykit-calibration-fresh_20261007`, branch
`codex/wgbs-default-fix-20261007`.

Implemented source and regression tests are committed as `f73e66d`. The tested
runtime snapshot has exactly the same package source hashes as this commit.

## Problem and change

The original chain-merge caller uses CpG p-values to select region boundaries,
then applies BH correction over only the surviving candidates. The testing
denominator therefore shrinks after examining the evidence; positions with no
selected regions do not contribute at all.

The revised caller keeps exactly the same CpG calculations, seed criteria,
candidate boundaries, effects and signed-Stouffer raw region scores. It counts
all possible contiguous spans meeting the existing minimum CpG count and
inclusive physical length, without crossing a consecutive-CpG gap exceeding
the existing merge limit. This family depends on positions and geometry only.
Every possible selected chain belongs to it. Unreported spans receive p=1.

For M possible spans, let H(M) be the harmonic sum. For candidate p-values
sorted increasingly, the BY adjusted value is
`min(1, min_{j >= i}(M * H(M) * p[j] / j))`.
All omitted p=1 terms contribute adjusted value 1, so this exact calculation
requires only the reported candidates. The harmonic factor uses digamma rather
than allocating M entries. Interval counting also avoids enumerating the
family. CpG inputs from every included chromosome contribute, even chromosomes
without candidates. Cache keys include the correction revision.

BY handles dependence between overlapping interval tests and data-dependent
selection with omitted scores set to one **provided the fixed-span raw
p-values are valid**. It does not repair the original CpG dispersion model or
the within-span independence assumption of signed Stouffer. This is a correction
to the search accounting, not a claim of established biological region FDR.
Reference: Benjamini and Yekutieli (2001),
https://doi.org/10.1214/aos/1013699998.

## Verification

Two regression tests were observed to fail against the original caller, then
pass after the change. They check that a chromosome with no selected regions
still changes multiplicity, and that adjusted values equal an explicit complete
BY vector with omitted spans assigned p=1. Further checks compare interval
counting with enumeration, cover ties/zeros/missing scores and verify cache
invalidation with the same input signature.

Focused final run: **43 passed** (`targeted_tests.log`). Full suite before the
last documentation-only revision: **750 passed, 29 failed, 6 skipped**. All 29
failures were reproduced using the original frozen package source, in
`baseline_failed_checks.log`; their names are in `full_suite_failures.json`.
Most involve absent bioframe; the list also includes the existing permutation
failure expectation, LR oracle, blacklist, sex-check and UMAP tests. These are
not silently treated as passes.

Independent calibrated-normal null diagnostic: seed 2026100719, 500 scans of
1,000 independent standard normal signed CpG statistics. Original caller had
12 scans with one false call; revised caller had none. This diagnostic checks
the analytic search mechanism, not realistic CpG calibration, biological
variability or spatial dependence.

## Full rerun and comparison

`run.py` recomputes epykit on the original full autosomal signal and null
simulations and all twelve prepared GSE64177 samples. It verifies the original
input and harness hash contract, freezes the revised package source, checks
that only `dmr.py`, `tl.py` and `_region_search.py` differ, and requires byte-for-
byte identical CpG TSVs and identical caller settings after completion.
No optional runner profiles or thresholds are changed.

Run: `python3 -u development/default_search_by_20261007/run.py` from the v4 root.
Analyze after completion:
`PYTHONPATH=vendor310:src python3 development/default_search_by_20261007/analyze.py`.

Candidate outputs: `results/default_search_by_20261007/`. All original six-tool
simulation measurements are reused for comparison, explicitly labeled.
Only original epykit and DSS have completed real results; the original overall
baseline remains failed at methylKit. The real comparison uses the published
GSE64177 list and is descriptive. It does not infer truth from paper overlap,
and preserves the original unpaired model despite the study's donor pairing.
The real paper is not used to choose the correction or tune parameters.

The final findings and tables belong in the candidate's `comparison/` directory.
Independent biological simulations with new seeds and spatial correlation are
still needed before a calibrated region-FDR claim or general-performance claim.
