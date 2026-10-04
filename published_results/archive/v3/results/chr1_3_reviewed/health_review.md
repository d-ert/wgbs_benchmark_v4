# Health review — 30 September 2026

The run completed successfully and is useful for diagnosing performance. Statistical calibration is not healthy enough to treat it as a final general comparison of the tools. The strongest concern is excess significant null CpGs; the clearest positive result is credible region detection by dmrseq, with substantial recovery improvement from DSS smoothing.

## Execution and integrity

All six manifest entries have exit code 0 and status ok. Total caller time was 2.632 hours (2 h 38 min), excluding scoring. Peak recorded memory was 12.49 GiB for BSmooth; all tools stayed below the 48 GiB target. DSS dropped from 61.98 GiB process-tree RSS in the original multicore run to 4.83 GiB with one core and smoothing. These are measured process-tree RSS values, not physical-memory comparisons corrected for shared pages.

Dataset metadata and saved runner-source hashes pass verification. The existing dataset validation passed. The health audit checks returned CpG keys for duplication and membership in the shared eligible universe and independently recomputes global BH adjusted p-values. All four native CpG p-value outputs agree with BH to within 5e-13. DSS returned 2,656,764 rows, 20 fewer than the common eligible universe; those sites are ranked last. No returned native p-values were nonfinite.

BSmooth logged 21 warnings without their text. Their contents cannot be recovered from this completed R process. This remains an unresolved execution qualification; future runs should print warnings immediately. No error was recorded. DMRcate's zero discoveries are explicitly logged as a completed analysis with no individually significant CpGs.

## Region comparison

Strict matching requires at least 50% reciprocal overlap in eligible CpGs and one-to-one matching. Recall divides by all 659 implanted regions. Fixed-window methylKit results have different boundary constraints from adaptive callers.

| Tool | Matched / called | Precision | Recall | Calls with no positive eligible CpG |
|---|---:|---:|---:|---:|
| epykit | 231 / 619 | 37.3% | 35.1% | 341 |
| DSS | 172 / 287 | 59.9% | 26.1% | 14 |
| methylKit | 65 / 744 | 8.7% | 9.9% | 620 |
| BSmooth | 81 / 204 | 39.7% | 12.3% | 27 |
| dmrseq | 216 / 241 | 89.6% | 32.8% | 5 |
| DMRcate | 0 / 0 | undefined | 0% | 0 |

The machine scorer represents precision with no calls as zero by convention; biologically there is no observed precision estimate for DMRcate.

dmrseq offers the strongest observed region precision–recall balance here (F1 0.480). Only five called regions contain no positive eligible CpG. Its 25 unmatched calls therefore should not all be described as entirely spurious regions: 20 overlap true signal but fail strict matching.

DSS improved from 3 matched regions to 172 after recommended WGBS smoothing. Of its 287 calls, 273 overlap at least one positive eligible CpG, but only 172 pass strict matching. BSmooth likewise has 177 of 204 calls overlapping signal, with 81 strict matches. Boundary accuracy and fragmentation are material parts of their apparent weakness. Mere overlap is a diagnostic, not a replacement success criterion: broad calls can overlap signal by chance.

BSmooth's documented size/effect filters reduced calls from 388 to 204 while retaining all 81 strict matches. epykit and dmrseq reproduced the original region counts. methylKit MN/BH reduced calls from 158,232 to 744, a substantial improvement, but most remaining calls still lack implanted signal. This comparison changes both overdispersion handling and adjustment method together, so it cannot attribute the improvement to either alone.

## CpG calibration concern

| Tool | True / false significant CpGs | Observed false-discovery proportion | False calls >2 kb from any implanted positive CpG |
|---|---:|---:|---:|
| epykit | 91 / 357 | 79.7% | 355 |
| DSS | 3,911 / 19,849 | 83.5% | 18,482 |
| methylKit | 412 / 3,536 | 89.6% | 3,522 |
| DMRcate | 0 / 0 | no discoveries | 0 |

These are realized error proportions on one dataset, not repeated-run estimates of FDR. Nevertheless, the discrepancy from the nominal 5% level is severe. Most false positives are far from any implanted positive CpG, including positives excluded by coverage filtering. Boundary leakage alone cannot explain the discrepancy. Global BH arithmetic is correct, so the next investigation should examine the underlying p-values, variance estimation, and data assumptions.

Among null CpGs more than 2 kb from truth, the proportions with p<0.001 are 0.289% for epykit, 1.027% for DSS, 0.435% for methylKit, and 0.134% for DMRcate. A calibrated continuous null would yield 0.1% on average. These observed excesses support concern about small-p-value calibration; spatial dependence and a single realization preclude treating these as independent binomial samples.

Existing small complete-null pilots already produced 1, 89, and 16 significant CpGs for epykit, smoothed DSS, and corrected methylKit, respectively. Those pilots contain roughly 16,000 eligible sites, so they are useful warnings but do not establish whole-chromosome error control. One all-null replicate cannot estimate FDR reliably.

## Ranking is better than thresholded performance

BSmooth has the strongest CpG average precision (0.484), followed by DSS (0.364). This is clear evidence that both detect signal in this simulation, despite weak strict region boundaries or poorly calibrated significance thresholds. Their AUROCs are 0.927 and 0.954, respectively. High AUROC alone can conceal poor precision because only 0.216% of eligible sites are positive.

DMRcate has AUROC 0.860 and average precision 0.0604. It ranks some positives above null sites but has insufficient evidence to cross the chosen multiple-testing threshold: minimum raw p=6.45e-8, minimum BH q=0.171. Zero regions are consistent with the initial CpG gate. This supports a threshold/power limitation under the current design; it does not prove DMRcate is intrinsically weak or justify loosening the cutoff until it wins.

## Dataset constraints

Only 45.0% of all simulated CpGs and 46.5% of positive CpGs survive coverage in every sample. This imposes a 46.5% ceiling on full-truth CpG recall for this shared-filter analysis, but is not itself a ceiling on region recall. Mean sample depth ranges from 13.65× to 24.76×. The generator also redraws sample library factors by chromosome. These are material limits of the retained dataset.

The simulation has independent per-site biological variation and dispersion calibrated across heterogeneous BLUEPRINT samples. It is not a complete model of tissue-matched WGBS replicates. These assumptions can interact with smoothing and variance estimation; the present audit does not identify one of them as the proven cause of miscalibration.

## Recommended next work

1. Keep this as a completed diagnostic comparison with its source snapshots and original outputs preserved.
2. Before interpreting nominal significance as reliable, run several full-size complete-null replicates with the frozen runner settings. Measure false discoveries and the probability of any discoveries, rather than choosing settings for higher recovery on this known seed.
3. Investigate small-p-value excess by coverage, baseline methylation and dispersion, with sample-label permutations or null simulations that respect the design. Capture BSmooth warning text in the next run.
4. For the next dataset, fix persistent sample-depth factors and validate missingness/dispersion before adding independent signal seeds and a separately labeled paper-like source. The current user-requested dataset was retained unchanged.

Reproducible audit: `analysis/check_reviewed_health.py`. Detailed diagnostic values: `health_diagnostics.json`. Full metric tables and figures: `report.md`, `comparison.tsv`, `ranking.tsv`, and `native_score_curves.png` in this directory. No new benchmark or simulation was launched for this review.
