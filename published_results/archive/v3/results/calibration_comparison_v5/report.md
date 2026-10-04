# epykit calibration repair and same-input comparison

## Readiness

Development comparison only. The fixes substantially reduce false discoveries, but calibrated CpG power and competitive region recovery must both be assessed. No five-point noninferiority claim is established.

## What changed and why

The LR path retains an F reference even when the estimated variance is clamped to one; the former chi-square switch discarded variance uncertainty. The artificial minimum of 50 dispersion degrees of freedom is removed. The new bb_site option estimates beta-binomial rho and computes the variance of the pooled contrast using sums of squared sample depths. Group-specific depth corrections also enter the approximate effect intervals.

Max-T region permutations repeat candidate discovery, retain successful empty scans in the denominator, abort on failed scans, and use the +1 Monte Carlo correction. The already scan-adjusted p-values are not BH-adjusted again. Under exchangeability, this procedure targets complete-null family-wise control; strong control under partial alternatives is not proved.

The experimental bb_eb moment shrinkage option was evaluated and lost substantial signal. It is not the selected configuration. Neighbour p-value combining, smoothing, and Fisher separation fallback were disabled in these comparisons.

## Mechanism evidence

On the same 100 compact v4 null datasets, the original site-mode adaptive switch produced discoveries in 100/100 datasets; retaining F reduced this to 5/100. Depth-aware site variance produced 0/100. The EB phi model remained at 14/100 with or without the switch. These are development results. Known-dispersion controls support the variance-mismatch diagnosis, but the Pearson moment and F-reference derivations remain approximations.

A subsequent 500 fresh-seed compact null check had 2 false CpG discoveries in 2 datasets: FDR 0.40%, exact 95% interval 0.05%–1.44%. Replicates are independent conditional on the fixed calibration library. This does not certify arbitrary tissues or full genomes.

The additional chr1–3 null contains 3,571,573 eligible CpGs. It produced 0 significant CpGs and 0 significant regions. This single null replicate is a scale check, not an FDR estimate.

## Unchanged chr1–3 signal dataset

All rows use the same existing seed 4101 count files, common eligible CpGs, immutable truth and current strict matching. Comparator outputs are reused from the reviewed run; epykit is rerun. Historical epykit region ends are converted before rescoring. No competitor was retuned.

| method | DML_TP | DML_FP | DML_FDP | DML_AP |
| --- | --- | --- | --- | --- |
| epykit | 91 | 357 | 0.7969 | 0.04497 |
| DSS | 3911 | 19849 | 0.8354 | 0.3639 |
| methylKit | 412 | 3536 | 0.8956 | 0.04247 |
| BSmooth | NA | NA | NA | NA |
| dmrseq | NA | NA | NA | NA |
| DMRcate | 0 | 0 | 0 | 0.06044 |
| epykit depth-aware, asymptotic regions | 0 | 0 | 0 | 0.05475 |
| epykit experimental shrinkage | 0 | 0 | 0 | 0.01152 |
| epykit depth-aware, permutation regions | 0 | 0 | 0 | 0.05475 |

| method | DMR_matched | DMR_called | DMR_precision | DMR_recall | DMR_F1 | DMR_wholly_null |
| --- | --- | --- | --- | --- | --- | --- |
| epykit | 231 | 619 | 0.3732 | 0.3505 | 0.3615 | 341 |
| DSS | 172 | 287 | 0.5993 | 0.261 | 0.3636 | 14 |
| methylKit | 65 | 744 | 0.08737 | 0.09863 | 0.09266 | 620 |
| BSmooth | 81 | 204 | 0.3971 | 0.1229 | 0.1877 | 27 |
| dmrseq | 216 | 241 | 0.8963 | 0.3278 | 0.48 | 5 |
| DMRcate | 0 | 0 | NA | 0 | 0 | 0 |
| epykit depth-aware, asymptotic regions | 185 | 352 | 0.5256 | 0.2807 | 0.366 | 115 |
| epykit experimental shrinkage | 52 | 103 | 0.5049 | 0.07891 | 0.1365 | 16 |
| epykit depth-aware, permutation regions | 22 | 22 | 1 | 0.03338 | 0.06461 | 0 |

Strict precision counts boundary/fragmentation failures as unmatched; wholly-null calls contain no true positive eligible CpG. These are different errors. No-call precision is undefined. The corrected test returns no significant CpGs on this difficult signal seed: error control has improved, but CpG power remains inadequate.

## Computation

| method | seconds | peak_RSS_GiB |
| --- | --- | --- |
| epykit | 67 | 1.929 |
| DSS | 521.5 | 4.833 |
| methylKit | 1784 | 7.548 |
| BSmooth | 1600 | 12.49 |
| dmrseq | 5387 | 4.848 |
| DMRcate | 115.2 | 6.818 |
| epykit depth-aware, asymptotic regions | 64.12 | 1.887 |
| epykit experimental shrinkage | 57.09 | 1.912 |
| epykit depth-aware, permutation regions | 3329 | 2.143 |

Permutation inference includes 100 complete scan repetitions and must be compared with that cost included. Earlier development timings overlapped diagnostic computations and are not clean speed estimates. Existing competitor timings are historical measurements under the same one-thread protocol, not randomized repeated timings from this session.

## Limits and reproducibility

Max-T targets scan-level family-wise error, a stricter criterion than the native region FDR reported by dmrseq. Different rows therefore do not establish matched-error noninferiority. The signal seed was already inspected, and one region comparison cannot support a general ranking. Correlated and real-data validation remain outstanding. The heterogeneous BLUEPRINT calibration and common-mask filtering limit biological generalization.

Source snapshots, dataset hashes and runner configurations accompany each new result directory. See comparison.tsv for source paths, evidence.json for null evidence, and analysis/diagnose_calibration_v5.py, analysis/validate_depth_aware_v5.py and analysis/oracle_single_site.py for reproducible controls. Rebuild this document with PYTHONPATH=vendor310:src python3 analysis/compare_calibration_v5.py.
