# Investigation sequence

Results reviewed first: complete autosomal signal/null six-tool outputs, region-only correction, experimental BB-F, development CpG strata, held-out BB validation, and completed uncapped GSE64177 real comparisons. The original overall genome run failed during real methylKit; its completed simulations remain usable. GSE263850 full manifest is running with DSS active; smoke is not a biological benchmark.

Literature searched next: DSS official guide; BSmooth original article; methylSig original article/manual; RADMeth original article; dmrseq original article/repository; Benjamini–Yekutieli original paper. These motivate the following hypotheses before detailed code tracing.

1. Original region false calls: selecting boundaries on the same CpG evidence and correcting only selected candidates understates search multiplicity. Status to test: compare frozen original and corrected source, reproduce correction and unchanged CpG results.
2. CpG null calibration varies with depth because a quasi-binomial scale is shrunk toward one chromosome mean although beta-binomial biological variation grows with depth. Status to test: inspect sufficient statistics/dispersion estimator and evaluate simple controlled count cases.
3. Small-sample F tails and lack of spatial borrowing can explain low site power despite acceptable ranking, especially after 16 million tests. Status to test: inspect actual reference branch, raw p/q distributions, theoretical threshold examples, and BB conditional nuisance handling. Do not assume stale comments describe behavior.
4. Region-only correction loses recall because complete-interval BY is conservative and seed/geometry rules miss weak, short or sparse signals. Status to test: inspect family counting, filtering and rejected-candidate stages, quantify geometry where feasible.
5. BB runtime rises because each site nests many dispersion searches and repeated group mean fits, rather than because of larger memory. Status to test: inspect compiled loops and run a bounded profile outside production code.
6. Real-paper recovery is reduced by minimum CpG/length rules, different thresholds and unpaired models. Status to test: audit real harness, paper morphology and matching definitions. Paper agreement is not truth.

No changes to the package, benchmark inputs or scoring are needed for this reporting task. Diagnostic scripts and report artifacts are isolated here.
