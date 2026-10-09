# Paper Table S2 versus real-data epykit and DSS

DSS agrees more closely with the reported regions under all tested overlap thresholds. This is concordance with a published call set, not a measurement of biological precision or FDR.

| Measure | epykit | DSS |
|---|---:|---:|
| Total calls | 431 | 253 |
| Calls overlapping any paper DMR, same direction | 158 (36.7%) | 180 (71.1%) |
| Distinct paper regions covered, same direction | 160 / 3,271 (4.9%) | 180 / 3,271 (5.5%) |
| ≥50% reciprocal overlap, same direction | 99 (23.0% of calls) | 130 (51.4% of calls) |
| ≥80% reciprocal overlap, same direction | 30 (7.0% of calls) | 60 (23.7% of calls) |
| Calls with no paper overlap | 273 | 73 |

Every observed overlapping pair agrees in methylation direction. Any-overlap counts permit one region to overlap several others; this explains 158 epykit calls covering 160 paper regions. Separate maximum-cardinality one-to-one match counts are available in the workbook. At 50% and 80% reciprocal overlap the displayed counts are also one-to-one matches.

## Direction-specific comparison

The workbook contains 1,714 hypo- and 1,557 hypermethylated paper regions.

- Epykit: 149/298 hypo calls overlap 151 paper hypo regions (8.8% of the paper hypo set). Nine of 133 hyper calls overlap paper hyper regions.
- DSS: 177/234 hypo calls overlap 177 paper hypo regions (10.3% of the paper hypo set). Three of 19 hyper calls overlap paper hyper regions.
- Using any overlap, 87 published regions are recovered by both tools, 73 by epykit alone and 93 by DSS alone: 253 published regions in the union. The corresponding hypo counts are 86 / 65 / 91, a union of 242 hypo regions.
- At ≥50% reciprocal overlap, both tools recover 53 published regions; epykit alone recovers 46 and DSS alone 77.

## Interpretation limits

The workbook explicitly warns that many reported hypermethylated regions may be false positives. The paper tested a small independent subset and validated the tested hypomethylated sites but not the tested hypermethylated sites; this does not individually validate all 1,714 hypo regions. Hypomethylated-region concordance is therefore the more useful comparison, while still not a ground-truth accuracy test.

The paper used BSmooth and regions with at least three CpGs. Epykit's current runner requires at least five CpGs and 50 bp, while DSS requires three CpGs and 50 bp. Of the 3,271 reported regions, 2,479 contain fewer than five CpGs and 746 are shorter than 50 bp; only 733 meet both ≥5 reported CpGs and ≥50 bp. These are structural differences, not a proof that the excluded short paper regions could never overlap a larger caller region. Paper CpG counts also need not equal the current common-coverage counts.

The workbook's log2FC and FDR columns describe associated-gene expression; they were neither used for methylation direction nor used to filter DMRs. These columns are renamed in the normalized output for clarity. No FDR threshold was applied to the reported reference list.

Current real-data runners omit the donor-pair term. Low paper overlap can reflect method, threshold, boundaries and input filtering differences. Unmatched calls must not automatically be labelled false positives.

## Inputs, methods and reproduction

- Input workbook: `/scratch/wgbs_benchmark_v4/data/real/paper_reported_dmrs.xlsx`, Sheet1, header on Excel row 3. Four trailing blank/footnote rows are excluded. All 3,271 genomic rows are retained without coordinate duplicates.
- Calls: `results/genome_baseline/real/{epykit,DSS}/dmr.tsv`.
- Genome: GRCh37/hg19 for paper and benchmark. Paper reported lengths exactly equal `end-start+1`; comparisons use 1-based inclusive intervals.
- Any overlap requires ≥1 shared base pair on the same chromosome. Reciprocal thresholds require the overlap to cover the specified fraction of **both** intervals. Counts require matching hypo/hyper direction, defined as MTB relative to NI.
- Primary article: https://pmc.ncbi.nlm.nih.gov/articles/PMC4665002/ (methods and Table S2 description).
- Reproduce: `python3 analysis/paper_dmr_comparison_20261006/compare_paper.py` from the v4 directory, or invoke the script by absolute path.
- `validation.txt` records an independent GenomicRanges cross-check of pair counts and unique interval counts for both tools at all four overlap thresholds.

## Files

`paper_vs_epykit_DSS.xlsx` contains summary tables, direction-specific results, every paper region with caller match flags, every caller region with matched paper gene labels and nearest-paper distances, all overlapping pairs, and one-to-one matches. `excel_row` links each paper region back to the original workbook. Paper gene labels are the workbook's annotations, not newly inferred gene assignments.

TSV/CSV counterparts support reuse in analysis. `provenance.json` records input hashes, methods, exclusions and source URL. The raw benchmark and original workbook are unchanged.
