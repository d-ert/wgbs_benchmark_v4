suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
root <- '/scratch/wgbs_benchmark_v4'
out <- file.path(root, 'results/beta_binomial_20261007_full/comparison')
paper <- fread(file.path(root, 'analysis/paper_dmr_comparison_20261006/paper_regions_with_matches.tsv'))
summary <- as.data.frame(fread(file.path(out, 'paper_comparison.csv')))
stopifnot(nrow(paper) == 3271L, all(paper$start >= 1L), all(paper$end >= paper$start))
pr <- GRanges(paper$chrom, IRanges(paper$start, paper$end))
for (tool in unique(summary$tool)) {
  calls <- fread(file.path(out, 'real', tool, 'dmr.tsv'))
  direction <- if ('mean_meth_diff' %in% names(calls)) {
    ifelse(calls$mean_meth_diff > 0, 'hyper', 'hypo')
  } else ifelse(calls$direction > 0, 'hyper', 'hypo')
  cr <- GRanges(calls$chrom, IRanges(calls$start, calls$end))
  hits <- findOverlaps(cr, pr, ignore.strand=TRUE)
  a <- queryHits(hits); b <- subjectHits(hits)
  overlap <- pmin(calls$end[a], paper$end[b]) - pmax(calls$start[a], paper$start[b]) + 1L
  same <- direction[a] == paper$direction[b]
  for (cutoff in c(0, .5, .8)) {
    eligible <- overlap/(calls$end[a]-calls$start[a]+1) >= cutoff &
                overlap/(paper$end[b]-paper$start[b]+1) >= cutoff
    kept <- eligible & same
    row <- summary[summary$tool == tool & summary$reciprocal_threshold == cutoff, ]
    stopifnot(nrow(row)==1L,
              row$calls_overlapping_paper == length(unique(a[kept])),
              row$paper_regions_covered == length(unique(b[kept])),
              row$discordant_pairs == sum(eligible & !same))
  }
}
cat('Independent GenomicRanges paper-overlap counts passed for all four call sets and all thresholds.\n')
