suppressPackageStartupMessages({library(data.table);library(GenomicRanges)})
root <- 'analysis/paper_dmr_comparison_20261006'
p <- fread(file.path(root,'paper_regions_with_matches.tsv'))
gp <- GRanges(p$chrom,IRanges(p$start,p$end))
s <- fread(file.path(root,'summary.csv'))
for (tool_name in c('epykit','DSS')) {
 calls <- fread(file.path(root,paste0(tool_name,'_calls_vs_paper.tsv')))
 gcalls <- GRanges(calls$chrom,IRanges(calls$start,calls$end))
 h <- findOverlaps(gcalls,gp,ignore.strand=TRUE)
 i <- queryHits(h);j <- subjectHits(h)
 ov <- pmin(calls$end[i],p$end[j])-pmax(calls$start[i],p$start[j])+1L
 for (cutoff in c(0,.25,.5,.8)) {
  keep <- ov/(calls$end[i]-calls$start[i]+1)>=cutoff & ov/(p$end[j]-p$start[j]+1)>=cutoff & calls$direction[i]==p$direction[j]
  rule_name <- if (cutoff==0) 'any_overlap' else paste0('reciprocal_',cutoff*100,'pct')
  index <- which(s[['tool']]==tool_name & s[['rule']]==rule_name)
  z <- s[index]
  stopifnot(nrow(z)==1,length(unique(i[keep]))==z$concordant_calls,length(unique(j[keep]))==z$concordant_paper,sum(keep)==z$concordant_pairs)
 }
}
writeLines('PASS: Independent GenomicRanges overlap counts agree for both tools and all four thresholds.',file.path(root,'validation.txt'))
cat('Independent overlap validation passed.\n')
