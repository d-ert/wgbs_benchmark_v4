#!/usr/bin/env Rscript
# Descriptive hg19 region/gene concordance across the DEFAULT tool outputs.
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: Rscript annotate_real_dmrs.R REAL_RESULTS OUT_DIR TOOL1,TOOL2,...")
root <- args[[1]]; out <- args[[2]]; tools <- strsplit(args[[3]], ",", fixed=TRUE)[[1]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (length(tools)<2L || anyDuplicated(tools)) stop("Need at least two distinct tools")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
gene_gr <- genes(txdb)
gene_id <- names(gene_gr)
gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=gene_id,
                                     keytype="ENTREZID", column="SYMBOL", multiVals="first")
prom <- promoters(gene_gr, upstream=2000L, downstream=500L)
names(prom) <- gene_id
empty <- function() data.table(chrom=character(),start=integer(),end=integer(),direction=numeric())
read_calls <- function(tool) {
  path <- file.path(root,tool,"dmr.tsv")
  if (!file.exists(path)) stop("Missing default DMR list: ",path)
  d <- fread(path)
  if (!nrow(d)) return(empty())
  if (!all(c("chrom","start","end") %in% names(d))) stop("Malformed DMR list: ",path)
  if (!"direction" %in% names(d)) {
    effect_col <- intersect(c("mean_meth_diff", "meth_diff"), names(d))
    d[,direction:=if (length(effect_col)) sign(as.numeric(d[[effect_col[[1]]]])) else NA_real_]
  }
  d <- d[, .(chrom=as.character(chrom),start=as.integer(start),end=as.integer(end),
              direction=if ("direction" %in% names(d)) as.numeric(direction) else rep(NA_real_,.N))]
  if (anyNA(d[,.(chrom,start,end)]) || any(d$start<1L | d$end<d$start)) stop("Invalid 1-based regions: ",path)
  d[,chrom:=ifelse(grepl("^chr",chrom),chrom,paste0("chr",chrom))]
  d[,dmr_id:=paste0(tool,"_",.I)]
  d
}
calls <- setNames(lapply(tools,read_calls),tools)
to_gr <- function(d) GRanges(seqnames=d$chrom,ranges=IRanges(start=d$start,end=d$end))
all_annot <- list(); gene_sets <- list()
for (tool in tools) {
  d <- calls[[tool]]
  if (!nrow(d)) {
    all_annot[[tool]] <- data.table(tool=character(),dmr_id=character(),chrom=character(),
      start=integer(),end=integer(),direction=numeric(),gene_id=character(),
      symbol=character(),relation=character())
    gene_sets[[tool]] <- character()
    next
  }
  g <- to_gr(d)
  prom_hit <- findOverlaps(g,prom,ignore.strand=TRUE)
  body_hit <- findOverlaps(g,gene_gr,ignore.strand=TRUE)
  make <- function(hit,relation) {
    if (!length(hit)) return(NULL)
    ii <- queryHits(hit);jj <- subjectHits(hit)
    data.table(tool=tool,dmr_id=d$dmr_id[ii],chrom=d$chrom[ii],start=d$start[ii],
      end=d$end[ii],direction=d$direction[ii],gene_id=gene_id[jj],
      symbol=unname(gene_symbol[jj]),relation=relation)
  }
  tab <- unique(rbindlist(list(make(prom_hit,"promoter"),make(body_hit,"gene_body")),fill=TRUE))
  if (!nrow(tab)) tab <- data.table(tool=tool,dmr_id=d$dmr_id,chrom=d$chrom,
    start=d$start,end=d$end,direction=d$direction,gene_id=NA_character_,
    symbol=NA_character_,relation="intergenic")
  missing <- setdiff(d$dmr_id,tab$dmr_id)
  if (length(missing)) {
    x <- d[match(missing,d$dmr_id)]
    tab <- rbind(tab,data.table(tool=tool,dmr_id=x$dmr_id,chrom=x$chrom,
      start=x$start,end=x$end,direction=x$direction,gene_id=NA_character_,
      symbol=NA_character_,relation="intergenic"),fill=TRUE)
  }
  all_annot[[tool]] <- tab
  gene_sets[[tool]] <- unique(na.omit(tab$gene_id))
}
annotation <- rbindlist(all_annot,use.names=TRUE,fill=TRUE)
fwrite(annotation,file.path(out,"dmr_gene_annotations.tsv"),sep="\t")
universe <- sort(unique(unlist(gene_sets)))
support <- data.table(gene_id=universe,symbol=unname(gene_symbol[match(universe,gene_id)]))
for (tool in tools) support[[tool]] <- universe %in% gene_sets[[tool]]
support[,n_tools:=rowSums(.SD),.SDcols=tools]
fwrite(support,file.path(out,"gene_support.tsv"),sep="\t")
pair_rows <- list();gene_rows <- list();idx <- 0L
for (i in seq_len(length(tools)-1L)) for (j in (i+1L):length(tools)) {
  a <- tools[[i]];b <- tools[[j]];da <- calls[[a]];db <- calls[[b]]
  ga <- gene_sets[[a]];gb <- gene_sets[[b]]
  inter_genes <- length(intersect(ga,gb));union_genes <- length(union(ga,gb))
  idx <- idx+1L
  gene_rows[[idx]] <- data.table(tool_a=a,tool_b=b,genes_a=length(ga),genes_b=length(gb),
    shared_genes=inter_genes,gene_jaccard=if (union_genes) inter_genes/union_genes else NA_real_)
  match_n <- 0L;agree <- 0L;known <- 0L
  if (nrow(da) && nrow(db)) {
    ha <- to_gr(da);hb <- to_gr(db);hit <- findOverlaps(ha,hb,ignore.strand=TRUE)
    if (length(hit)) {
      qa <- queryHits(hit);qb <- subjectHits(hit)
      overlap <- pmax(0L,pmin(da$end[qa],db$end[qb])-pmax(da$start[qa],db$start[qb])+1L)
      wa <- da$end[qa]-da$start[qa]+1L;wb <- db$end[qb]-db$start[qb]+1L
      keep <- overlap/wa>=.5 & overlap/wb>=.5
      candidates <- data.table(qa=qa[keep],qb=qb[keep],score=overlap[keep]/(wa[keep]+wb[keep]-overlap[keep]))
      if (nrow(candidates)) {
        setorder(candidates,-score,qa,qb)
        used_a <- rep(FALSE,nrow(da));used_b <- rep(FALSE,nrow(db))
        accepted <- logical(nrow(candidates))
        for (k in seq_len(nrow(candidates))) {
          x <- candidates$qa[k];y <- candidates$qb[k]
          if (!used_a[x] && !used_b[y]) {accepted[k] <- TRUE;used_a[x] <- TRUE;used_b[y] <- TRUE}
        }
        selected <- candidates[accepted];match_n <- nrow(selected)
        x <- da$direction[selected$qa];y <- db$direction[selected$qb]
        valid <- is.finite(x) & is.finite(y) & x!=0 & y!=0
        known <- sum(valid);agree <- sum(sign(x[valid])==sign(y[valid]))
      }
    }
  }
  pair_rows[[idx]] <- data.table(tool_a=a,tool_b=b,dmrs_a=nrow(da),dmrs_b=nrow(db),
    reciprocal_50pct_matches=match_n,
    fraction_a_matched=if (nrow(da)) match_n/nrow(da) else NA_real_,
    fraction_b_matched=if (nrow(db)) match_n/nrow(db) else NA_real_,
    matched_direction_known=known,
    matched_direction_concordance=if (known) agree/known else NA_real_)
}
fwrite(rbindlist(pair_rows),file.path(out,"dmr_concordance.tsv"),sep="\t")
fwrite(rbindlist(gene_rows),file.path(out,"gene_concordance.tsv"),sep="\t")
fwrite(data.table(tool=tools,dmr_count=vapply(calls,nrow,integer(1)),
                  gene_count=vapply(gene_sets,length,integer(1))),
       file.path(out,"tool_counts.tsv"),sep="\t")
writeLines(c("Genome build: hg19/GRCh37; gene model: TxDb.Hsapiens.UCSC.hg19.knownGene.",
  "Promoter: 2 kb upstream and 500 bp downstream of each gene TSS; gene body uses gene boundaries.",
  "DMR matches: >=50% reciprocal base-pair overlap, greedy one-to-one by interval Jaccard.",
  "Gene concordance: Jaccard of Entrez gene IDs touched by promoter or gene body DMRs.",
  "These descriptive overlaps do not estimate biological truth, FDR, or independent replication.",
  "Only default tool directories are compared; exploratory epykit profiles are excluded."),
  file.path(out,"methods.txt"))
sink(file.path(out,"session_info.txt"));print(sessionInfo());sink()
