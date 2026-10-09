#!/usr/bin/env Rscript
# One frozen, count-identical runner for the five R tool families.
options(warn=1) # Preserve every warning in the caller log, including failures.
suppressPackageStartupMessages({
  library(data.table)
  library(bsseq)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Usage: Rscript run_r_tool.R TOOL DATASET OUTPUT_DIR")
tool <- args[[1]]; dataset <- args[[2]]; out <- args[[3]]
if (!tool %in% c("DSS", "methylKit", "BSmooth", "dmrseq", "DMRcate")) stop("Unknown tool")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(20260929L)
data.table::setDTthreads(1L)
sheet <- fread(file.path(dataset, "samples.tsv"))
assembly <- if (file.exists(file.path(dataset, "manifest.json"))) {
  manifest <- jsonlite::fromJSON(file.path(dataset, "manifest.json"))
  if (!is.null(manifest$assembly)) manifest$assembly else "hg38"
} else "hg38"
stopifnot(identical(sort(unique(sheet$group)), c("control", "treatment")))
stopifnot(nrow(sheet) >= 4L, all(file.exists(sheet$path)))
setorder(sheet, group, sample_id)
elapsed_start <- proc.time()[["elapsed"]]
phase <- list()
mark <- function(name) { phase[[name]] <<- proc.time()[["elapsed"]] - elapsed_start }
empty_regions <- function() data.frame(chrom=character(), start=integer(), end=integer(),
                                       pvalue=numeric(), qvalue=numeric(), direction=integer())
write_tab <- function(x, name) fwrite(x, file.path(out, name), sep="\t")
write_tab(empty_regions(), "dmr.tsv")

if (tool == "methylKit") {
  suppressPackageStartupMessages(library(methylKit))
  overdispersion <- Sys.getenv("WGBS_METHYLKIT_OVERDISPERSION", "MN")
  adjustment <- Sys.getenv("WGBS_METHYLKIT_ADJUST", "BH")
  message("methylKit configuration: overdispersion=", overdispersion,
          ", adjust=", adjustment)
  obj <- methRead(as.list(sheet$path), sample.id=as.list(sheet$sample_id),
                  treatment=as.integer(sheet$group == "treatment"),
                  assembly=assembly, pipeline="bismarkCoverage", context="CpG", mincov=1L)
  united <- unite(obj, destrand=FALSE)
  mark("read_unite")
  rm(obj); gc(FALSE) # Release inputs retained by neither downstream stage.
  result <- calculateDiffMeth(united, overdispersion=overdispersion,
                              adjust=adjustment, mc.cores=1)
  dt <- as.data.table(getData(result))
  write_tab(data.table(chrom=as.character(dt$chr), pos=as.integer(dt$start),
                       pvalue=dt$pvalue, qvalue=dt$qvalue,
                       meth_diff=dt$`meth.diff` / 100), "dml.tsv")
  mark("dml")
  rm(result, dt); gc(FALSE) # CpG results are exported; the window stage uses united.
  # The native fixed-window method is reported separately from de novo callers.
  tiles <- tileMethylCounts(united, win.size=1000, step.size=1000)
  tiled <- calculateDiffMeth(tiles, overdispersion=overdispersion,
                             adjust=adjustment, mc.cores=1)
  rd <- as.data.table(getData(tiled))
  regions <- data.table(chrom=as.character(rd$chr), start=as.integer(rd$start),
                        end=as.integer(rd$end), pvalue=rd$pvalue,
                        qvalue=rd$qvalue, direction=sign(rd$`meth.diff`))
  write_tab(regions, "dmr_candidates.tsv")
  write_tab(regions[is.finite(qvalue) & qvalue <= .05], "dmr.tsv")
  mark("tiles")
} else {
  tables <- lapply(sheet$path, function(file) {
    d <- fread(file, header=FALSE, select=c(1L, 2L, 5L, 6L),
               col.names=c("chrom", "pos", "m", "u"))
    setkey(d, chrom, pos)
    if (anyDuplicated(d, by=c("chrom", "pos"))) stop("Duplicate CpG in ", file)
    d
  })
  keys <- tables[[1]][, .(chrom, pos)]
  for (i in 2:length(tables)) keys <- tables[[i]][keys, on=.(chrom, pos), nomatch=0L][, .(chrom, pos)]
  setkey(keys, chrom, pos)
  if (!nrow(keys)) stop("No CpG has coverage in every sample")
  M <- matrix(0L, nrow(keys), nrow(sheet)); Cov <- M
  for (i in seq_along(tables)) {
    d <- tables[[i]][keys, on=.(chrom, pos)]
    if (anyNA(d$m) || anyNA(d$u)) stop("Common CpG join failed")
    M[, i] <- d$m; Cov[, i] <- d$m + d$u
  }
  rm(tables); gc(FALSE)
  if (any(Cov < 1L) || any(M > Cov)) stop("Invalid common count matrix")
  bs <- BSseq(chr=keys$chrom, pos=keys$pos, M=M, Cov=Cov, sampleNames=sheet$sample_id)
  pData(bs)$group <- factor(sheet$group, levels=c("control", "treatment"))
  mark("read_filter")
  if (tool == "DSS") {
    suppressPackageStartupMessages(library(DSS))
    dss_smoothing <- identical(Sys.getenv("WGBS_DSS_SMOOTHING", "true"), "true")
    dss_span <- as.numeric(Sys.getenv("WGBS_DSS_SMOOTHING_SPAN", "500"))
    dss_p_threshold <- as.numeric(Sys.getenv("WGBS_DSS_P_THRESHOLD", "1e-5"))
    if (!is.finite(dss_span) || dss_span <= 0 || !is.finite(dss_p_threshold) ||
        dss_p_threshold <= 0 || dss_p_threshold >= 1) stop("Invalid DSS parameters")
    message("DSS configuration: smoothing=", dss_smoothing,
            ", span=", dss_span, ", callDMR p.threshold=", dss_p_threshold)
    test <- DMLtest(bs, group1=sheet$sample_id[sheet$group == "treatment"],
                    group2=sheet$sample_id[sheet$group == "control"],
                    equal.disp=FALSE, smoothing=dss_smoothing,
                    smoothing.span=dss_span, ncores=1L)
    message("DSS DML columns: ", paste(names(test), collapse=", "))
    pv <- if ("pvals" %in% names(test)) test$pvals else test$pval
    if (is.null(pv)) stop("DSS produced no p-value column")
    q <- p.adjust(pv, method="BH")
    write_tab(data.table(chrom=as.character(test$chr), pos=as.integer(test$pos),
                         pvalue=pv, qvalue=q, meth_diff=test$diff), "dml.tsv")
    mark("dml")
    regions <- callDMR(test, delta=0, p.threshold=dss_p_threshold, minlen=50,
                       minCG=3, dis.merge=100, pct.sig=0.5)
    if (!is.null(regions) && nrow(regions)) {
      rd <- as.data.table(regions)
      message("DSS DMR columns: ", paste(names(rd), collapse=", "))
      effect_col <- intersect(c("diff.Methy", "meanDiff"), names(rd))
      if (!length(effect_col)) stop("DSS DMR output lacks an effect direction column")
      region_effect <- as.numeric(rd[[effect_col[1]]])
      region_p <- if ("pval" %in% names(rd)) rd$pval else rep(NA_real_, nrow(rd))
      write_tab(data.table(chrom=as.character(rd$chr), start=as.integer(rd$start),
                           end=as.integer(rd$end), pvalue=region_p,
                           qvalue=NA_real_, direction=sign(region_effect)), "dmr.tsv")
    }
    mark("dmr")
  } else if (tool == "BSmooth") {
    threshold <- Sys.getenv("WGBS_BSMOOTH_THRESHOLD", "fixed")
    if (!threshold %in% c("fixed", "quantile")) stop("Invalid BSmooth threshold mode")
    smooth_ns <- as.integer(Sys.getenv("WGBS_BSMOOTH_NS", "70"))
    smooth_h <- as.numeric(Sys.getenv("WGBS_BSMOOTH_H", "1000"))
    fixed_cutoff <- as.numeric(Sys.getenv("WGBS_BSMOOTH_FIXED_CUTOFF", "4.6"))
    if (is.na(smooth_ns) || smooth_ns < 2 || !is.finite(smooth_h) || smooth_h < 1)
      stop("Invalid BSmooth smoothing parameters")
    if (!is.finite(fixed_cutoff) || fixed_cutoff <= 0) stop("Invalid BSmooth fixed cutoff")
    message("BSmooth region threshold: ", threshold, "; ns: ", smooth_ns,
            "; h: ", smooth_h)
    smooth <- BSmooth(bs, h=smooth_h, ns=smooth_ns,
                     BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    stat <- BSmooth.tstat(smooth,
                         group1=sheet$sample_id[sheet$group == "treatment"],
                         group2=sheet$sample_id[sheet$group == "control"],
                         estimate.var="same", local.correct=TRUE, verbose=FALSE)
    write_tab(data.table(chrom=as.character(seqnames(stat)), pos=as.integer(start(stat)),
                         score=abs(as.numeric(getStats(stat)[, "tstat.corrected"]))), "cpg_scores.tsv")
    mark("smooth_tstat")
    if (threshold == "quantile") {
      regions <- dmrFinder(stat, cutoff=NULL, qcutoff=c(.025, .975),
                           maxGap=300, verbose=FALSE)
      if (!is.null(regions) && nrow(regions)) {
        regions <- subset(regions, abs(meanDiff) > .1 & n >= 3)
      }
    } else {
      regions <- dmrFinder(stat, cutoff=c(-fixed_cutoff, fixed_cutoff), qcutoff=NULL,
                           maxGap=300, verbose=FALSE)
    }
    if (!is.null(regions) && nrow(regions)) {
      regions <- subset(regions, n >= 3 & abs(meanDiff) >= .1)
      rd <- as.data.table(regions)
      write_tab(data.table(chrom=as.character(rd$chr), start=as.integer(rd$start),
                           end=as.integer(rd$end), pvalue=NA_real_,
                           qvalue=NA_real_, direction=sign(rd$meanDiff)), "dmr.tsv")
    }
    mark("dmr")
  } else if (tool == "dmrseq") {
    suppressPackageStartupMessages({library(BiocParallel); register(SerialParam(), default=TRUE); library(dmrseq)})
    regions <- dmrseq(bs=bs, testCovariate="group", cutoff=.1,
                      minNumRegion=5, maxPerms=100,
                      BPPARAM=SerialParam(), verbose=FALSE)
    if (length(regions)) {
      rd <- as.data.table(as.data.frame(regions))
      pv <- if ("pval" %in% names(rd)) rd$pval else rd$pvalue
      qv <- if ("qval" %in% names(rd)) rd$qval else rd$qvalue
      native <- data.table(chrom=as.character(rd$seqnames), start=as.integer(rd$start),
                           end=as.integer(rd$end), pvalue=pv,
                           qvalue=qv, direction=sign(rd$beta))
      write_tab(native, "dmr_candidates.tsv")
      write_tab(native[is.finite(qvalue) & qvalue <= .05], "dmr.tsv")
    }
    mark("dmr")
  } else if (tool == "DMRcate") {
    suppressPackageStartupMessages(library(DMRcate))
    initial_fdr <- as.numeric(Sys.getenv("WGBS_DMRCATE_INITIAL_FDR", "0.05"))
    kernel_lambda <- as.numeric(Sys.getenv("WGBS_DMRCATE_LAMBDA", "1000"))
    if (!is.finite(initial_fdr) || initial_fdr <= 0 || initial_fdr >= 1) stop("Invalid DMRcate initial FDR")
    if (!is.finite(kernel_lambda) || kernel_lambda < 1) stop("Invalid DMRcate kernel lambda")
    cutoff_text <- Sys.getenv("WGBS_DMRCATE_PCUTOFF", "fdr")
    region_cutoff <- if (cutoff_text == "fdr") "fdr" else as.numeric(cutoff_text)
    if (cutoff_text != "fdr" && (!is.finite(region_cutoff) || region_cutoff <= 0 || region_cutoff >= 1))
      stop("Invalid DMRcate region cutoff")
    message("DMRcate initial CpG FDR: ", initial_fdr, "; region pcutoff: ", cutoff_text,
            "; lambda: ", kernel_lambda)
    design <- model.matrix(~ group, data=pData(bs))
    methdesign <- edgeR::modelMatrixMeth(design)
    annotated <- sequencing.annotate(bs, methdesign=methdesign, all.cov=FALSE,
                                     fdr=initial_fdr, coef=tail(colnames(methdesign), 1))
    ar <- annotated@ranges
    write_tab(data.table(chrom=as.character(seqnames(ar)), pos=as.integer(start(ar)),
                         pvalue=ar$rawpval, qvalue=ar$ind.fdr), "dml.tsv")
    regions <- tryCatch(dmrcate(annotated, lambda=kernel_lambda, C=2,
                                pcutoff=region_cutoff, min.cpgs=2),
                        error=function(e) {
                          if (grepl("returned no significant CpGs", conditionMessage(e), fixed=TRUE)) {
                            message("DMRcate: no significant CpGs; zero DMR discoveries")
                            return(NULL)
                          }
                          stop(e)
                        })
    if (!is.null(regions)) {
      message("DMRcate candidate regions: ", length(regions@coord),
              "; min smoothed FDR <=0.05: ", sum(is.finite(regions@min_smoothed_fdr) & regions@min_smoothed_fdr <= .05))
      # extractRanges() additionally invokes ExperimentHub for gene annotation.
      # Genomic coordinates and scores are already present in DMResults.
      coords <- as.data.table(DMRcate:::extractCoords(regions@coord))
      if (nrow(coords)) {
        native <- data.table(chrom=as.character(coords$chrom),
                             start=as.integer(as.character(coords$chromStart)),
                             end=as.integer(as.character(coords$chromEnd)),
                             pvalue=NA_real_, qvalue=regions@min_smoothed_fdr,
                             direction=sign(regions@meandiff))
        write_tab(native, "dmr_candidates.tsv")
      write_tab(native[is.finite(qvalue) & qvalue <= .05], "dmr.tsv")
      }
    }
    mark("dmr")
  }
}
write_tab(data.table(phase=names(phase), cumulative_wall_seconds=unlist(phase)), "timing.tsv")
sink(file.path(out, "session_info.txt")); print(sessionInfo()); sink()
