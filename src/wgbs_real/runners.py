import json
from pathlib import Path

PAPER_SETTINGS={
    'DSS':dict(model='DMLfit.multiFactor + DMLtest.multiFactor',formula='~group',smoothing=True,
        smoothing_span=500,p_threshold=1e-5,delta=0,minlen=50,minCG=3,dis_merge_requested=100,
        dis_merge_effective=50,pct_sig=.5,region_significance='native raw-p heuristic'),
    'epykit':dict(test='lr',dispersion='eb',smoothing=True,smoothing_span=500,raw_seed_alpha=1e-5,
        min_abs_meth_diff=0,minlen=51,min_cpgs=4,dis_merge=100,pct_sig=.5,
        region_significance='accepted complete-interval-family BY q<=0.05 retained'),
    'methylKit':dict(test='native MN overdispersion',CpG_adjust='BH',raw_seed_alpha=1e-5,
        smoothing=False,minlen_strict=50,min_cpgs_strict=3,dis_merge=100,pct_sig=.5,
        region_significance='raw-p CpG aggregation; no nominal region FDR'),
    'BSmooth':dict(smoothing=True,h=500,ns=70,maxGap=100,min_cpgs=4,minlen=51,
        min_abs_meth_diff=0,score_cutoff=4.417173413469023,
        difference='normal-score equivalent of raw p=1e-5; native statistic has no calibrated CpG p'),
    'dmrseq':dict(smooth=True,bpSpan=500,maxGap=100,minNumRegion=4,cutoff=.01,maxPerms=100,
        region_q=.05,difference='zero candidate cutoff is rejected by native API; 1% candidate cutoff, no raw CpG p gate'),
    'DMRcate':dict(raw_CpG_gate=1e-5,lambda_bp=100,C=2,min_cpgs=4,minlen=51,
        native_region_q=.05,difference='raw-p candidate count drives native smoothed-rank threshold; lambda controls both merging and smoothing')}

def replace_once(text,old,new):
    if text.count(old)!=1:raise ValueError(f'Unexpected runner source around {old[:60]!r}')
    return text.replace(old,new)

def build_r_runner(original,target):
    text=Path(original).read_text()
    text=replace_once(text,'set.seed(20260929L)',
        'paper_aligned <- identical(Sys.getenv("WGBS_REAL_PROFILE", "paper_aligned"), "paper_aligned")\nset.seed(20260929L)')
    text=replace_once(text,'  united <- unite(obj, destrand=FALSE)',
        '  united <- unite(obj, destrand=FALSE)\n  rm(obj); gc(FALSE)')
    first=text.index('  # The native fixed-window')
    last=text.index('  mark("tiles")',first)+len('  mark("tiles")')
    block=text[first:last]
    text=text[:first]+'  if (!paper_aligned) {\n'+block+'\n  } else {\n    message("methylKit: raw-p CpG regions are aggregated by the pipeline adapter")\n  }\n'+text[last:]
    first=text.index('    test <- DMLtest(bs,')
    last=text.index('    message("DSS DML columns:',first)
    original_fit=text[first:last]
    replacement='''    if (paper_aligned) {
      design <- data.frame(group=pData(bs)$group)
      fit <- DMLfit.multiFactor(bs, design=design, formula=~group,
                               smoothing=TRUE, smoothing.span=500)
      test <- DMLtest.multiFactor(fit, coef=2)
      raw <- getMeth(bs, type="raw")
      test$diff <- rowMeans(raw[, sheet$group=="treatment", drop=FALSE]) -
                   rowMeans(raw[, sheet$group=="control", drop=FALSE])
      rm(fit, raw); gc(FALSE)
      message("DSS: multi-factor ~group; smoothing=TRUE; span=500; raw p=1e-5; delta=0; requested merge=100, native effective merge=50")
    } else {
'''+original_fit+'''    }
'''
    text=text[:first]+replacement+text[last:]
    text=replace_once(text,'      if (!length(effect_col)) stop("DSS DMR output lacks an effect direction column")\n      region_effect <- as.numeric(rd[[effect_col[1]]])',
        '      if (inherits(test, "DMLtest.multiFactor")) {\n        region_effect <- as.numeric(rd$areaStat)\n      } else {\n        if (!length(effect_col)) stop("DSS DMR output lacks an effect direction column")\n        region_effect <- as.numeric(rd[[effect_col[1]]])\n      }')
    text=replace_once(text,'    smooth <- BSmooth(bs, h=smooth_h, ns=smooth_ns,',
        '    if (paper_aligned) { smooth_h <- 500; fixed_cutoff <- qnorm(1-1e-5/2) }\n    smooth <- BSmooth(bs, h=smooth_h, ns=smooth_ns,')
    text=text.replace('maxGap=300','maxGap=if (paper_aligned) 100 else 300')
    text=replace_once(text,'      regions <- subset(regions, n >= 3 & abs(meanDiff) >= .1)',
        '      regions <- if (paper_aligned) subset(regions, n > 3 & end-start+1 > 50) else subset(regions, n >= 3 & abs(meanDiff) >= .1)')
    text=replace_once(text,'    regions <- dmrseq(bs=bs, testCovariate="group", cutoff=.1,\n                      minNumRegion=5, maxPerms=100,',
        '    regions <- dmrseq(bs=bs, testCovariate="group", cutoff=if (paper_aligned) .01 else .1,\n                      minNumRegion=if (paper_aligned) 4 else 5, maxPerms=100,\n                      bpSpan=if (paper_aligned) 500 else 1000, maxGap=if (paper_aligned) 100 else 1000,')
    text=replace_once(text,'    design <- model.matrix(~ group, data=pData(bs))',
        '    if (paper_aligned) kernel_lambda <- 100\n    design <- model.matrix(~ group, data=pData(bs))')
    text=replace_once(text,'    ar <- annotated@ranges',
        '    if (paper_aligned) annotated@ranges$is.sig <- is.finite(annotated@ranges$rawpval) & annotated@ranges$rawpval <= 1e-5\n    ar <- annotated@ranges')
    text=text.replace('pcutoff=region_cutoff, min.cpgs=2','pcutoff=region_cutoff, min.cpgs=if (paper_aligned) 4 else 2')
    text=replace_once(text,'if (grepl("returned no significant CpGs", conditionMessage(e), fixed=TRUE))',
        'if (grepl("returned no significant CpGs|No signficant regions found", conditionMessage(e)))')
    old='      write_tab(native[is.finite(qvalue) & qvalue <= .05], "dmr.tsv")'
    if text.count(old)!=2:raise ValueError('Unexpected native region export branches')
    text=text.replace(old,'      write_tab(native[is.finite(qvalue) & qvalue <= .05 & (!paper_aligned | end-start+1 > 50)], "dmr.tsv")')
    Path(target).write_text(text)

def build_epykit_runner(original,target):
    text=Path(original).read_text()
    text=replace_once(text,'"permutation_max_t", "permutation_region"}',
        '"permutation_max_t", "permutation_region", "paper_aligned"}')
    text=replace_once(text,'    if profile == "historical_site_smooth":',
        '    if profile == "paper_aligned":\n        dmc_kwargs.update(smoothing=True, smoothing_span_bp=500)\n        dmr_kwargs.update(alpha=1e-5, min_abs_meth_diff=0., minlen_bp=51, min_cpgs=4, dis_merge_bp=100, min_mean_qvalue=None)\n    if profile == "historical_site_smooth":')
    text=replace_once(text,'    strand_collapsed = dataset_manifest.get("kind") == "real" and dataset_manifest.get("study") == "GSE64177"',
        '    strand_collapsed = bool(dataset_manifest.get("strand_collapsed")) or (dataset_manifest.get("kind") == "real" and dataset_manifest.get("study") == "GSE64177")')
    Path(target).write_text(text)
