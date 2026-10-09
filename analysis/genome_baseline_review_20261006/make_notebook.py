from pathlib import Path
import contextlib
import io
import nbformat
from nbconvert import HTMLExporter

OUT=Path(__file__).resolve().parent
nb=nbformat.v4.new_notebook()
md=nbformat.v4.new_markdown_cell
code=nbformat.v4.new_code_cell
nb.cells=[
md('''# Genome baseline review — 6 October 2026

## tl;dr
The overall run failed on real-data methylKit (48 GiB memory limit). All twelve simulation runs completed; only epykit and DSS completed on real data.

On the signal simulation, dmrseq has the strongest 50% reciprocal-CpG region F1 (0.559), with no regions called in the null run. Epykit recovers almost the same number of true regions (1,202 versus 1,213), at 6.18 versus 485.09 minutes, but has 33.2% region precision and 1,869 null regions. Calibration of region inference is the main epykit development priority.

## Context & Methods
Decision: prioritize epykit improvements from the fixed v4 benchmark.

Source: `/scratch/wgbs_benchmark_v4/results/genome_baseline/`. Runs used 25,923,286 calibrated autosomal CpGs, five samples per group, nominal 20× depth, and one seed (20261003). Signal truth: 2,893 regions / 54,619 CpGs, with 34,889 positive CpGs in the 16,070,989-site common coverage mask. Null common coverage: 16,070,397 sites.

Region precision/recall require a one-to-one match with at least 50% reciprocal overlap in eligible CpGs; a stricter 80% version is also retained. Unmatched calls include boundary/localization failures as well as biologically null calls. `dmr_null_region_fdp` counts calls touching no positive eligible CpG; it is not the same as one minus strict region precision. Direction accuracy applies only to matched regions.

### Key Assumptions
This is one development simulation pair with no across-seed uncertainty estimate. The simulator's spatial-correlation option was disabled. At 20×, thinning can make signal/null realised coverage differ. Real GSE64177 uses six paired donors but the inherited runners omit pairing; concordance has no truth-based accuracy interpretation. Competitors use different statistical units and settings, including 1 kb methylKit windows and heuristic BSmooth/DSS regions. All comparisons are for these recorded configurations.

## Data'''),
code('''from pathlib import Path
import pandas as pd
ROOT=Path('/scratch/wgbs_benchmark_v4')
OUT=ROOT/'analysis/genome_baseline_review_20261006'
regions=pd.read_csv(OUT/'region_comparison.csv').set_index('tool')
dml=pd.read_csv(OUT/'dml_comparison.csv').set_index('tool')
real=pd.read_csv(OUT/'real_status.csv')
thresholds=pd.read_csv(OUT/'epykit_exploratory_threshold_tradeoff.csv')
assert len(regions)==6 and len(real)==6
print('Six tools; twelve completed simulation runs; two completed real runs.')'''),
md('''## Results
### Region detection on the signal simulation
Precision, recall and F1 use the 50% matching rule. Null discoveries are all false in this one null dataset.'''),
code('''table=regions[['dmr_called_regions','dmr_matched_regions','dmr_region_precision','dmr_region_recall','dmr_region_f1','null_dmr_calls']].rename(columns={'dmr_called_regions':'Calls','dmr_matched_regions':'Matched','dmr_region_precision':'Precision','dmr_region_recall':'Recall','dmr_region_f1':'F1','null_dmr_calls':'Null_calls'})
print(table.round(3).to_string())'''),
md('''Dmrseq has the best overall balance at this matching threshold. Epykit's 2,137 signal calls with no positive eligible CpG show that its precision problem is not explained just by boundary mismatch. DSS has 128 such calls, and dmrseq 56. DMRcate's one signal call overlaps a real signal but fails the reciprocal matching rule; its zero null calls accompany almost no power here.

### Runtime, memory, and stricter boundary matching'''),
code('''print(regions[['runtime_minutes','peak_rss_gib','dmr_80_region_f1']].rename(columns={'runtime_minutes':'Minutes','peak_rss_gib':'Peak_GiB','dmr_80_region_f1':'F1_at_80pct'}).round(3).to_string())'''),
md('''These are measured caller times, excluding scoring and preparation. Memory is peak sampled process-tree RSS and can double-count shared pages. At 80% reciprocal matching, epykit's F1 is 0.211 versus dmrseq's 0.188: boundary sensitivity changes the ranking. At 50%, matched-region direction accuracy is 100% for the five callers with matches.

### Single-CpG testing at q ≤ 0.05'''),
code('''cols=['dml_tp','dml_fp','dml_precision','dml_recall_eligible','dml_false_discovery_proportion','null_dml_calls']
print(dml[cols].dropna(subset=['dml_tp']).rename(columns={'dml_tp':'TP','dml_fp':'FP','dml_precision':'Precision','dml_recall_eligible':'Recall','dml_false_discovery_proportion':'Observed_FDP','null_dml_calls':'Null_calls'}).round(4).to_string())'''),
md('''Epykit recovers only 75 of 34,889 eligible positive CpGs at the site-level q threshold, while returning 323 false-positive sites. Its region recall comes from raw-p-value region seeding and aggregation, not from strong BH-significant CpG detection. DSS has much higher CpG recall but high site-level FDP under this exact CpG truth definition; smoothing near region boundaries may contribute, and this analysis does not establish the cause. BSmooth and dmrseq do not emit comparable CpG q-value calls in these runners. The observed FDPs are not estimates of across-replicate FDR.

### Exploratory filtering of existing epykit regions'''),
code('''print(thresholds[['q_cutoff','signal_calls','matched_regions','precision','recall','f1','null_calls']].round(5).to_string(index=False))'''),
md('''Changing the retained combined-q cutoff from 0.05 to 0.0001 leaves 1,137 matched regions (94.6% of the original matched count), increases precision to 74.6%, and leaves eight null calls. This is post-hoc exploration of the same data, not a validated new default or a claim of calibrated q-values. The independent compact overlap calculation reproduces the stored default epykit matched count and no-positive-CpG count before applying stricter thresholds.

### Partial real-data results'''),
code('''print(real[['tool','status','dmr_count','runtime_minutes','peak_rss_gib']].round(3).to_string(index=False))
print('\\nRegion agreement:')
print(pd.read_csv(OUT/'real_epykit_DSS/dmr_concordance.tsv',sep='\\t')[['reciprocal_50pct_matches','fraction_a_matched','fraction_b_matched','matched_direction_concordance']].round(3).to_string(index=False))
print('\\nGene agreement:')
print(pd.read_csv(OUT/'real_epykit_DSS/gene_concordance.tsv',sep='\\t')[['genes_a','genes_b','shared_genes','gene_jaccard']].round(3).to_string(index=False))'''),
md('''Epykit calls 431 regions (133 hyper / 298 hypo); DSS calls 253 (19 hyper / 234 hypo). There are 65 greedy one-to-one matches at ≥50% reciprocal base-pair overlap (15.1% of epykit, 25.7% of DSS), all with consistent direction. Gene overlap is 74 shared genes out of 320 epykit / 201 DSS genes. These are descriptive results only. MethylKit was terminated at 48.24 GiB after 222.5 minutes. Its partial files are excluded; the empty DMR file is not evidence of zero discoveries. BSmooth, dmrseq and DMRcate never ran on the real input.

## Takeaways
1. Preserve this run as an incomplete baseline; complete the missing real-data runs with memory-aware methylKit execution and a safe resume path. The current epykit-only wrapper requires a completed baseline.
2. Prioritize epykit region inference: selection from raw p-values followed by Stouffer aggregation and BH on selected regions needs empirical evaluation of the complete discovery pipeline. Investigate permutation or another selection-aware calibration method. The current code labels the asymptotic region statistic as a ranking signal.
3. Diagnose CpG calibration and power separately by methylation level, depth and dispersion, using additional seeds. Strong region recovery does not establish calibrated or powerful site-level inference.
4. Retain the speed/memory advantage while using held-out signal/null datasets to select and validate improved settings. The stricter threshold result is promising evidence for prioritization, not independent confirmation.

### Reproduction and checks
`analyze.py` produces the summary CSVs and checks the 12 raw DMR counts against score files, independently recomputes precision/recall/F1 arithmetic, and checks snapshotted harness files against the run's recorded hashes. `region_thresholds.py` computes threshold tradeoffs using saved eligibility and truth. Original and both simulation epykit source snapshots agree with the initial recorded source hashes.

Real annotations were generated with:
`Rscript analysis/annotate_real_dmrs.R results/genome_baseline/real analysis/genome_baseline_review_20261006/real_epykit_DSS epykit,DSS`

The full benchmark was not rerun or modified for this analysis.''')]
# Execute the simple, print-only companion cells in order without opening a kernel/server.
namespace={}
for count,cell in enumerate([c for c in nb.cells if c.cell_type=='code'],1):
    captured=io.StringIO()
    with contextlib.redirect_stdout(captured),contextlib.redirect_stderr(captured):
        exec(compile(cell.source,'genome_baseline_review.ipynb','exec'),namespace)
    cell.execution_count=count
    cell.outputs=[nbformat.v4.new_output('stream',name='stdout',text=captured.getvalue())]
nb.metadata.kernelspec={'display_name':'Python 3','language':'python','name':'python3'}
nb.metadata.execution_note='Code cells executed sequentially in a shared Python namespace; print-only outputs captured.'
nbformat.validate(nb)
nbformat.write(nb,OUT/'genome_baseline_review.ipynb')
html,_=HTMLExporter().from_notebook_node(nb)
(OUT/'genome_baseline_review.html').write_text(html)
for c in nb.cells:
    if c.cell_type=='code':print(''.join(o.get('text','') for o in c.outputs))
print('Notebook saved and all code cells executed successfully.')
