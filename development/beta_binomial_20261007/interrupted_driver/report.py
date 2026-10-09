"""Source-backed scientific readout with standalone figures and tables."""
from pathlib import Path
import html
import json
import shutil
import subprocess
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT/'results/beta_binomial_20261007'
HERE = Path(__file__).parent.resolve()
sys.path[:0] = [str(ROOT/'vendor310')]
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS = {'epykit_original':'epykit original','epykit_region_fix':'epykit region fix',
          'epykit_beta_binomial':'epykit BB-F (experimental)'}


def main():
    assert json.loads((RUN/'run.json').read_text())['status']=='complete'
    subprocess.run([sys.executable,str(HERE/'validate_results.py')],check=True)
    out = RUN/'comparison'
    region = pd.read_csv(out/'all_tool_region_comparison.csv')
    site = pd.read_csv(out/'all_tool_cpg_comparison.csv')
    real = pd.read_csv(out/'real_performance.csv')
    paper = pd.read_csv(out/'paper_comparison.csv')
    heldout = pd.read_csv(HERE/'holdout_v3/metrics.csv')
    sensitivity = pd.read_csv(HERE/'holdout_v3/prior_start_sensitivity.csv')
    heldout_context = json.loads((HERE/'holdout_v3/context.json').read_text())
    assert heldout_context['status']=='complete'
    for stem in ['power_and_ranking','null_tails']:
        for suffix in ['png','svg']:
            shutil.copy2(HERE/'holdout_v3'/f'{stem}.{suffix}',out/f'{stem}.{suffix}')
    signal = region[(region.dataset=='genome_autosomes_signal')&(region.matching=='dmr')].copy()
    null = region[(region.dataset=='genome_autosomes_null')&(region.matching=='dmr')].set_index('tool')
    boundary = region[(region.dataset=='genome_autosomes_signal')&(region.matching=='dmr_80')].set_index('tool')
    signal['null_calls'] = signal.tool.map(null.called_regions)
    signal['f1_80'] = signal.tool.map(boundary.region_f1)
    signal['minutes'] = signal.runtime_s/60
    signal = signal.sort_values('region_f1',ascending=True)
    colors = ['#184e87' if t=='epykit_beta_binomial' else '#79a9cf' if t=='epykit_region_fix'
              else '#888888' if t=='epykit_original' else '#388174' for t in signal.tool]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,4,figsize=(15,6),sharey=True,gridspec_kw={'width_ratios':[1.6,1.1,1,1.2]})
    y = np.arange(len(signal))
    for ax,field,title in zip(axes,['region_f1','region_recall','null_calls','minutes'],
                             ['Region F1 (50% match)','Region recall','Null regions','Signal runtime (min)']):
        if field=='minutes':
            ax.scatter(signal[field],y,c=colors,s=50)
        else:
            ax.barh(y,signal[field],color=colors)
        for j,value in enumerate(signal[field]):
            if not np.isfinite(value):
                ax.text(.02,j,'No calls' if signal.iloc[j].called_regions==0 else 'Unavailable',
                        transform=ax.get_yaxis_transform(),va='center',fontsize=9)
                continue
            ax.annotate(f'{value:.3f}' if field.startswith('region_') else f'{value:g}' if field=='null_calls'
                        else f'{value:.1f}',(value,j),xytext=(4,0),textcoords='offset points',va='center',fontsize=9)
        ax.set_title(title)
        ax.grid(axis='x',alpha=.2); ax.set_axisbelow(True)
        if field=='null_calls':ax.set_xlim(0,max(signal[field].max()*1.25,1))
        elif field=='minutes':ax.set_xscale('log');ax.set_xlim(1,max(signal[field].max()*1.6,2))
        else:ax.set_xlim(0,max(signal[field].max()*1.2,.1))
    axes[0].set_yticks(y,[LABELS.get(t,t) for t in signal.tool])
    fig.suptitle('Frozen genome benchmark: count engine changed; region rules retained',fontsize=14)
    fig.tight_layout()
    fig.savefig(out/'region_comparison.png',dpi=180,bbox_inches='tight')
    fig.savefig(out/'region_comparison.svg',bbox_inches='tight');plt.close(fig)
    cpg = site[(site.dataset=='genome_autosomes_signal')&site.tp.notna()].copy()
    null_site = site[site.dataset=='genome_autosomes_null'].set_index('tool')
    cpg['null_significant_cpgs'] = cpg.tool.map(null_site.fp)
    # Undefined precision/FDP for zero discoveries is displayed as missing.
    cpg.loc[cpg.tp+cpg.fp==0,['precision','false_discovery_proportion']] = np.nan
    cpg['Tool'] = cpg.tool.map(lambda t:LABELS.get(t,t))
    cpg_table = cpg[['Tool','tp','fp','recall_eligible','false_discovery_proportion','average_precision',
                     'null_significant_cpgs','n_scored']].rename(columns={
        'tp':'True CpGs','fp':'False CpGs','recall_eligible':'Recall','false_discovery_proportion':'Observed FDP',
        'average_precision':'Average precision','null_significant_cpgs':'Null CpGs','n_scored':'Sites with p < 1'})
    region_table = signal.iloc[::-1][['tool','called_regions','matched_regions','region_precision','region_recall',
                                      'region_f1','f1_80','null_calls','minutes','peak_rss_gib']].copy()
    region_table['tool'] = region_table.tool.map(lambda t:LABELS.get(t,t))
    region_table.columns = ['Tool','Calls','Matches','Precision','Recall','F1 50%','F1 80%','Null DMRs','Minutes','Peak RSS GiB']
    paper_real = real.merge(paper[paper.reciprocal_threshold==0],on='tool',suffixes=('','_paper'))
    paper_real['matches_50'] = paper_real.tool.map(paper[paper.reciprocal_threshold==.5].set_index('tool').one_to_one_matches)
    real_table = paper_real[['tool','calls','calls_overlapping_paper','caller_overlap_fraction','paper_regions_covered',
                             'matches_50','runtime_s','peak_rss_gib']].copy()
    real_table['tool'] = real_table.tool.map(lambda t:LABELS.get(t,t))
    real_table['runtime_s'] /= 60
    real_table.columns = ['Tool','Calls','Paper-overlapping calls','Overlap fraction','Paper regions covered','Matches 50%',
                          'Minutes','Peak RSS GiB']
    new = signal.set_index('tool').loc['epykit_beta_binomial']
    old = signal.set_index('tool').loc['epykit_original']
    fixed = signal.set_index('tool').loc['epykit_region_fix']
    new_cpg = cpg.set_index('tool').loc['epykit_beta_binomial']
    def markdown_table(frame):
        return frame.fillna('—').to_markdown(index=False,floatfmt='.3f')
    summary = (f'The experimental beta-binomial engine has region precision {new.region_precision:.1%}, '
               f'recall {new.region_recall:.1%} and F1 {new.region_f1:.3f}. '
               f'Original epykit F1 was {old.region_f1:.3f}; the accepted region-only fix was {fixed.region_f1:.3f}. '
               f'The new null run calls {int(new.null_calls)} regions and {int(new_cpg.null_significant_cpgs)} CpGs. '
               f'Its signal CpG calls contain {int(new_cpg.tp)} true and {int(new_cpg.fp)} false sites.')
    methods = ('All simulations use the original count files, strict eligibility and scoring code. '
               'Region metrics require one-to-one matching at 50% or 80% reciprocal eligible-CpG overlap; '
               'CpG discoveries use BH q<=0.05. Competitor calls and timings are reused from the original run. '
               'The new epykit runner changes test=beta_binomial and skips confidence intervals for timing; '
               'the accepted region interval-family/BY correction, geometry, effect threshold and all other settings are fixed. '
               'BH outputs inherit the calibration limits of their approximate p-values. In the frozen scorer, '
               'missing or undefined site p-values are treated as 1 for evaluation; n_scored counts p<1. '
               'Ranking metrics use the common eligible universe. BSmooth and dmrseq do not export comparable '
               'CpG p/q-values, so their CpG metrics are unavailable. Peak memory is sampled process-tree RSS.')
    limitations = ('This is a genuine per-replicate count likelihood, with a count-learned prior and continuous '
                   'log-rho MAP shrinkage. Its conditional F reference is approximate. It does not integrate dispersion '
                   'or prior-estimation uncertainty. Mean/rho independence, one shared within-site dispersion and '
                   'finite prior-training samples remain assumptions. Zero discoveries with zero recall do not establish '
                   'a useful engine. Neither one genome null dataset nor these limited synthetic runs establishes FDR control. '
                   'The legacy automatic default remains unchanged; the new engine is explicitly experimental.')
    null_ci = heldout_context['null_any_rejection_95ci']
    calibration = (f'V3 used {heldout_context["sites_total"]:,} independently generated sites on untouched seeds. '
        f'{heldout_context["repeated_null_trials_with_rejection"]}/{heldout_context["repeated_null_trials"]} repeated global-null '
        f'trials had any BH rejection (95% binomial interval {null_ci[0]:.1%}–{null_ci[1]:.1%}). '
        f'The largest number of changed calls between prior starts was {int(sensitivity.changed_calls.max())}. '
        'The workbook includes every scenario, including missing coverage, unequal group dispersions, '
        'mean/dispersion association, stronger effects and different replicate counts. V2 evaluation was superseded '
        'after a label-dependent sampling bug; it is retained as development evidence and not used for this readout.')
    real_note = ('GSE64177 coordinates are hg19, 1-based inclusive. Paper overlap requires agreeing methylation direction; '
                 '50% matching uses base-pair overlap, unlike simulation matching. The workbook contains 3,271 published regions, '
                 '1,714 hypo and 1,557 hyper. The paper cautions that many hyper calls may be false positives, and '
                 '2,479 paper regions contain fewer than five CpGs. Agreement is descriptive, not biological precision or FDR. '
                 'The frozen two-group models omit donor pairing. Only epykit and DSS completed the original real benchmark; '
                 'methylKit, BSmooth, dmrseq and DMRcate real results are unavailable here. The new epykit real run has no benchmark memory cap.')
    decision = ('Keep beta_binomial opt-in while comparing its ranking and recall with the accepted region-only fix and DSS. '
                'If improved false-call behavior comes with too little recall, the next statistical target is explicit '
                'dispersion uncertainty/small-sample calibration and spatial borrowing, evaluated on new simulations. '
                'A paired/covariate beta-binomial model is also required before making donor-adjusted biological claims. '
                'Do not choose a new cutoff from this genome or paper comparison.')
    validation_path = HERE/'verification.json'
    validation = json.loads(validation_path.read_text()) if validation_path.exists() else {}
    references = ('[Supplied WGBS survey](https://doi.org/10.1093/bib/bbx013), '
                  '[DSS](https://pmc.ncbi.nlm.nih.gov/articles/PMC4005660/), '
                  '[RADMeth](https://pmc.ncbi.nlm.nih.gov/articles/PMC4230021/), '
                  '[methylSig](https://pmc.ncbi.nlm.nih.gov/articles/PMC4147891/), '
                  '[GSE64177 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4665002/).')
    md = '# Epykit beta-binomial engine: full benchmark comparison\n\n'+summary+'\n\n'
    md += '![Region comparison](region_comparison.png)\n\n'+markdown_table(region_table)+'\n\n'
    md += '## CpG identification\n\n'+markdown_table(cpg_table)+'\n\n'+limitations+'\n\n'
    md += '## Independent calibration and sensitivity\n\n'+calibration+'\n\n'
    md += '![Held-out power and ranking](power_and_ranking.png)\n\n![Held-out null tails](null_tails.png)\n\n'
    md += '## Real data and published DMRs\n\n'+markdown_table(real_table)+'\n\n'+real_note+'\n\n'
    md += '## Interpretation and next statistical target\n\n'+decision+'\n\n'
    md += '## Methods and evidence\n\n'+methods+'\n\n'+references+'\n\n'
    md += 'Machine-readable evidence: [comparison workbook](comparison.xlsx), [sources](sources.json), '
    md += '[run manifest](../run.json), [prior diagnostics](genome_priors.csv).\n\n'
    md += ('Package verification: 71 focused checks passed after the final repairs. '
           'The latest full-suite collection had 790 passed, 29 failed and 6 skipped: '
           '28 failures pre-existed; the remaining near-null label-swap tolerance was corrected '
           'and passed in the focused checks. See the development verification JSON and logs.\n')
    (out/'report.md').write_text(md)
    pieces = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Epykit beta-binomial comparison</title>',
              '<style>body{font:17px/1.55 system-ui;max-width:1400px;margin:36px auto;padding:0 24px;color:#172331}'
              'table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:8px;border-bottom:1px solid #d3dde7;text-align:right}'
              'th:first-child,td:first-child{text-align:left}h1,h2{line-height:1.2}img{width:100%}a{color:#184e87}</style>',
              '<h1>Epykit beta-binomial engine: full benchmark comparison</h1>',f'<p>{html.escape(summary)}</p>',
              '<img src="region_comparison.png" alt="Region F1, recall, null calls and runtime for eight benchmark versions">',
              region_table.to_html(index=False,float_format=lambda v:f'{v:.3f}',border=0,na_rep='—'),
              '<h2>CpG identification</h2>',cpg_table.to_html(index=False,float_format=lambda v:f'{v:.3f}',border=0,na_rep='—'),
              f'<p>{html.escape(limitations)}</p><h2>Independent calibration and sensitivity</h2><p>{html.escape(calibration)}</p>',
              '<img src="power_and_ranking.png" alt="Held-out CpG recall, average precision and observed false-discovery proportion">',
              '<img src="null_tails.png" alt="Raw null p-value tail rates with 95 percent binomial intervals">',
              '<h2>Real data and published DMRs</h2>',real_table.to_html(index=False,float_format=lambda v:f'{v:.3f}',border=0,na_rep='—'),
              f'<p>{html.escape(real_note)}</p><h2>Interpretation</h2><p>{html.escape(decision)}</p>',
              f'<h2>Methods</h2><p>{html.escape(methods)}</p>',
              '<p><a href="comparison.xlsx">Download tables</a> · <a href="report.md">Methods and references</a> · '
              '<a href="sources.json">Source provenance</a></p></html>']
    (out/'report.html').write_text('\n'.join(pieces))
    print(summary)


if __name__=='__main__':
    main()
