"""Write the completed benchmark review and scientific comparison figure."""
from pathlib import Path
import base64
import hashlib
import html
import json
import subprocess
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT / 'results/default_search_by_20261007'
OUT = RUN / 'comparison'
sys.path[:0] = [str(ROOT / 'vendor310'), str(ROOT / 'src'), str(ROOT)]
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def percentage(value):
    return f'{100*value:.1f}%' if pd.notna(value) else 'N/A'


def main():
    run = json.loads((RUN / 'run.json').read_text())
    assert run['status'] == 'complete'
    checks = json.loads((RUN / 'regression_checks.json').read_text())
    regions = pd.read_csv(OUT / 'epykit_region_comparison.csv')
    all_tools = pd.read_csv(OUT / 'all_tool_simulation_metrics.csv')
    signal = all_tools[all_tools.dataset == 'genome_autosomes_signal'].copy()
    null = all_tools[all_tools.dataset == 'genome_autosomes_null'].copy()
    label = lambda r: ('epykit original' if r.source_run == 'genome_baseline' else 'epykit revised') if r.tool == 'epykit' else r.tool
    signal['label'] = signal.apply(label, axis=1); null['label'] = null.apply(label, axis=1)
    signal = signal.merge(null[['label','dmr_called_regions']], on='label', suffixes=('','_null'))
    signal = signal.sort_values('dmr_region_f1', ascending=False)
    table = pd.DataFrame({
        'Tool': signal.label, 'Precision': signal.dmr_region_precision.map(percentage),
        'Recall': signal.dmr_region_recall.map(percentage), 'F1': signal.dmr_region_f1.map(lambda x:f'{x:.3f}'),
        'Matched true regions': signal.dmr_matched_regions.astype(int),
        'Signal calls': signal.dmr_called_regions.astype(int), 'Null calls': signal.dmr_called_regions_null.astype(int),
        'Runtime (min)': (signal.wall_seconds/60).map(lambda x:f'{x:.2f}'),
        'Peak RSS (GiB)': (signal.peak_rss_bytes/2**30).map(lambda x:f'{x:.2f}'),
        'Measurement': signal.reused.map({True:'recorded original',False:'fresh rerun'}),
    })
    table.to_csv(OUT / 'all_tool_summary.csv', index=False)
    old = regions[(regions.dataset == 'genome_autosomes_signal') & (regions.version == 'original') & (regions.matching == 'dmr')].iloc[0]
    new = regions[(regions.dataset == 'genome_autosomes_signal') & (regions.version == 'revised') & (regions.matching == 'dmr')].iloc[0]
    null_new = int(regions[(regions.dataset == 'genome_autosomes_null') & (regions.version == 'revised') & (regions.matching == 'dmr')].called_regions.iloc[0])
    boundary_rows = []; site_rows = []
    for row in signal.itertuples():
        score = json.loads((ROOT / 'results' / row.source_run / 'simulated/genome_autosomes_signal' / row.tool / 'score.json').read_text())
        boundary_rows.append(dict(Tool=row.label, **score['dmr_80'], source_run=row.source_run))
        site_rows.append(dict(Tool=row.label, **(score.get('dml') or {}), source_run=row.source_run))
    boundary_metrics = pd.DataFrame(boundary_rows)
    site_metrics = pd.DataFrame(site_rows)
    boundary_metrics.to_csv(OUT/'all_tool_boundary80.csv',index=False)
    site_metrics.to_csv(OUT/'all_tool_cpg_metrics.csv',index=False)
    sensitivity = boundary_metrics[['Tool','called_regions','matched_regions','region_precision','region_recall','region_f1']].copy()
    for col in ['region_precision','region_recall']:
        sensitivity[col] = sensitivity[col].map(percentage)
    real = pd.read_csv(OUT / 'real_performance.csv')
    paper = pd.read_csv(OUT / 'paper_comparison.csv')
    directions = pd.read_csv(OUT / 'paper_comparison_by_direction.csv')
    concise_paper = paper[['tool','reciprocal_threshold','calls','calls_overlapping_paper','caller_overlap_fraction','paper_regions_covered','one_to_one_matches','discordant_pairs']].copy()
    concise_paper['caller_overlap_fraction'] = concise_paper.caller_overlap_fraction.map(percentage)
    cpg = pd.read_csv(OUT / 'epykit_cpg_comparison.csv')
    cpg_fields = ['dataset','version','tp','fp','fn','precision','recall_eligible','false_discovery_proportion']
    cpg_table = cpg[cpg_fields].copy()
    for column in ['precision','recall_eligible','false_discovery_proportion']:
        cpg_table[column] = cpg_table[column].map(percentage)
    site_table = site_metrics[['Tool','tp','fp','precision','recall_eligible','false_discovery_proportion']].copy()
    for column in ['tp','fp']:
        site_table[column] = site_table[column].map(lambda x:f'{int(x):,}' if pd.notna(x) else 'N/A')
    for column in ['precision','false_discovery_proportion']:
        site_table[column] = site_table[column].map(percentage)
    site_table['recall_eligible'] = site_table.recall_eligible.map(lambda x:f'{100*x:.3f}%' if pd.notna(x) else 'N/A')

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5,5.4), gridspec_kw={'width_ratios':[1.8,1]})
    ypos = np.arange(len(signal))
    ax.barh(ypos-.17, 100*signal.dmr_region_precision, height=.32, label='Precision',color='#145f73')
    ax.barh(ypos+.17, 100*signal.dmr_region_recall, height=.32, label='Recall',color='#b77a35')
    ax.set_yticks(ypos, signal.label); ax.invert_yaxis(); ax.set_xlim(0,100)
    ax.set_xlabel('Percent'); ax.set_title('Signal: ≥50% reciprocal eligible-CpG overlap'); ax.legend(loc='lower right')
    colors = ['#267b67' if x == 'epykit revised' else '#8898a0' for x in signal.label]
    values = signal.dmr_called_regions_null.to_numpy()
    bx.barh(ypos, values, color=colors)
    bx.set_yticks(ypos, signal.label); bx.invert_yaxis(); bx.set_xscale('symlog',linthresh=1); bx.set_xlim(0,max(values)*2.8)
    bx.set_xlabel('Null DMR count (linear below 1, log above)'); bx.set_title('Complete-null discoveries')
    for y, value in zip(ypos, values):
        bx.text(value*1.08 if value else .06,y,str(int(value)),va='center')
    fig.suptitle('Fresh default epykit region correction vs recorded original tools',fontsize=13)
    fig.text(.02,.01,'2,893 true signal regions; one development seed. Competitor runs are reused.',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,.94]); fig.savefig(OUT/'comparison.png',dpi=180); plt.close(fig)

    findings = (
        f'Epykit region precision changed from {percentage(old.region_precision)} to {percentage(new.region_precision)}, '
        f'recall from {percentage(old.region_recall)} to {percentage(new.region_recall)}, '
        f'and F1 from {old.region_f1:.3f} to {new.region_f1:.3f}. '
        f'The revised caller matched {int(new.matched_regions):,} true regions versus {int(old.matched_regions):,} originally '
        f'({100*new.matched_regions/old.matched_regions:.1f}% of the original match count). Null calls fell from 1,869 to {null_new:,}.'
    )
    sections = [
        ('Result', findings, None),
        ('What changed', 'The default chain-merge caller now uses Benjamini–Yekutieli (BY) adjustment over every contiguous interval allowed by its existing CpG count, length and gap constraints. Unreported intervals receive p=1. The original caller applied BH only after p-value-driven selection. The correction was derived from the search geometry; the significance threshold remains 0.05. This branch starts from the frozen original source and excludes yesterday’s permutation and window-screen changes.', None),
        ('All-tool simulation comparison', 'Both original full autosomal inputs were reused, with five samples per group and 2,893 true signal regions. The primary scorer requires ≥50% reciprocal overlap in eligible CpGs with maximum-cardinality one-to-one matching. Region definitions and statistical units differ between tools. Runtime includes caller execution; peak memory is sampled process-tree RSS. Competitor timings describe the recorded original runs and are not contemporaneous measurements.', table),
        ('Boundary sensitivity', 'The same scorer also evaluates ≥80% reciprocal eligible-CpG overlap. Retained calls keep their original boundaries; this patch changes significance, not segmentation. Matching precision penalizes localization errors and does not equal a calibrated biological false-discovery rate.', sensitivity),
        ('CpG testing', 'All three recomputed CpG TSVs are byte-for-byte identical to the original run. At site q ≤ 0.05 on the signal input, the 398 significant CpGs include 75 true positives and 323 false positives: 81.2% observed false discoveries and 0.22% recall among 34,889 eligible positive CpGs. The raw site test, dispersion model and site-level multiple testing are unchanged. These site-level metrics differ from the accuracy of CpGs covered by predicted regions: region aggregation can recover weak site evidence without any of those sites passing the individual-site q threshold. The CpG power and false-positive problems remain unresolved.', cpg_table),
        ('CpG comparison across tools', 'Site-level q ≤ 0.05 is evaluated against the same 34,889 eligible positive CpGs. DSS recovers substantially more sites than epykit but also has high observed false discoveries. DMRcate makes just one true-positive site call, so its zero observed false discoveries accompany almost no recovery. BSmooth and dmrseq do not expose comparable site-level calls in this runner; their missing values are N/A, not zero. Competitor results remain the recorded original measurements.', site_table),
        ('Real-data performance', 'All twelve GSE64177 samples were recomputed with the same coverage intersection, hg19 coordinates and original unpaired model. The donor pairing remains unmodeled. DSS is reused from the original run; methylKit failed and the other real competitors were not run in that baseline.', real),
        ('Published DMR agreement', 'The reference workbook has 3,271 hg19 regions. Counts below require matching methylation direction and use inclusive base-pair overlap. Any-overlap and reciprocal counts allow many-to-many relationships; one-to-one matches are shown separately. Paper agreement is descriptive, not biological precision or FDR. Published hypermethylated calls carry the paper’s false-positive caveat, and 2,479 published regions contain fewer than five CpGs, while this runner requires five.', concise_paper),
        ('Paper agreement by direction', 'Hypomethylated and hypermethylated regions are retained separately so the less reliable published hyper calls do not obscure the hypo comparison.', directions),
        ('Verification and limits', 'The original input and harness hash contract matched before and after the run. Caller settings match exactly, optional empirical_fdr remains false, source snapshots are frozen, and every surviving region has an original boundary. Focused checks: 43 passed. Full suite: 750 passed, 29 failed, 6 skipped; all 29 failures reproduced with the original package. BY accounts for search and overlapping interval dependence only if each fixed-span raw p-value is valid. The original CpG calibration and within-span Stouffer independence assumptions remain limitations. This is one already-inspected simulation seed; independent realistic seeds with spatial correlation are needed for general performance and calibrated region-FDR claims.', None),
    ]
    md = ['# Epykit default region correction — full benchmark review', '', findings, '', '![Simulation comparison](comparison.png)', '']
    html_sections = []
    for title, text, data in sections:
        md.extend([f'## {title}', '', text, ''])
        body = f'<section><h2>{html.escape(title)}</h2><p>{html.escape(text)}</p>'
        if data is not None:
            md.extend([data.to_markdown(index=False), ''])
            body += '<div class="table">'+data.to_html(index=False,escape=True,border=0)+'</div>'
        html_sections.append(body+'</section>')
    reference = 'BY method: https://doi.org/10.1214/aos/1013699998. GSE64177 paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC4665002/.'
    md += ['## Reproducibility', '', 'Source checkout: `/scratch/epykit-calibration-fresh_20261007`. Baseline: `results/genome_baseline`. Candidate: `results/default_search_by_20261007`. Run, analysis and report scripts: `development/default_search_by_20261007`. Raw tables, overlap pairs, source snapshots and hash checks accompany this report.', '', reference]
    (OUT/'report.md').write_text('\n'.join(md)+'\n')
    encoded = base64.b64encode((OUT/'comparison.png').read_bytes()).decode()
    page = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Epykit benchmark review</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#18323b;max-width:1200px;margin:40px auto;padding:0 24px;background:#fafcfc}h1{font-size:30px}h2{font-size:22px;margin-top:32px}section{margin:22px 0}.table{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px 10px;text-align:left;border-bottom:1px solid #d7e1e4}th{background:#eaf1f2}img{width:100%;height:auto}footer{font-size:13px;margin-top:36px}</style><h1>Epykit default region correction</h1><p>'+html.escape(findings)+'</p><img alt="Precision and recall on the signal simulation; null discovery counts for each tool" src="data:image/png;base64,'+encoded+'">'+''.join(html_sections)+'<footer>'+html.escape(reference)+'</footer></html>'
    (OUT/'report.html').write_text(page)
    with pd.ExcelWriter(OUT/'comparison.xlsx',engine='openpyxl') as workbook:
        for name,data in [('All tools',table),('All tools at 80pct',boundary_metrics),('All tools CpGs',site_metrics),('Regions',regions),('CpGs',cpg),('Real performance',real),('Paper overlap',paper),('Paper directions',directions)]:
            data.to_excel(workbook,sheet_name=name,index=False)
    provenance = dict(candidate=str(RUN),baseline=str(ROOT/'results/genome_baseline'),
                      source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd='/scratch/epykit-calibration-fresh_20261007',text=True).strip(),
                      script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.iterdir() if p.suffix in {'.py','.R'}},
                      paper_workbook_sha256=hashlib.sha256((ROOT/'data/real/paper_reported_dmrs.xlsx').read_bytes()).hexdigest(),
                      independent_paper_validation=(Path(__file__).parent/'paper_validation.log').read_text().strip(),
                      search_families=json.loads((RUN/'search_families.json').read_text()),
                      scope='Only epykit freshly rerun; original six simulation tools and real DSS reused',
                      regression_checks=checks)
    (OUT/'analysis_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(findings)


if __name__ == '__main__':
    main()
