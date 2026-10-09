"""Read-only review of completed genome_baseline outputs; run from any directory."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
RUN=ROOT/'results/genome_baseline'
metrics=pd.read_csv(RUN/'simulated/simulation_metrics.csv')
assert len(metrics)==12 and not metrics.duplicated(['dataset','tool']).any()
assert metrics.status.eq('ok').all()
checks=[]
for row in metrics.itertuples(index=False):
    folder=RUN/'simulated'/row.dataset/row.tool
    score=json.loads((folder/'score.json').read_text())
    raw=pd.read_csv(folder/'dmr.tsv',sep='\t')
    assert len(raw)==score['dmr']['called_regions']==row.dmr_called_regions
    matched=score['dmr']['matched_regions']; called=len(raw); truth=score['dmr']['true_regions']
    if called: assert np.isclose(matched/called,row.dmr_region_precision)
    if truth: assert np.isclose(matched/truth,row.dmr_region_recall)
    if called+truth: assert np.isclose(2*matched/(called+truth),row.dmr_region_f1)
    checks.append(dict(dataset=row.dataset,tool=row.tool,raw_dmr_rows=len(raw),score_matches=True))
metrics['runtime_minutes']=metrics.wall_seconds/60
metrics['peak_rss_gib']=metrics.peak_rss_bytes/2**30
metrics.to_csv(OUT/'all_simulation_metrics.csv',index=False)
signal=metrics[metrics.dataset.str.endswith('_signal')].set_index('tool')
null=metrics[metrics.dataset.str.endswith('_null')].set_index('tool')
region=signal[['dmr_called_regions','dmr_matched_regions','dmr_region_precision','dmr_region_recall','dmr_region_f1','dmr_80_region_precision','dmr_80_region_recall','dmr_80_region_f1','dmr_calls_without_positive_eligible_cpg','dmr_null_region_fdp','dmr_cpg_precision','dmr_cpg_recall','dmr_matched_direction_accuracy','runtime_minutes','peak_rss_gib']].copy()
region['null_dmr_calls']=null.dmr_called_regions
region.to_csv(OUT/'region_comparison.csv')
dml=signal[[c for c in signal if c.startswith('dml_')]].copy()
dml['null_dml_calls']=null.dml_tp+null.dml_fp
dml.to_csv(OUT/'dml_comparison.csv')
records=json.loads((RUN/'real/run_manifest.json').read_text())
real=[]
for tool in ['epykit','DSS','methylKit','BSmooth','dmrseq','DMRcate']:
    record=records.get(tool,{})
    status=record.get('status','not_run')
    m=record.get('monitoring',{})
    row=dict(tool=tool,status=status,termination=m.get('termination_reason'),
        runtime_minutes=m.get('wall_seconds',np.nan)/60,peak_rss_gib=m.get('peak_process_tree_rss_bytes',np.nan)/2**30)
    if status=='ok':
        calls=pd.read_csv(RUN/'real'/tool/'dmr.tsv',sep='\t')
        direction=np.sign(calls.mean_meth_diff) if 'mean_meth_diff' in calls else calls.direction
        row.update(dmr_count=len(calls),hypermethylated=int((direction>0).sum()),hypomethylated=int((direction<0).sum()),
                   median_length_bp=float((calls.end-calls.start+1).median()))
    real.append(row)
pd.DataFrame(real).to_csv(OUT/'real_status.csv',index=False)
epy=[]
for dataset in ['genome_autosomes_signal','genome_autosomes_null']:
    calls=pd.read_csv(RUN/'simulated'/dataset/'epykit/dmr.tsv',sep='\t')
    for cutoff in [.05,.01,.001,.0001,.00001]:
        sub=calls[calls.combined_qvalue<=cutoff]
        epy.append(dict(dataset=dataset,q_cutoff=cutoff,calls=len(sub),
            mean_effect_below_10pct=int((sub.mean_meth_diff.abs()<.1).sum())))
pd.DataFrame(epy).to_csv(OUT/'epykit_region_threshold_diagnostics.csv',index=False)
# Verify saved benchmark source agrees with the run's recorded harness.
run_meta=json.loads((RUN/'run.json').read_text())
for rel,expected in run_meta['contract']['harness'].items():
    p=RUN/'source_snapshot/benchmark'/rel
    checks.append(dict(source_file=rel,snapshot_matches=p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==expected))
# Real concordance is written separately by the existing R annotation script.
(OUT/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print('REGIONS\n',region[['dmr_called_regions','dmr_matched_regions','dmr_region_precision','dmr_region_recall','dmr_region_f1','null_dmr_calls','runtime_minutes','peak_rss_gib']].round(4).to_string())
print('\nREAL\n',pd.DataFrame(real).to_string(index=False))
print('\nEPYKIT THRESHOLDS\n',pd.DataFrame(epy).to_string(index=False))
print('\nSnapshot mismatches:',[c for c in checks if c.get('snapshot_matches') is False])
