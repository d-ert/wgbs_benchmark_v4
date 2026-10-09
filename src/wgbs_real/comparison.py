import json
from pathlib import Path
import subprocess

import pandas as pd
from .regions import paper_metrics

def compare(run,annotate=True):
    run=Path(run);record=json.loads((run/'run.json').read_text());dataset=Path(record['dataset'])
    paper=pd.read_csv(dataset/'paper_dmrs.tsv',sep='\t')
    statuses=json.loads((run/'real/run_manifest.json').read_text())
    out=run/'comparison';out.mkdir(exist_ok=True)
    rows=[];pairs=[];performance=[];ok=[]
    for tool,entry in statuses.items():
        if entry['status']!='ok':continue
        calls=pd.read_csv(run/'real'/tool/'dmr.tsv',sep='\t')
        if len(calls) and not ((calls.start>=1)&(calls.end>=calls.start)).all():raise ValueError('Invalid output coordinates')
        metrics,overlaps=paper_metrics(calls,paper)
        rows.extend(dict(tool=tool,**m) for m in metrics)
        overlaps['tool']=tool;pairs.append(overlaps)
        monitor=entry['monitoring'];ok.append(tool)
        performance.append(dict(tool=tool,calls=len(calls),minutes=monitor['wall_seconds']/60,
            peak_rss_gib=monitor['peak_process_tree_rss_bytes']/2**30))
    if not ok:return
    metrics=pd.DataFrame(rows);perf=pd.DataFrame(performance)
    metrics.to_csv(out/'paper_comparison.csv',index=False)
    perf.to_csv(out/'real_performance.csv',index=False)
    overlaps=pd.concat(pairs,ignore_index=True);overlaps.to_csv(out/'paper_overlap_pairs.csv',index=False)
    with pd.ExcelWriter(out/'comparison.xlsx',engine='openpyxl') as writer:
        metrics.to_excel(writer,sheet_name='Paper comparison',index=False)
        perf.to_excel(writer,sheet_name='Performance',index=False)
        paper.to_excel(writer,sheet_name='Published DMRs',index=False)
        overlaps.to_excel(writer,sheet_name='Overlap pairs',index=False)
    table=perf.merge(metrics[metrics.reciprocal_threshold==.5],on=['tool','calls'])
    text='# '+record['study']+' benchmark\n\n'+table[['tool','calls','one_to_one_matches','one_to_one_precision','one_to_one_recall','minutes','peak_rss_gib']].to_markdown(index=False,floatfmt='.3f')
    text+='\n\nPrecision/recall treat the published list as the reference; matching is one-to-one, same direction, 50% reciprocal base-pair overlap. Any and 80% overlap results are in the workbook. Paper agreement is not independent biological validation.\n'
    text+='\nRunner mappings and remaining differences: '+json.dumps(record['settings'],indent=2)+'\n'
    (out/'README.md').write_text(text)
    if annotate and len(ok)>=2:
        subprocess.run(['Rscript',str(run/'source_snapshot/annotate_hg38.R'),str(run/'real'),
            str(out/'concordance'),','.join(ok)],check=True)
    print(table[['tool','calls','one_to_one_matches','one_to_one_precision','one_to_one_recall']].to_string(index=False),flush=True)
