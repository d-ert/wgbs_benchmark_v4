"""Freeze a self-consistent comparison of completed GSE263850 callers."""
from pathlib import Path
import sys,json,datetime,hashlib
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'vendor310'),str(ROOT/'src')]
import polars as pl
import pandas as pd
import networkx as nx
base=ROOT/'results/GSE263850_paper_20261009'
manifest=json.load(open(base/'real/run_manifest.json'));state=json.load(open(base/'run.json'))
perf=pd.read_csv(base/'comparison/real_performance.csv');metrics=pd.read_csv(base/'comparison/paper_comparison.csv')
completed=[t for t in perf.tool if manifest[t]['status']=='ok' and manifest[t]['monitoring']['exit_code']==0]
perf=perf[perf.tool.isin(completed)];metrics=metrics[metrics.tool.isin(completed)]
paper=pd.read_csv(base/'input/paper_dmrs.tsv',sep='\t');assert len(paper)==813
checks=[]
for tool in completed:
 calls=pd.read_csv(base/'real'/tool/'dmr.tsv',sep='\t');assert len(calls)==int(perf[perf.tool==tool].calls.iloc[0])
 g=nx.Graph()
 for i,c in calls.iterrows():
  direction=int(__import__('numpy').sign(c.mean_meth_diff)) if 'mean_meth_diff' in calls else int(c.direction)
  for j,p in paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)].iterrows():
   overlap=min(c.end,p.end)-max(c.start,p.start)+1
   if direction==(1 if p.direction=='hyper' else -1) and overlap/(c.end-c.start+1)>=.5 and overlap/(p.end-p.start+1)>=.5:g.add_edge(('c',int(i)),('p',int(j)))
 nmatch=len(nx.max_weight_matching(g,maxcardinality=True));expected=int(metrics[(metrics.tool==tool)&(metrics.reciprocal_threshold==.5)].one_to_one_matches.iloc[0]);assert nmatch==expected
 checks.append(dict(tool=tool,calls=len(calls),independent_50pct_matches=nmatch))
perf.to_csv(OUT/'gse263850_performance.csv',index=False);metrics.to_csv(OUT/'gse263850_paper.csv',index=False)
dml=pl.scan_csv(base/'real/epykit/dml.tsv',separator='\t').select(pl.len().alias('sites'),(pl.col('pvalue')<1e-5).sum().alias('raw_seed_sites'),(pl.col('qvalue')<=.05).sum().alias('bh_calls')).collect().to_dicts()[0]
caches=list((base/'real/epykit/store').rglob('.dmr_chain_merge*.parquet'));assert len(caches)==1
c=pl.read_parquet(caches[0]);assert len(c)==51;assert int((c['combined_qvalue']<=.05).sum())==51
summary=dict(captured_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),run_state={k:state.get(k) for k in ['status','stage','active_tool','updated_at']},completed_tools=completed,completion_records={t:manifest[t] for t in completed},independent_match_checks=checks,epykit=dict(dml,pre_q_candidates=len(c),after_q_candidates=int((c['combined_qvalue']<=.05).sum())),settings=state['settings'],note='Frozen partial real-study evidence. Remaining tools were not complete at capture; no pending output treated as zero.')
(OUT/'gse263850_snapshot.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'completed':completed,'matches':checks,'epykit':summary['epykit']}))
