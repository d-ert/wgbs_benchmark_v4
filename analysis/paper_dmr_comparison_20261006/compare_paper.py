"""Compare real epykit/DSS intervals with paper Table S2; no benchmark mutations."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import networkx as nx

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
XLSX=ROOT/'data/real/paper_reported_dmrs.xlsx'
RUN=ROOT/'results/genome_baseline/real'
raw=pd.read_excel(XLSX,header=2)
valid=raw.Chr.astype(str).str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)') & raw.Start.notna() & raw.End.notna()
paper=raw.loc[valid].copy().rename(columns={'Chr':'chrom','Start':'start','End':'end','Direction':'direction','# CpGs':'n_cpgs','Length (bp)':'length_bp','Hugo gene ID':'paper_gene','log2FCa':'gene_expression_log2FC','FDRb':'gene_expression_FDR'})
paper.insert(0,'paper_id',[f'paper_{i:04}' for i in range(1,len(paper)+1)])
paper.insert(1,'excel_row',paper.index+4)
paper=paper.reset_index(drop=True)
for col in ['start','end','n_cpgs','length_bp']:
 assert np.isfinite(paper[col]).all() and (paper[col]%1==0).all()
 paper[col]=paper[col].astype(int)
assert len(paper)==3271 and not paper.duplicated(['chrom','start','end']).any()
assert (paper.start>=1).all() and (paper.end>=paper.start).all()
assert (paper.end-paper.start+1==paper.length_bp).all()
assert set(paper.direction)=={'hypo','hyper'}
assert json.loads((ROOT/'data/real/manifest.json').read_text())['assembly']=='hg19'
summary=[]; direction_summary=[]; all_pairs=[]; one_to_one=[]; call_tables={}
rules={'any_overlap':0,'reciprocal_25pct':.25,'reciprocal_50pct':.5,'reciprocal_80pct':.8}
for tool in ['epykit','DSS']:
 calls=pd.read_csv(RUN/tool/'dmr.tsv',sep='\t')
 calls.insert(0,'call_id',[f'{tool}_{i:04}' for i in range(1,len(calls)+1)])
 assert not calls.duplicated(['chrom','start','end']).any()
 calls['length_bp']=calls.end-calls.start+1
 if tool=='epykit': calls['direction']=np.where(calls.mean_meth_diff>0,'hyper',np.where(calls.mean_meth_diff<0,'hypo','unknown'))
 else:calls['direction']=calls.direction.map({1:'hyper',-1:'hypo',0:'unknown'})
 pairs=[]
 for i,c in calls.iterrows():
  candidates=paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)]
  for j,p in candidates.iterrows():
   overlap=min(c.end,p.end)-max(c.start,p.start)+1
   pairs.append(dict(tool=tool,call_index=i,paper_index=j,call_id=c.call_id,paper_id=p.paper_id,
      chrom=c.chrom,call_start=int(c.start),call_end=int(c.end),paper_start=int(p.start),paper_end=int(p.end),
      call_direction=c.direction,paper_direction=p.direction,same_direction=c.direction==p.direction,
      overlap_bp=int(overlap),fraction_call=overlap/c.length_bp,fraction_paper=overlap/p.length_bp,
      interval_jaccard=overlap/(c.length_bp+p.length_bp-overlap),paper_gene=p.paper_gene))
 pairs=pd.DataFrame(pairs)
 all_pairs.append(pairs)
 for rule,cutoff in rules.items():
  chosen=pairs[(pairs.fraction_call>=cutoff)&(pairs.fraction_paper>=cutoff)]
  same=chosen[chosen.same_direction]
  graph=nx.Graph()
  graph.add_weighted_edges_from((('c',int(r.call_index)),('p',int(r.paper_index)),float(r.interval_jaccard)) for r in same.itertuples())
  matched=nx.max_weight_matching(graph,maxcardinality=True,weight='weight')
  for a,b in matched:
   c,p=(a,b) if a[0]=='c' else (b,a)
   one_to_one.append(dict(tool=tool,rule=rule,call_id=calls.loc[c[1],'call_id'],paper_id=paper.loc[p[1],'paper_id']))
  summary.append(dict(tool=tool,rule=rule,caller_total=len(calls),paper_total=len(paper),
    calls_with_overlap=chosen.call_id.nunique(),paper_with_overlap=chosen.paper_id.nunique(),overlap_pairs=len(chosen),
    concordant_calls=same.call_id.nunique(),concordant_paper=same.paper_id.nunique(),concordant_pairs=len(same),
    discordant_pairs=int((~chosen.same_direction).sum()),
    caller_overlap_fraction=same.call_id.nunique()/len(calls),paper_recovered_fraction=same.paper_id.nunique()/len(paper),
    one_to_one_concordant_matches=len(matched)))
  paper[f'{tool}_{rule}_same_direction']=paper.paper_id.isin(same.paper_id)
  calls[f'{rule}_same_direction']=calls.call_id.isin(same.call_id)
  for direction in ['hypo','hyper']:
   sub=same[same.paper_direction==direction]
   total_p=int((paper.direction==direction).sum());total_c=int((calls.direction==direction).sum())
   direction_summary.append(dict(tool=tool,rule=rule,direction=direction,caller_total=total_c,paper_total=total_p,
       concordant_calls=sub.call_id.nunique(),concordant_paper=sub.paper_id.nunique(),
       caller_overlap_fraction=sub.call_id.nunique()/total_c if total_c else np.nan,
       paper_recovered_fraction=sub.paper_id.nunique()/total_p))
  hit_genes=same.groupby('call_id').paper_gene.apply(lambda s:';'.join(sorted(set(s.dropna().astype(str)))))
  calls[f'{rule}_paper_genes']=calls.call_id.map(hit_genes).fillna('')
 for i,c in calls.iterrows():
  cp=paper[paper.chrom==c.chrom]
  gap=np.maximum(np.maximum(cp.start.to_numpy()-c.end-1,c.start-cp.end.to_numpy()-1),0)
  j=int(np.argmin(gap));nearest=cp.iloc[j]
  calls.loc[i,'nearest_paper_id']=nearest.paper_id
  calls.loc[i,'gap_to_nearest_paper_bp']=int(gap[j])
 call_tables[tool]=calls
 calls.to_csv(OUT/f'{tool}_calls_vs_paper.tsv',sep='\t',index=False)

summary=pd.DataFrame(summary);directions=pd.DataFrame(direction_summary);pairs=pd.concat(all_pairs,ignore_index=True)
summary.to_csv(OUT/'summary.csv',index=False)
directions.to_csv(OUT/'by_direction.csv',index=False)
pairs.to_csv(OUT/'overlap_pairs.tsv',sep='\t',index=False)
pd.DataFrame(one_to_one).to_csv(OUT/'one_to_one_matches.tsv',sep='\t',index=False)
paper.to_csv(OUT/'paper_regions_with_matches.tsv',sep='\t',index=False)
sets=[]
for rule in rules:
 for direction in ['all','hypo','hyper']:
  p=paper if direction=='all' else paper[paper.direction==direction]
  e=p[f'epykit_{rule}_same_direction'];d=p[f'DSS_{rule}_same_direction']
  sets.append(dict(rule=rule,direction=direction,total=len(p),both=int((e&d).sum()),epykit_only=int((e&~d).sum()),DSS_only=int((~e&d).sum()),neither=int((~e&~d).sum())))
sets=pd.DataFrame(sets);sets.to_csv(OUT/'paper_recovery_sets.csv',index=False)
shape=[]
for direction in ['all','hypo','hyper']:
 p=paper if direction=='all' else paper[paper.direction==direction]
 shape.append(dict(direction=direction,total=len(p),median_length_bp=float(p.length_bp.median()),
      below_50bp=int((p.length_bp<50).sum()),fewer_than_5_cpgs=int((p.n_cpgs<5).sum()),
      length_ge50_and_cpg_ge5=int(((p.length_bp>=50)&(p.n_cpgs>=5)).sum())))
pd.DataFrame(shape).to_csv(OUT/'paper_region_sizes.csv',index=False)
with pd.ExcelWriter(OUT/'paper_vs_epykit_DSS.xlsx',engine='openpyxl') as writer:
 summary.to_excel(writer,sheet_name='Summary',index=False)
 directions.to_excel(writer,sheet_name='By direction',index=False)
 sets.to_excel(writer,sheet_name='Paper recovery sets',index=False)
 paper.to_excel(writer,sheet_name='Paper regions',index=False)
 for tool,calls in call_tables.items():calls.to_excel(writer,sheet_name=tool+' calls',index=False)
 pairs.to_excel(writer,sheet_name='Overlap pairs',index=False)
 pd.DataFrame(one_to_one).to_excel(writer,sheet_name='One-to-one matches',index=False)
 for ws in writer.book.worksheets:
  ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
  for cells in ws.columns:
   width=min(42,max(12,max(len(str(c.value or '')) for c in cells[:40])+2))
   ws.column_dimensions[cells[0].column_letter].width=width
(OUT/'provenance.json').write_text(json.dumps(dict(paper_sha256=hashlib.sha256(XLSX.read_bytes()).hexdigest(),
    caller_sha256={t:hashlib.sha256((RUN/t/'dmr.tsv').read_bytes()).hexdigest() for t in call_tables},
    assembly='hg19',coordinates='1-based inclusive; paper lengths exactly match end-start+1',
    paper_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4665002/',
    excluded_rows=raw.loc[~valid,'Chr'].dropna().astype(str).tolist(),
    interpretation='Paper concordance, not truth-based precision/FDR. Expression columns are not DMR significance.',
    same_direction='MTB relative to NI; use paper Direction, never expression log2FCa',
    matching='Any overlap and reciprocal base-pair overlap; unique call/paper counts permit many-to-many. Separate maximum-cardinality one-to-one matches maximize total interval Jaccard as a tie-break.'),indent=2)+'\n')
print(summary.to_string(index=False));print('\nDIRECTION\n',directions[directions.rule.isin(['any_overlap','reciprocal_50pct'])].to_string(index=False))
print('\nRECOVERY SETS\n',sets[sets.rule.isin(['any_overlap','reciprocal_50pct'])].to_string(index=False))
print('\nSIZES\n',pd.DataFrame(shape).to_string(index=False))
