"""Build the report's reviewed tables, workbook and visual-app snapshot."""
from pathlib import Path
import sys, json, re, datetime
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor310'))
import numpy as np
import pandas as pd

LABELS={'epykit_original':'epykit original','epykit_region_fix':'epykit corrected','epykit_beta_binomial':'epykit BB-F','epykit':'epykit corrected','DSS':'DSS','dmrseq':'dmrseq','DMRcate':'DMRcate','methylKit':'methylKit','BSmooth':'BSmooth'}
ORDER=['dmrseq','epykit_region_fix','DSS','epykit_beta_binomial','epykit_original','BSmooth','methylKit','DMRcate']
def clean(df):return json.loads(df.to_json(orient='records'))
def table_md(df):
    def fmt(v):
        if v is None or (isinstance(v,float) and np.isnan(v)):return '—'
        if isinstance(v,(int,np.integer)):return f'{v:,}'
        if isinstance(v,float):return f'{v:.3f}' if abs(v)<1 else f'{v:.2f}'
        return str(v).replace('|',' / ')
    return '| '+' | '.join(df.columns)+' |\n| '+ ' | '.join(['---']*len(df.columns))+' |\n'+ '\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))

def main():
    derived={'region_comparison','cpg_comparison','real_comparison','candidate_comparison','timing_comparison','gse263850_comparison','roadmap'}
    raw={p.stem:pd.read_csv(p,keep_default_na=False,na_values=['','NA','NaN','nan']) for p in OUT.glob('*.csv') if p.stem not in derived}
    r=raw['regions'];signal=r[(r.dataset=='genome_autosomes_signal')&(r.matching=='dmr')].set_index('tool').loc[ORDER]
    r80=r[(r.dataset=='genome_autosomes_signal')&(r.matching=='dmr_80')].set_index('tool')
    null=r[(r.dataset=='genome_autosomes_null')&(r.matching=='dmr')].set_index('tool')
    region_table=pd.DataFrame([{'Tool':LABELS[t],'Calls':int(x.called_regions),'Matches 50%':int(x.matched_regions),'Precision (%)':100*x.region_precision,'Recall (%)':100*x.region_recall,'F1 50%':x.region_f1,'F1 80%':r80.loc[t,'region_f1'],'Null regions':int(null.loc[t,'called_regions']),'Minutes':x.runtime_s/60,'Peak RSS (GiB)':x.peak_rss_gib} for t,x in signal.iterrows()])
    c=raw['cpg'];cs=c[(c.dataset=='genome_autosomes_signal')&(c.available==True)].set_index('tool');cn=c[(c.dataset=='genome_autosomes_null')&(c.available==True)].set_index('tool')
    cpg_table=pd.DataFrame([{'Tool':LABELS[t],'True sites':int(x.tp),'False sites':int(x.fp),'Recall (%)':100*x.recall_eligible,'Observed FDP (%)':100*x.false_discovery_proportion if x.tp+x.fp else None,'Average precision':x.average_precision,'AUROC':x.auroc,'Null sites':int(cn.loc[t,'fp'])} for t,x in cs.iterrows()])
    p=raw['paper_agreement'];p0=p[p.reciprocal_threshold==0].set_index('tool');p50=p[p.reciprocal_threshold==.5].set_index('tool');real=raw['real_performance'].set_index('tool')
    real_table=pd.DataFrame([{'Tool':LABELS[t],'Calls':int(x.calls),'Paper-overlapping calls':int(x.calls_overlapping_paper),'Agreement (%)':100*x.overlap_fraction if x.calls else None,'Paper regions covered':int(x.paper_regions_covered),'Matches 50%':int(p50.loc[t,'one_to_one_matches']),'Minutes':real.loc[t,'minutes'],'Peak RSS (GiB)':real.loc[t,'peak_rss_gib']} for t,x in p0.iterrows()])
    gp=raw['gse263850_performance'].set_index('tool');gm=raw['gse263850_paper'];g0=gm[gm.reciprocal_threshold==0].set_index('tool');g50=gm[gm.reciprocal_threshold==.5].set_index('tool')
    gse_table=pd.DataFrame([{'Tool':LABELS[t],'Calls':int(x.calls),'Any-overlap calls':int(x.calls_overlapping_paper),'Matches 50%':int(g50.loc[t,'one_to_one_matches']),'Reference recovery 50% (%)':100*g50.loc[t,'one_to_one_recall'],'Minutes':gp.loc[t,'minutes'],'Peak RSS (GiB)':gp.loc[t,'peak_rss_gib']} for t,x in g0.iterrows()])
    cand=raw['candidate_ceiling'];cand=cand[cand.threshold==.5]
    candidate_table=pd.DataFrame([{'Engine':{'original':'original','region_fix':'corrected','BB_F':'BB-F'}[x.tool],'Stage':{'before_q':'Before q filter','after_q':'After q filter'}[x.stage],'Calls':int(x.calls),'Matches':int(x.matched),'Recall (%)':100*x.recall} for x in cand.itertuples()])
    st=raw['stage_timings'];st=st[st.dataset=='genome_autosomes_signal']
    timing_table=pd.DataFrame([{'Engine':LABELS[x.tool],'Read/filter (s)':x.read_filter_s,'CpG analysis (s)':x.dml_s,'Region analysis (s)':x.dmr_s,'Phase total (s)':x.read_filter_s+x.dml_s+x.dmr_s} for x in st.itertuples()])
    roadmap_table=pd.DataFrame([
        {'Priority':'1 · small','Change':'Document actual F tails; export pre-q candidates; expose inference/source revision','Expected benefit':'Reliable diagnostics and reproducibility','Evidence':'Confirmed comments/output mismatch','Validation':'Metadata and candidate-cache/output agreement'},
        {'Priority':'1 · statistical','Change':'Fit biological dispersion with depth/mean-aware count modeling and valid small-sample inference','Expected benefit':'Better calibrated site detection and usable power','Evidence':'Depth-stratified null excess; BB extreme-tail collapse','Validation':'Untouched factorial signal/null studies; tails and power together'},
        {'Priority':'2 · numerical','Change':'Reduce nested rho/mean refits; adaptive derivatives, warm starts and bounded parallelism','Expected benefit':'Reduce the principal BB compute cost','Evidence':'29× DMC phase increase; 47 objectives/site in bounded probe','Validation':'Count-likelihood oracle, boundary cases; repeated warm/cold timings'},
        {'Priority':'2 · regions','Change':'Coordinate-only multiscale family or discovery-aware region resampling','Expected benefit':'Recover weak coherent signals at valid error rates','Evidence':'41.8% candidate ceiling; 776M-interval penalty','Validation':'Selection replay; 50/80% matching; correlated nulls'},
        {'Priority':'2 · design','Change':'Support paired/covariate count models and appropriate resampling','Expected benefit':'Use paired real design; control residual/confounding variation','Evidence':'Six donor pairs omitted; BB rejects formulas','Validation':'Paired simulations and predeclared donor-adjusted real comparisons'},
        {'Priority':'3 · memory','Change':'Load count blocks; avoid unnecessary full-table materialization; exact disk-backed correction','Expected benefit':'Scale to larger cohorts and reduce memory/disk cost','Evidence':'Full chromosome×sample arrays; full p-value sort vector','Validation':'Peak RSS/disk versus sites and replicates; unchanged values'},
        {'Priority':'3 · major','Change':'Spatial beta-binomial regime/change-point backend','Expected benefit':'Joint spatial detection, boundaries and uncertainty','Evidence':'Literature-supported research direction; not benchmarked here','Validation':'Matching hypotheses, spatial misspecification, particle/state sensitivity'},
        {'Priority':'Conditional · major','Change':'Count mixed model with kinship/repeated-subject covariance','Expected benefit':'Handle structured cohorts when required','Evidence':'MACAU precedent; no kinship cause established in current simulation','Validation':'Structured-population simulations, identifiability, runtime scaling'},
    ])
    tables={'REGION_TABLE':('region_comparison',region_table),'CPG_TABLE':('cpg_comparison',cpg_table),'REAL_TABLE':('real_comparison',real_table),'GSE_TABLE':('gse263850_comparison',gse_table),'CANDIDATE_TABLE':('candidate_comparison',candidate_table),'TIMING_TABLE':('timing_comparison',timing_table),'ROADMAP_TABLE':('roadmap',roadmap_table)}
    template=OUT/'report.template.md'
    if not template.exists():template.write_text((OUT/'report.md').read_text())
    prose=template.read_text();rendered=prose
    for marker,(key,df) in tables.items():
        df.to_csv(OUT/(key+'.csv'),index=False)
        rendered=rendered.replace('<!-- '+marker+' -->',table_md(df))
    assert not re.search(r'<!-- [A-Z_]+ -->',rendered)
    (OUT/'report.md').write_text(rendered)
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    queries={}
    for name,df in raw.items():
        queries[name]={'rows':clean(df),'source':{'title':name.replace('_',' ').capitalize(),'files':[str((OUT/(name+'.csv')).relative_to(ROOT))],'executedAt':now,'caveats':['Saved benchmark/development measurements; each source run and metric unit retained.'], 'evidenceFlow':[{'title':'Read local evidence','detail':str((OUT/(name+'.csv')).relative_to(ROOT))},{'title':'Reproduce','detail':'analyze.py and candidate_ceiling.py; read validation.json for checks and authoritative input hashes.'}]}}
    for marker,(name,df) in tables.items():
        queries[name]={'rows':clean(df),'source':{'title':name.replace('_',' ').capitalize(),'files':[str((OUT/(name+'.csv')).relative_to(ROOT))],'executedAt':now,'caveats':['Comparison reflects recorded workflow settings, not equal realized error rates. Null and undefined values are retained.'],'evidenceFlow':[{'title':'Derived comparison','detail':'build_report.py selects signal or real rows, converts seconds to minutes and fractions to percentages; raw measures retained in accompanying CSVs.'}]}}
    qmap=[['regions','cpg','real_performance'],['regions','cpg','truth_geometry'],['regions'],['cpg','pvalue_tails'],['paper_agreement','real_performance','gse263850_paper','gse263850_performance'],[],['search_family','candidate_ceiling'],['candidate_ceiling','truth_geometry'],['strata'],['heldout','pvalue_tails','reference_tails'],[],['stage_timings'],['roadmap'],[],[],[]]
    pieces=re.split(r'^## ',prose,flags=re.M)
    sections=[]
    for i,piece in enumerate(pieces[1:]):
        text='## '+piece
        # Source anchor table is better represented as prose in the app;
        # the Markdown report retains the exact table and local evidence links.
        if i==15:
            text=re.sub(r'\| Finding \| Source anchor \|.*?(?=\nRun the report diagnostics)', 'Detailed source anchors are recorded in the companion Markdown report.\n',text,flags=re.S)
        parts=re.split(r'(<!-- [A-Z_]+ -->)',text)
        blocks=[]
        for j,part in enumerate(parts):
            match=re.fullmatch(r'<!-- ([A-Z_]+) -->',part)
            if match:
                marker=match.group(1);name,df=tables[marker]
                blocks.append({'kind':'table','id':name,'queryId':name,'columns':list(df.columns)})
            elif part.strip():
                blocks.append({'kind':'narrative','id':f'section-{i}-prose-{j}','value':part.strip()})
        sections.append({'id':f'section-{i}','title':piece.split('\n')[0],'queryIds':qmap[i] if i<len(qmap) else [],'blocks':blocks})
    queries['report_sections']={'rows':sections,'payloadColumns':['blocks'],'source':{'title':'Technical narrative','files':['analysis/epykit_final_report_20261009/report.md'],'caveats':['Observed mechanisms, illustrative calculations and untested improvements are distinguished in the text.']}}
    snapshot={'title':'Epykit is fast, but calibrated inference still limits its power','status':'reviewed','buildStatus':'creating','generatedAt':now,'surface':'report','report':{'asOf':'2026-10-09'},'filters':[],'queries':queries}
    (OUT/'reviewed.json').write_text(json.dumps(snapshot,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    with pd.ExcelWriter(OUT/'evidence.xlsx',engine='openpyxl') as w:
        for name,df in {**raw,**{name:df for name,df in tables.values()}}.items():df.to_excel(w,sheet_name=name[:31],index=False)
    print(json.dumps({'sections':len(sections),'queries':len(queries),'report_words':len(rendered.split()),'snapshot':str(OUT/'reviewed.json'),'workbook':str(OUT/'evidence.xlsx')}))

if __name__=='__main__':main()
