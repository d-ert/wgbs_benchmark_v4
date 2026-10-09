import React from 'react';
import { DataComponent, DataTable, EvidenceChart, ReportSection, RichNarrative, useDataApp } from '../../data-app-public.jsx';
const regionSpec={type:'horizontalBar',x:'Tool',y:'F1 50%',fields:['F1 50%','F1 80%'],stackable:false,valueDecimals:3,legend:{labels:{'F1 50%':'50% overlap','F1 80%':'80% overlap'}}};
const timingSpec={type:'horizontalStackedBar',x:'Engine',y:'Read/filter (s)',fields:['Read/filter (s)','CpG analysis (s)','Region analysis (s)'],valueDecimals:1,xLabel:'Seconds'};
const candidateSpec={type:'horizontalBar',x:'Engine',y:'Recall (%)',series:'Stage',stackable:false,valueDecimals:1};
export function ReportContent(){
  const {snapshot,appTitle,reviewedRows,visible,canEdit,mode,setAppTitle}=useDataApp();
  const sourcePreviews={
    'https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html':{title:'DSS official guide',summary:'Describes count-level dispersion shrinkage, smoothing and multifactor differential methylation analysis.',approvedForReport:true},
    'https://github.com/kdkorthauer/dmrseq':{title:'Dmrseq official repository',summary:'Describes de novo region discovery and permutation-based region inference.',approvedForReport:true},
    'https://arxiv.org/pdf/2211.07311':{title:'Hirt et al., methylome change points',summary:'A joint case/control regime model with particle inference and hypothesis-specific decision rules.',approvedForReport:true},
    'https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005650':{title:'MACAU original paper',summary:'A count mixed model incorporating sample covariance and overdispersion; the reported gain is specific to the studied cohort and threshold.',approvedForReport:true},
  };
  function narrative(b){
    const value=b.value.replace(/\[([^\]]+)\]\((?!https?:)([^)]+)\)/g,'$1');
    return <RichNarrative key={b.id} id={b.id} value={value} sourcePreviews={sourcePreviews} label="Edit report text"/>;
  }
  function table(b){
    const rows=reviewedRows(b.queryId);
    const columns=b.columns.map(field=>({field,label:field,renderCell:value=>typeof value==='number'?(field.startsWith('F1')||field==='Average precision'||field==='AUROC'?value.toFixed(3):Number.isInteger(value)?value.toLocaleString():value.toFixed(2)):value??'—'}));
    return <DataComponent key={b.id} id={b.id} title={b.queryId.replaceAll('_',' ')} queryId={b.queryId} kind="table" displayRows={rows} sourceRows={rows}><DataTable rows={rows} columns={columns} caption={b.queryId.replaceAll('_',' ')} searchable={false} compactNumbers={false}/></DataComponent>;
  }
  function section(s){
    const blocks=s.blocks.map(b=>b.kind==='table'?table(b):narrative(b));
    const chart=s.id==='section-2'?<EvidenceChart id="region-boundaries-chart" queryId="region_comparison" title="The region ranking changes with the overlap requirement" spec={regionSpec} rows={reviewedRows('region_comparison')} sourceRows={reviewedRows('region_comparison')} height={430}/>:
      s.id==='section-7'?<EvidenceChart id="candidate-recall-chart" queryId="candidate_comparison" title="Candidate detection limits recall before correction" spec={candidateSpec} rows={reviewedRows('candidate_comparison')} sourceRows={reviewedRows('candidate_comparison')} height={310}/>:
      s.id==='section-11'?<EvidenceChart id="runtime-stage-chart" queryId="timing_comparison" title="CpG analysis accounts for most of the BB runtime" spec={timingSpec} rows={reviewedRows('timing_comparison')} sourceRows={reviewedRows('timing_comparison')} height={300}/>:null;
    const qs=s.queryIds,inside=<>{blocks}{chart}</>;
    return qs.length?<ReportSection key={s.id} id={s.id} title={s.title} queryId={qs[0]} queryIds={qs} sourceRowsByQuery={Object.fromEntries(qs.map(q=>[q,reviewedRows(q)]))} showHeading={false}>{inside}</ReportSection>:<section key={s.id} id={s.id} className="report-section">{inside}</section>;
  }
  return <article className="report-content epykit-report" aria-label="Epykit comprehensive performance report"><header className="report-hero"><h1 data-data-app-title contentEditable={canEdit&&mode==='edit'} suppressContentEditableWarning onBlur={e=>{if(canEdit&&mode==='edit')setAppTitle(e.currentTarget.textContent.trim()||appTitle)}}>{appTitle}</h1><RichNarrative id="report-evidence-date" value="Technical performance review · completed benchmark evidence available 9 October 2026" className="report-deck"/></header>{snapshot.queries.report_sections.rows.filter(s=>visible(s.id)).map(section)}</article>;
}
