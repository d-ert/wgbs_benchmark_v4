"""Compare completed real calls with the paper and the recorded epykit versions."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT/'results/real_uncapped_20261008'
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor310'))
import numpy as np
import pandas as pd
import networkx as nx

def main():
    out = RUN/'comparison';out.mkdir(exist_ok=True)
    real = out/'real';real.mkdir(exist_ok=True)
    baseline = ROOT/'results/genome_baseline'
    fixed = ROOT/'results/default_search_by_20261007'
    bb = ROOT/'results/beta_binomial_20261007_full'
    sources = {'epykit':(fixed,'epykit'),'epykit_original':(baseline,'epykit'),
               'epykit_beta_binomial':(bb,'epykit'),'DSS':(baseline,'DSS')}
    status = json.loads((RUN/'real/run_manifest.json').read_text())
    sources.update({tool:(RUN,tool) for tool,entry in status.items() if entry['status']=='ok'})
    paper = pd.read_csv(ROOT/'analysis/paper_dmr_comparison_20261006/paper_regions_with_matches.tsv',sep='\t')
    performance,summary,all_pairs,matched,directions = [],[],[],[],[]
    for label,(origin,tool) in sources.items():
        target = real/label;target.mkdir(exist_ok=True)
        shutil.copy2(origin/'real'/tool/'dmr.tsv',target/'dmr.tsv')
        calls = pd.read_csv(target/'dmr.tsv',sep='\t')
        monitor = json.loads((origin/'real/run_manifest.json').read_text())[tool]['monitoring']
        performance.append(dict(tool=label,calls=len(calls),minutes=monitor['wall_seconds']/60,
            peak_rss_gib=monitor['peak_process_tree_rss_bytes']/2**30,source_run=origin.name,reused=origin!=RUN))
        calls['direction'] = (np.where(calls.mean_meth_diff>0,'hyper','hypo') if 'mean_meth_diff' in calls
                              else calls.direction.map({1:'hyper',-1:'hypo',0:'unknown'}))
        pairs = []
        for i,c in calls.iterrows():
            subset = paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)]
            for j,p in subset.iterrows():
                overlap = min(c.end,p.end)-max(c.start,p.start)+1
                pairs.append(dict(tool=label,call_index=i,paper_index=j,same_direction=c.direction==p.direction,
                    fraction_call=overlap/(c.end-c.start+1),fraction_paper=overlap/(p.end-p.start+1),overlap_bp=overlap))
        pairs = pd.DataFrame(pairs,columns=['tool','call_index','paper_index','same_direction','fraction_call','fraction_paper','overlap_bp'])
        all_pairs.append(pairs)
        for cutoff in [0.,.5,.8]:
            selected = pairs[(pairs.fraction_call>=cutoff)&(pairs.fraction_paper>=cutoff)]
            concordant = selected[selected.same_direction.astype(bool)]
            graph = nx.Graph()
            graph.add_edges_from((('c',int(r.call_index)),('p',int(r.paper_index))) for r in concordant.itertuples())
            matching = nx.max_weight_matching(graph,maxcardinality=True)
            for a,b in matching:
                c,p = (a,b) if a[0]=='c' else (b,a)
                matched.append(dict(tool=label,reciprocal_threshold=cutoff,call_index=c[1],paper_index=p[1]))
            summary.append(dict(tool=label,reciprocal_threshold=cutoff,calls=len(calls),paper_total=len(paper),
                calls_overlapping_paper=concordant.call_index.nunique(),paper_regions_covered=concordant.paper_index.nunique(),
                overlap_fraction=concordant.call_index.nunique()/len(calls) if len(calls) else None,
                one_to_one_matches=len(matching),discordant_pairs=int((~selected.same_direction.astype(bool)).sum())))
            for direction in ['hypo','hyper']:
                hit = concordant[concordant.paper_index.isin(paper.index[paper.direction==direction])]
                directions.append(dict(tool=label,reciprocal_threshold=cutoff,direction=direction,
                    calls=int((calls.direction==direction).sum()),calls_overlapping_paper=hit.call_index.nunique(),
                    paper_regions_covered=hit.paper_index.nunique()))
        calls['paper_overlap'] = calls.index.isin(pairs.loc[pairs.same_direction.astype(bool),'call_index'])
        calls.to_csv(out/f'{label}_calls_vs_paper.tsv',sep='\t',index=False)
    tables = {'real_performance':pd.DataFrame(performance),'paper_comparison':pd.DataFrame(summary),
              'paper_directions':pd.DataFrame(directions),'paper_overlap_pairs':pd.concat(all_pairs,ignore_index=True),
              'paper_one_to_one_matches':pd.DataFrame(matched)}
    for name,frame in tables.items():frame.to_csv(out/f'{name}.csv',index=False)
    subprocess.run(['Rscript',str(HERE/'annotate_real_dmrs.R'),str(real),str(real/'concordance'),','.join(sources)],check=True)
    subprocess.run(['Rscript',str(HERE/'validate_paper.R')],check=True)
    with pd.ExcelWriter(out/'real_comparison.xlsx') as writer:
        for name,frame in tables.items():frame.to_excel(writer,sheet_name=name[:31],index=False)
    view = pd.DataFrame(performance).merge(pd.DataFrame(summary).query('reciprocal_threshold==0'),on=['tool','calls'])
    text = '# GSE64177 real-data comparison\n\n'+view[['tool','calls','calls_overlapping_paper','overlap_fraction','paper_regions_covered','minutes','peak_rss_gib']].to_markdown(index=False,floatfmt='.3f')
    text += '\n\nEpykit uses the region-only correction and legacy CpG engine. BB remains experimental. '
    text += 'New R calls use the recorded statistical settings and all twelve samples with no benchmark memory or time cap. '
    text += 'The models omit donor pairing. Paper agreement is descriptive, not biological precision or FDR. '
    text += 'Matches use direction-concordant, 1-based inclusive hg19 intervals at any, 50%, and 80% reciprocal base-pair overlap. '
    text += 'Four-call-set historical comparisons are retained as separate source runs. Failed calls are excluded and recorded in the run manifest.\n'
    (out/'README.md').write_text(text)
    (out/'sources.json').write_text(json.dumps({k:{'run':str(p),'tool':t} for k,(p,t) in sources.items()},indent=2)+'\n')
    print(text,flush=True)

if __name__=='__main__':main()
