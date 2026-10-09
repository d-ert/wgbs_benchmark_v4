"""Compare new epykit, both earlier epykit versions and all frozen tools."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
BASE = ROOT/'results/genome_baseline'
REGION = ROOT/'results/default_search_by_20261007'
RUN = ROOT/'results/beta_binomial_20261007_full'
HERE = Path(__file__).parent.resolve()
sys.path[:0] = [str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)]
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import chi2, f


def main():
    assert json.loads((RUN/'run.json').read_text())['status']=='complete'
    out = RUN/'comparison'; out.mkdir(exist_ok=True)
    analysis_snapshot = out/'source_snapshot'
    analysis_snapshot.mkdir(exist_ok=True)
    analysis_hashes = {}
    for name in ['analyze.py','report.py','plot_validation.py','validate_paper.R','validate_results.py','finish.py']:
        path = HERE/name
        shutil.copy2(path,analysis_snapshot/name)
        analysis_hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    versions = [('epykit_original',BASE,'epykit'),('epykit_region_fix',REGION,'epykit'),
                ('epykit_beta_binomial',RUN,'epykit')]
    versions += [(tool,BASE,tool) for tool in ['DSS','methylKit','BSmooth','dmrseq','DMRcate']]
    frames = []
    for origin,label in [(BASE,'original'),(REGION,'region_fix'),(RUN,'beta_binomial')]:
        frame = pd.read_csv(origin/'simulated/simulation_metrics.csv')
        frame['source_run'] = origin.name
        frame['version'] = np.where(frame.tool=='epykit',f'epykit_{label}',frame.tool)
        frame['reused'] = origin!=RUN
        frames.append(frame)
    pd.concat(frames,ignore_index=True).to_csv(out/'all_tool_simulation_metrics.csv',index=False)
    regions,sites,perf = [],[],[]
    for dataset in ['genome_autosomes_signal','genome_autosomes_null']:
        for label,origin,tool in versions:
            score = json.loads((origin/'simulated'/dataset/tool/'score.json').read_text())
            monitor = json.loads((origin/'simulated'/dataset/'manifest.json').read_text())['tools'][tool]['monitoring']
            for rule in ['dmr','dmr_80']:
                regions.append(dict(dataset=dataset,tool=label,matching=rule,source_run=origin.name,
                    reused=origin!=RUN,**score[rule],runtime_s=monitor['wall_seconds'],
                    peak_rss_gib=monitor['peak_process_tree_rss_bytes']/2**30))
            sites.append(dict(dataset=dataset,tool=label,source_run=origin.name,reused=origin!=RUN,
                              available='dml' in score,**score.get('dml',{})))
            perf.append(dict(dataset=dataset,tool=label,runtime_s=monitor['wall_seconds'],
                             peak_rss_gib=monitor['peak_process_tree_rss_bytes']/2**30))
    pd.DataFrame(regions).to_csv(out/'all_tool_region_comparison.csv',index=False)
    pd.DataFrame(sites).to_csv(out/'all_tool_cpg_comparison.csv',index=False)
    pd.DataFrame(perf).to_csv(out/'simulation_performance.csv',index=False)
    # Fixed-statistic diagnostic of the reference tails; no inference retuning.
    reference_rows = []
    for per_group in [2,3,5,10]:
        for statistic in [10.,20.,30.,50.]:
            p_chi = float(chi2.sf(statistic,1))
            p_f = float(f.sf(statistic,1,2*per_group-2))
            reference_rows.append(dict(replicates_per_group=per_group,
                reference_denominator_df=2*per_group-2,lr_statistic=statistic,
                chi2_reference_p=p_chi,F_reference_p=p_f,F_to_chi2_p_ratio=p_f/p_chi))
    pd.DataFrame(reference_rows).to_csv(out/'reference_tail_diagnostic.csv',index=False)
    # Other real tools never completed in the preserved baseline.
    sources = {label:(origin,tool) for label,origin,tool in versions[:3]+[('DSS',BASE,'DSS')]}
    real = out/'real'; real.mkdir(exist_ok=True)
    real_perf,paper_rows,direction_rows,all_pairs,matches = [],[],[],[],[]
    paper_path = ROOT/'analysis/paper_dmr_comparison_20261006/paper_regions_with_matches.tsv'
    paper = pd.read_csv(paper_path,sep='\t')
    assert len(paper)==3271 and int((paper.direction=='hypo').sum())==1714
    for label,(origin,tool) in sources.items():
        dest = real/label; dest.mkdir(exist_ok=True)
        shutil.copy2(origin/'real'/tool/'dmr.tsv',dest/'dmr.tsv')
        calls = pd.read_csv(dest/'dmr.tsv',sep='\t')
        monitor = json.loads((origin/'real/run_manifest.json').read_text())[tool]['monitoring']
        real_perf.append(dict(tool=label,source_run=origin.name,reused=origin!=RUN,calls=len(calls),
                             runtime_s=monitor['wall_seconds'],peak_rss_gib=monitor['peak_process_tree_rss_bytes']/2**30))
        calls['direction'] = (np.where(calls.mean_meth_diff>0,'hyper','hypo') if label.startswith('epykit')
                              else calls.direction.map({1:'hyper',-1:'hypo',0:'unknown'}))
        pairs = []
        for i,c in calls.iterrows():
            for j,p in paper[(paper.chrom==c.chrom)&(paper.start<=c.end)&(paper.end>=c.start)].iterrows():
                overlap = min(c.end,p.end)-max(c.start,p.start)+1
                union = c.end-c.start+1+p.end-p.start+1-overlap
                pairs.append(dict(tool=label,call_index=i,paper_index=j,same_direction=c.direction==p.direction,
                    overlap_bp=overlap,fraction_call=overlap/(c.end-c.start+1),
                    fraction_paper=overlap/(p.end-p.start+1),jaccard=overlap/union))
        pairs = pd.DataFrame(pairs,columns=['tool','call_index','paper_index','same_direction',
                                           'overlap_bp','fraction_call','fraction_paper','jaccard'])
        all_pairs.append(pairs)
        for cutoff in [0.,.5,.8]:
            subset = pairs[(pairs.fraction_call>=cutoff)&(pairs.fraction_paper>=cutoff)]
            concordant = subset[subset.same_direction.astype(bool)]
            graph = nx.Graph()
            graph.add_weighted_edges_from((('c',int(r.call_index)),('p',int(r.paper_index)),float(r.jaccard))
                                           for r in concordant.itertuples())
            matching = nx.max_weight_matching(graph,maxcardinality=True,weight='weight')
            for a,b in matching:
                c,p = (a,b) if a[0]=='c' else (b,a)
                matches.append(dict(tool=label,reciprocal_threshold=cutoff,call_index=c[1],paper_index=p[1]))
            paper_rows.append(dict(tool=label,reciprocal_threshold=cutoff,calls=len(calls),paper_total=len(paper),
                calls_overlapping_paper=concordant.call_index.nunique(),
                caller_overlap_fraction=concordant.call_index.nunique()/len(calls) if len(calls) else None,
                paper_regions_covered=concordant.paper_index.nunique(),one_to_one_matches=len(matching),
                discordant_pairs=int((~subset.same_direction.astype(bool)).sum())))
            for direction in ['hypo','hyper']:
                selected = concordant[concordant.paper_index.isin(paper.index[paper.direction==direction])]
                direction_rows.append(dict(tool=label,reciprocal_threshold=cutoff,direction=direction,
                    calls=int((calls.direction==direction).sum()),paper_total=int((paper.direction==direction).sum()),
                    calls_overlapping_paper=selected.call_index.nunique(),paper_regions_covered=selected.paper_index.nunique()))
        calls['paper_overlap'] = calls.index.isin(pairs.loc[pairs.same_direction.astype(bool),'call_index'])
        calls.to_csv(out/f'{label}_calls_vs_paper.tsv',sep='\t',index=False)
    pd.DataFrame(real_perf).to_csv(out/'real_performance.csv',index=False)
    pd.DataFrame(paper_rows).to_csv(out/'paper_comparison.csv',index=False)
    pd.DataFrame(direction_rows).to_csv(out/'paper_comparison_by_direction.csv',index=False)
    pd.concat(all_pairs,ignore_index=True).to_csv(out/'paper_overlap_pairs.tsv',sep='\t',index=False)
    pd.DataFrame(matches).to_csv(out/'paper_one_to_one_matches.tsv',sep='\t',index=False)
    subprocess.run(['Rscript',str(ROOT/'analysis/annotate_real_dmrs.R'),str(real),
                    str(real/'concordance'),','.join(sources)],check=True)
    # Inspect per-chromosome learned priors and check manifest integrity.
    prior_rows = []
    for dataset,target in [(d,RUN/'simulated'/d/'epykit') for d in
                           ['genome_autosomes_signal','genome_autosomes_null']]+[('GSE64177',RUN/'real/epykit')]:
        manifests = list(target.rglob('.epykit_dmc_manifest.json'))
        assert len(manifests)==1, (dataset,manifests)
        manifest = json.loads(manifests[0].read_text())
        verified_paths = set()
        for entry in manifest['chroms']:
            path = manifests[0].parent/entry['prior_file']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['prior_sha256']
            verified_paths.add(path.resolve())
        paths = sorted(target.rglob('priors/*/*.json'))
        assert {p.resolve() for p in paths}==verified_paths
        for path in paths:
            prior = json.loads(path.read_text())
            grid = np.asarray(prior['rho_grid'])
            log_grid = np.log(np.where(grid>0,grid,prior['binomial_representative']))
            centers = [float(np.asarray(w)@log_grid) for w in prior['start_rho_weights']]
            prior_rows.append(dict(dataset=dataset,prior_file=str(path.relative_to(RUN)),
                prior_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                n_groups=prior['n_groups'],converged=prior['converged'],
                log_rho_mean=prior['log_rho_mean'],log_rho_sd=prior['log_rho_sd'],
                start_loglik_difference=max(prior['start_logliks'])-min(prior['start_logliks']),
                start_log_rho_mean_difference=max(centers)-min(centers)))
    pd.DataFrame(prior_rows).to_csv(out/'genome_priors.csv',index=False)
    provenance = dict(status='complete',runs={k:str(p) for k,p in
        [('original',BASE),('region_fix',REGION),('beta_binomial',RUN)]},
        real_completed_tools=list(sources),real_missing_tools=['methylKit','BSmooth','dmrseq','DMRcate'],
        paper_workbook_sha256=hashlib.sha256((ROOT/'data/real/paper_reported_dmrs.xlsx').read_bytes()).hexdigest(),
        normalized_paper_sha256=hashlib.sha256(paper_path.read_bytes()).hexdigest(),
        simulation_match='one-to-one, 50% and 80% reciprocal eligible-CpG overlap',
        paper_match='direction-concordant, 1-based inclusive base-pair overlap',
        paper_is_ground_truth=False,genome_results_used_to_tune_engine=False,
        reference_tail_diagnostic='Fixed statistics only; does not establish chi-square calibration',
        analysis_source_sha256=analysis_hashes,
        heldout_context_sha256=hashlib.sha256((HERE/'holdout_v3/context.json').read_bytes()).hexdigest(),
        heldout_metrics_sha256=hashlib.sha256((HERE/'holdout_v3/metrics.csv').read_bytes()).hexdigest())
    (out/'sources.json').write_text(json.dumps(provenance,indent=2)+'\n')
    with pd.ExcelWriter(out/'comparison.xlsx') as writer:
        for filename,name in [('all_tool_region_comparison.csv','Simulation regions'),
                              ('all_tool_cpg_comparison.csv','Simulation CpGs'),
                              ('real_performance.csv','Real performance'),('paper_comparison.csv','Paper agreement'),
                              ('paper_comparison_by_direction.csv','Paper directions'),
                              ('paper_overlap_pairs.tsv','Paper overlap pairs'),('genome_priors.csv','Learned priors'),
                              ('reference_tail_diagnostic.csv','Reference tails')]:
            pd.read_csv(out/filename,sep='\t' if filename.endswith('.tsv') else ',').to_excel(writer,sheet_name=name,index=False)
        pd.read_csv(HERE/'holdout_v3/metrics.csv').to_excel(writer,sheet_name='Held-out calibration',index=False)
        pd.read_csv(HERE/'holdout_v3/prior_start_sensitivity.csv').to_excel(writer,sheet_name='Prior sensitivity',index=False)
    print(pd.DataFrame(regions)[['dataset','tool','matching','called_regions','matched_regions',
                                'region_precision','region_recall','region_f1','runtime_s','peak_rss_gib']].to_string(index=False))
    print(pd.DataFrame(paper_rows).to_string(index=False))


if __name__=='__main__':
    main()
