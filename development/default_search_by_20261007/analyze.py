"""Compare a fresh default epykit benchmark with all recorded original tools."""
from pathlib import Path
import json
import shutil
import subprocess
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
BASE = ROOT / 'results/genome_baseline'
RUN = ROOT / 'results/default_search_by_20261007'
sys.path[:0] = [str(ROOT / 'vendor310'), str(ROOT / 'src'), str(ROOT)]
import numpy as np
import pandas as pd
import networkx as nx


def main():
    assert json.loads((RUN / 'run.json').read_text())['status'] == 'complete'
    out = RUN / 'comparison'
    out.mkdir(exist_ok=True)
    old = pd.read_csv(BASE / 'simulated/simulation_metrics.csv')
    new = pd.read_csv(RUN / 'simulated/simulation_metrics.csv')
    old['source_run'] = BASE.name; old['reused'] = True
    new['source_run'] = RUN.name; new['reused'] = False
    pd.concat([old, new], ignore_index=True).to_csv(out / 'all_tool_simulation_metrics.csv', index=False)
    rows = []
    details = []
    for dataset in ['genome_autosomes_signal', 'genome_autosomes_null']:
        for label, origin in [('original', BASE), ('revised', RUN)]:
            score = json.loads((origin / 'simulated' / dataset / 'epykit/score.json').read_text())
            monitor = json.loads((origin / 'simulated' / dataset / 'manifest.json').read_text())['tools']['epykit']['monitoring']
            for rule in ['dmr', 'dmr_80']:
                rows.append(dict(dataset=dataset, version=label, matching=rule,
                                 **score[rule], runtime_s=monitor['wall_seconds'],
                                 peak_rss_gib=monitor['peak_process_tree_rss_bytes'] / 2**30))
            details.append(dict(dataset=dataset, version=label, **score['dml']))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out / 'epykit_region_comparison.csv', index=False)
    pd.DataFrame(details).to_csv(out / 'epykit_cpg_comparison.csv', index=False)
    changes = []
    for dataset in new.dataset:
        a = old[(old.dataset == dataset) & (old.tool == 'epykit')].iloc[0]
        b = new[new.dataset == dataset].iloc[0]
        for field in old.select_dtypes(include='number').columns.intersection(new.select_dtypes(include='number').columns):
            changes.append(dict(dataset=dataset, metric=field, baseline=a[field], candidate=b[field], delta=b[field]-a[field]))
    pd.DataFrame(changes).to_csv(out / 'epykit_changes.csv', index=False)

    real = out / 'real'; real.mkdir(exist_ok=True)
    sources = {'epykit_original': (BASE, 'epykit'), 'epykit_revised': (RUN, 'epykit'), 'DSS': (BASE, 'DSS')}
    real_perf = []
    for label, (origin, tool) in sources.items():
        dest = real / label; dest.mkdir(exist_ok=True)
        shutil.copy2(origin / 'real' / tool / 'dmr.tsv', dest / 'dmr.tsv')
        calls = pd.read_csv(dest / 'dmr.tsv', sep='\t')
        m = json.loads((origin / 'real/run_manifest.json').read_text())[tool]['monitoring']
        real_perf.append(dict(tool=label, source_run=origin.name, reused=origin == BASE, calls=len(calls),
                              runtime_s=m['wall_seconds'], peak_rss_gib=m['peak_process_tree_rss_bytes'] / 2**30))
    pd.DataFrame(real_perf).to_csv(out / 'real_performance.csv', index=False)
    subprocess.run(['Rscript', str(ROOT / 'analysis/annotate_real_dmrs.R'), str(real),
                    str(real / 'concordance'), ','.join(sources)], check=True)
    paper = pd.read_csv(ROOT / 'analysis/paper_dmr_comparison_20261006/paper_regions_with_matches.tsv', sep='\t')
    paper_rows = []; direction_rows = []; all_pairs = []; matches = []
    for label in sources:
        calls = pd.read_csv(real / label / 'dmr.tsv', sep='\t')
        calls['direction'] = np.where(calls.mean_meth_diff > 0, 'hyper', 'hypo') if label.startswith('epykit') else calls.direction.map({1:'hyper', -1:'hypo', 0:'unknown'})
        pairs = []
        for i, c in calls.iterrows():
            p = paper[(paper.chrom == c.chrom) & (paper.start <= c.end) & (paper.end >= c.start)]
            for j, r in p.iterrows():
                overlap = min(c.end, r.end) - max(c.start, r.start) + 1
                union = c.end-c.start+1 + r.end-r.start+1 - overlap
                pairs.append(dict(tool=label, call_index=i, paper_index=j, same_direction=c.direction == r.direction,
                                  overlap_bp=overlap, fraction_call=overlap/(c.end-c.start+1),
                                  fraction_paper=overlap/(r.end-r.start+1), jaccard=overlap/union))
        pairs = pd.DataFrame(pairs, columns=['tool','call_index','paper_index','same_direction','overlap_bp','fraction_call','fraction_paper','jaccard'])
        all_pairs.append(pairs)
        for cutoff in [0., .5, .8]:
            subset = pairs[(pairs.fraction_call >= cutoff) & (pairs.fraction_paper >= cutoff)]
            concordant = subset[subset.same_direction.astype(bool)]
            graph = nx.Graph()
            graph.add_weighted_edges_from((('c', int(r.call_index)), ('p', int(r.paper_index)), float(r.jaccard)) for r in concordant.itertuples())
            matching = nx.max_weight_matching(graph, maxcardinality=True, weight='weight')
            for a, b in matching:
                c, p = (a, b) if a[0] == 'c' else (b, a)
                matches.append(dict(tool=label, reciprocal_threshold=cutoff, call_index=c[1], paper_index=p[1]))
            paper_rows.append(dict(tool=label, reciprocal_threshold=cutoff, calls=len(calls), paper_total=len(paper),
                                   calls_overlapping_paper=concordant.call_index.nunique(),
                                   caller_overlap_fraction=concordant.call_index.nunique()/len(calls) if len(calls) else np.nan,
                                   paper_regions_covered=concordant.paper_index.nunique(), one_to_one_matches=len(matching),
                                   discordant_pairs=int((~subset.same_direction.astype(bool)).sum())))
            for direction in ['hypo', 'hyper']:
                selected = concordant[concordant.paper_index.isin(paper.index[paper.direction == direction])]
                direction_rows.append(dict(tool=label, reciprocal_threshold=cutoff, direction=direction,
                                            calls=int((calls.direction == direction).sum()),
                                            paper_total=int((paper.direction == direction).sum()),
                                            calls_overlapping_paper=selected.call_index.nunique(), paper_regions_covered=selected.paper_index.nunique()))
        calls['paper_overlap'] = calls.index.isin(pairs.loc[pairs.same_direction.astype(bool), 'call_index'])
        calls.to_csv(out / f'{label}_calls_vs_paper.tsv', sep='\t', index=False)
    pd.DataFrame(paper_rows).to_csv(out / 'paper_comparison.csv', index=False)
    pd.DataFrame(direction_rows).to_csv(out / 'paper_comparison_by_direction.csv', index=False)
    pd.concat(all_pairs, ignore_index=True).to_csv(out / 'paper_overlap_pairs.tsv', sep='\t', index=False)
    pd.DataFrame(matches).to_csv(out / 'paper_one_to_one_matches.tsv', sep='\t', index=False)
    # Region-only fix must retain the original boundaries for every surviving call.
    subsets = []
    for relative in ['simulated/genome_autosomes_signal/epykit', 'simulated/genome_autosomes_null/epykit', 'real/epykit']:
        a = pd.read_csv(BASE / relative / 'dmr.tsv', sep='\t')
        b = pd.read_csv(RUN / relative / 'dmr.tsv', sep='\t')
        keys = ['chrom','start','end']
        old_set = set(map(tuple, a[keys].to_numpy())); new_set = set(map(tuple, b[keys].to_numpy()))
        subsets.append(dict(dataset=relative, original=len(a), revised=len(b), preserved_boundaries=new_set <= old_set))
    assert all(row['preserved_boundaries'] for row in subsets)
    (out / 'boundary_regression.json').write_text(json.dumps(subsets, indent=2)+'\n')
    print(metrics[['dataset','version','matching','called_regions','matched_regions','region_precision','region_recall','region_f1','runtime_s','peak_rss_gib']].to_string(index=False))
    print(pd.DataFrame(paper_rows).to_string(index=False))


if __name__ == '__main__':
    main()
