"""Reproduce report evidence without altering callers, inputs or scoring."""
from pathlib import Path
import os, sys, json, hashlib, importlib.util, time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'vendor310'))
os.environ.setdefault('POLARS_MAX_THREADS', '1')
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/epykit-report-numba')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/epykit-report-mpl')
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def main():
    sources = {}
    def csv(path, name):
        p = ROOT / path
        df = pd.read_csv(p, keep_default_na=False, na_values=['', 'NA', 'NaN', 'nan'])
        sources[name] = {'file': path, 'sha256': sha(p), 'rows': len(df)}
        df.to_csv(OUT / (name+'.csv'), index=False)
        return df
    regions = csv('results/beta_binomial_20261007_full/comparison/all_tool_region_comparison.csv', 'regions')
    cpg = csv('results/beta_binomial_20261007_full/comparison/all_tool_cpg_comparison.csv', 'cpg')
    real = csv('results/real_uncapped_20261008/comparison/real_performance.csv', 'real_performance')
    paper = csv('results/real_uncapped_20261008/comparison/paper_comparison.csv', 'paper_agreement')
    heldout = csv('development/beta_binomial_20261007/holdout_v3/metrics.csv', 'heldout')
    strata = csv('development/epykit_calibration/reports/cpg_strata.csv', 'strata')
    assert len(strata[(strata.dataset=='null') & (strata.dimension=='coverage_bin')]) == 6
    csv('results/beta_binomial_20261007_full/comparison/genome_priors.csv','genome_priors')
    csv('development/epykit_calibration/reports/false_positive_distances.csv','false_positive_distances')
    csv('analysis/paper_dmr_comparison_20261006/paper_region_sizes.csv','paper_morphology')
    csv('analysis/genome_baseline_review_20261006/epykit_exploratory_threshold_tradeoff.csv','exploratory_thresholds')
    for row in regions.itertuples():
        if row.true_regions and row.called_regions:
            assert np.isclose(row.region_precision, row.matched_regions/row.called_regions)
            assert np.isclose(row.region_recall, row.matched_regions/row.true_regions)
            assert np.isclose(row.region_f1, 2*row.matched_regions/(row.called_regions+row.true_regions))
    assert (cpg.loc[cpg.available == True, ['tp','fp','fn','tn']].sum(axis=1).values == cpg.loc[cpg.available == True, 'n_eligible'].values).all()
    for r in paper.itertuples():
        if r.calls:
            assert np.isclose(r.overlap_fraction, r.calls_overlapping_paper/r.calls)
    timings, tails, equality, rawcounts = [], [], [], []
    runs = [('epykit_original','genome_baseline'), ('epykit_region_fix','default_search_by_20261007'), ('epykit_beta_binomial','beta_binomial_20261007_full')]
    for tool, run in runs:
        for dataset in ['genome_autosomes_signal','genome_autosomes_null','real']:
            base = ROOT/'results'/run/(Path('real/epykit') if dataset=='real' else Path('simulated')/dataset/'epykit')
            tpath = base/'timing.json'
            if tpath.exists():
                t = json.load(open(tpath))
                timings.append(dict(tool=tool,dataset=dataset,**t))
            dml = base/'dml.tsv'
            if tool != 'epykit_original' and dml.exists():
                summary = pl.scan_csv(dml, separator='\t').select(
                    pl.len().alias('sites'), pl.col('pvalue').min().alias('min_p'),
                    pl.col('qvalue').min().alias('min_q'),
                    (pl.col('qvalue')<=.05).sum().alias('bh_calls'),
                    (pl.col('pvalue')==0).sum().alias('p_zero'),
                    *[(pl.col('pvalue')<a).sum().alias('p_below_'+str(a)) for a in [.05,.001,1e-5,1e-7]]
                ).collect().to_dicts()[0]
                tails.append(dict(tool=tool,dataset=dataset,**summary))
            for kind in ['dmr','dmr_candidates']:
                path=base/(kind+'.tsv')
                if path.exists():
                    frame=pl.read_csv(path,separator='\t')
                    rawcounts.append(dict(tool=tool,dataset=dataset,kind=kind,rows=len(frame)))
            if tool == 'epykit_region_fix':
                original = ROOT/'results/genome_baseline'/(Path('real/epykit') if dataset=='real' else Path('simulated')/dataset/'epykit')/'dml.tsv'
                if original.exists():
                    h1,h2=sha(original),sha(dml)
                    equality.append(dict(dataset=dataset,original_sha256=h1,revised_sha256=h2,identical=h1==h2))
                    assert h1==h2
    pd.DataFrame(timings).to_csv(OUT/'stage_timings.csv',index=False)
    pd.DataFrame(tails).to_csv(OUT/'pvalue_tails.csv',index=False)
    pd.DataFrame(rawcounts).to_csv(OUT/'raw_counts.csv',index=False)
    rs = load_module('report_region_search',ROOT/'results/beta_binomial_20261007_full/source_snapshot/epykit/_region_search.py')
    family_rows=[]
    for dataset in ['genome_autosomes_signal','genome_autosomes_null','real']:
        base=ROOT/'results/default_search_by_20261007'/(Path('real/epykit') if dataset=='real' else Path('simulated')/dataset/'epykit')
        paths=list((base/'store/.cache/dmc/lr').glob('chrom=*.parquet'))
        if not paths: continue
        family=sum(rs.count_intervals(pl.read_parquet(p,columns=['pos'])['pos'].to_numpy(),5,50,500) for p in paths)
        # TSV named dmr_candidates is already filtered by tl.dmr; the hidden
        # native cache is the actual pre-q candidate family.
        caches=list((base/'store/.cache/dmc/lr').glob('.dmr_chain_merge*.parquet'))
        assert len(caches)==1
        candidate=pl.read_parquet(caches[0])
        q=rs.adjust_selected(candidate['combined_pvalue'].to_numpy(),family)
        assert np.allclose(q,candidate['combined_qvalue'].to_numpy().astype(float),rtol=1e-12,atol=1e-15,equal_nan=True)
        family_rows.append(dict(dataset=dataset,family_size=family,harmonic=float(__import__('scipy').special.digamma(family+1)+np.euler_gamma),candidates=len(candidate),calls=int((q<=.05).sum())))
    pd.DataFrame(family_rows).to_csv(OUT/'search_family.csv',index=False)
    # Geometry limits computed from the common eligibility catalogue; these
    # limits are about exact truth spans, not a ceiling on overlap matching.
    truth=pd.read_csv(ROOT/'data/genome_autosomes_signal/dmr_truth.tsv',sep='\t')
    geometry=[]
    for chrom,g in truth.groupby('chrom'):
        e=pl.read_parquet(ROOT/'data/genome_autosomes_signal/eligibility_parts'/(chrom+'.parquet'))
        pos=e.filter(pl.col('eligible'))['pos'].to_numpy()
        for r in g.itertuples():
            a=np.searchsorted(pos,r.start);b=np.searchsorted(pos,r.end,side='right')
            p=pos[a:b]
            geometry.append(dict(region_id=r.region_id,chrom=chrom,truth_cpgs=r.n_cpg,eligible_cpgs=len(p),eligible_span_bp=int(p[-1]-p[0]+1) if len(p) else 0,max_gap=int(np.diff(p).max()) if len(p)>1 else 0))
    pd.DataFrame(geometry).to_csv(OUT/'truth_geometry.csv',index=False)
    th=[]
    for nrep in [3,5,10]:
        df=2*nrep-2
        for lr in [10,20,30,50,100]:
            th.append(dict(nrep=nrep,df=df,lr=lr,p_F=float(stats.f.sf(lr,1,df)),p_chi2=float(stats.chi2.sf(lr,1))))
    pd.DataFrame(th).to_csv(OUT/'reference_tails.csv',index=False)
    states={}
    for path in (ROOT/'results').glob('*/run.json'):
        d=json.load(open(path));states[path.parent.name]={k:d.get(k) for k in ['status','stage','active_tool','started_at','updated_at']}
    source_files=['dmc.py','dmr.py','_beta_binomial.py','_region_search.py','_dmc_stages.py']
    hashes={}
    for fname in source_files:
        paths={label:root/fname for label,root in [
            ('original',ROOT/'results/genome_baseline/source_snapshot/epykit'),
            ('region_fix',ROOT/'results/default_search_by_20261007/source_snapshot/epykit'),
            ('bb_snapshot',ROOT/'results/beta_binomial_20261007_full/source_snapshot/epykit'),
            ('fresh_checkout',Path('/scratch/epykit-calibration-fresh_20261007/src/epykit')),
            ('epykit3_checkout',Path('/scratch/epykit3/src/epykit'))]}
        hashes[fname]={label:sha(p) if p.exists() else None for label,p in paths.items()}
    validation=dict(status='verified',sources=sources,run_states=states,cpg_byte_equality=equality,source_hashes=hashes,checks=['region metric arithmetic','CpG confusion totals','real overlap fractions','raw candidate row counts','complete-family BY reproduced','original vs region-fix CpG TSV equality'],script_sha256=sha(__file__))
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(dict(status='verified',timings=timings,tails=tails,search_families=family_rows,geometry_summary={'truth_regions':len(geometry),'fewer_than_5_eligible_cpgs':sum(r['eligible_cpgs']<5 for r in geometry),'eligible_span_below_50bp':sum(r['eligible_span_bp']<50 for r in geometry)},byte_equality=equality),indent=2))

if __name__=='__main__': main()
