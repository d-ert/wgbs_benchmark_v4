"""Three sequential, immutable epykit runs with the new count engine."""
from pathlib import Path
import fcntl
import hashlib
import json
import shutil
import subprocess
import sys
import time

ROOT = Path('/scratch/wgbs_benchmark_v4')
SOURCE = Path('/scratch/epykit-calibration-fresh_20261007/src')
BASELINE = ROOT/'results/genome_baseline'
REGION_ONLY = ROOT/'results/default_search_by_20261007'
OUTPUT = ROOT/'results/beta_binomial_20261007_full'
HERE = Path(__file__).parent.resolve()
sys.path[:0] = [str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)]
import pandas as pd
from benchmark import contract, preflight
from analysis.run_publication_workflow import environment, flatten
from wgbs_v3.benchmark import run_monitored
from wgbs_v3.provenance import snapshot, source_inventory, write_json
from wgbs_v3.score import score


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if OUTPUT.exists():
        raise ValueError(f'Preserving existing output: {OUTPUT}')
    datasets = [ROOT/'data/genome_autosomes_signal',ROOT/'data/genome_autosomes_null']
    real = ROOT/'data/real'
    baseline = json.loads((BASELINE/'run.json').read_text())
    with (ROOT/'results/.workflow.lock').open('a') as lock, (ROOT/'results/.timing.lock').open('a') as timing:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        fcntl.flock(timing,fcntl.LOCK_EX|fcntl.LOCK_NB)
        print('Verifying frozen input and original harness hashes...',flush=True)
        frozen = contract(datasets,real)
        if frozen!=baseline['contract']:
            raise ValueError('Original inputs or harness changed')
        preflight(SOURCE)
        OUTPUT.mkdir()
        snapshot(ROOT,SOURCE,real,OUTPUT)
        source = OUTPUT/'source_snapshot'
        for name in ['dmr.py','_region_search.py']:
            assert digest(source/'epykit'/name)==digest(REGION_ONLY/'source_snapshot/epykit'/name),name
        shutil.copy2(HERE/'runner.py',source/'runner_beta_binomial.py')
        shutil.copy2(Path(__file__),source/'run_candidate.py')
        current = source_inventory(source/'epykit')
        previous = source_inventory(BASELINE/'source_snapshot/epykit')
        changed = sorted(k for k in set(current)|set(previous) if current.get(k)!=previous.get(k))
        commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=SOURCE.parent,text=True).strip()
        record = dict(status='running',kind='epykit',baseline=str(BASELINE),region_only=str(REGION_ONLY),
            tools=['epykit'],datasets=list(map(str,datasets)),real_input=str(real),contract=frozen,
            epykit_source=str(source),epykit_hashes=current,package_commit=commit,
            changed_package_files=changed,runner_sha256=digest(source/'runner_beta_binomial.py'),
            orchestrator_sha256=digest(source/'run_candidate.py'),
            region_correction='complete-interval-family-by-v1',
            model='per-replicate beta-binomial; count-likelihood EB log-rho MAP; approximate F reference',
            confidence_intervals=False,real_memory_limit_gib=None,
            note='Only epykit rerun. Competitor simulations reused; original real run partial.')
        write_json(OUTPUT/'run.json',record)
        print('Source frozen. Waiting for held-out validation and package checks before timing...',flush=True)
        record['stage'] = 'waiting_for_validation'
        write_json(OUTPUT/'run.json',record)
        while True:
            verification_path = HERE/'verification.json'
            if verification_path.exists() and json.loads(verification_path.read_text()).get('ready_for_benchmark'):
                break
            time.sleep(10)
        validation = json.loads((HERE/'holdout_v3/context.json').read_text())
        assert validation['status']=='complete'
        assert validation['source_sha256']==digest(source/'epykit/_beta_binomial.py')
        record.update(stage='benchmarking',heldout_validation=validation)
        write_json(OUTPUT/'run.json',record)
        env = environment(source)
        env['MPLCONFIGDIR'] = '/tmp/epykit_bb_benchmark_mpl'
        rows,checks = [],[]
        def run_one(dataset,target,memory_gib):
            target.mkdir(parents=True)
            cmd = [sys.executable,str(source/'runner_beta_binomial.py'),str(dataset),str(target)]
            monitoring = run_monitored(cmd,target/'run.log',env,memory_gib=memory_gib)
            if memory_gib==float('inf'):
                monitoring['memory_limit_gib'] = None
            result = dict(status='ok' if monitoring['exit_code']==0 and (target/'dmr.tsv').exists() else 'failed',
                          monitoring=monitoring,profile='beta_binomial_eb_F')
            if result['status']=='ok':
                a = json.loads((BASELINE/('real' if dataset==real else 'simulated/'+dataset.name)/
                                'epykit/runner_config.json').read_text())
                b = json.loads((target/'runner_config.json').read_text())
                original_dmc = dict(a['dmc'],test='beta_binomial',beta_binomial_ci=False)
                assert original_dmc==b['dmc'], 'Unexpected CpG configuration change'
                assert a['dmr']==b['dmr'], 'Region configuration changed'
                assert a['coverage_contract']==b['coverage_contract']
                assert a['assembly']==b['assembly']
                checks.append(dict(dataset=dataset.name,dmc_change=['test','beta_binomial_ci'],
                                   dmr_settings_identical=True,coverage_settings_identical=True,
                                   empirical_fdr=b['dmr']['empirical_fdr'],assembly=b['assembly']))
            return result
        try:
            for dataset in datasets:
                print(f'Running full {dataset.name}...',flush=True)
                target = OUTPUT/'simulated'/dataset.name
                result = run_one(dataset,target/'epykit',48.)
                manifest = dict(dataset=str(dataset),coordinate_system='1-based inclusive',
                                epykit_source=str(source),thread_budget=1,tools={'epykit':result})
                write_json(target/'manifest.json',manifest)
                if result['status']!='ok':
                    raise RuntimeError(f'{dataset.name}: epykit failed; inspect run.log')
                try:
                    result['score'] = score(dataset,target/'epykit/dml.tsv',target/'epykit/dmr.tsv',target/'epykit/score.json')
                except Exception as exc:
                    result.update(status='scoring_failed',scoring_error=f'{type(exc).__name__}: {exc}')
                    write_json(target/'manifest.json',manifest)
                    raise
                write_json(target/'manifest.json',manifest)
                metrics = flatten(result['score']); metrics.pop('dataset',None)
                row = dict(dataset=dataset.name,tool='epykit',status=result['status'],profile=result['profile'],**metrics,
                           wall_seconds=result['monitoring']['wall_seconds'],
                           peak_rss_bytes=result['monitoring']['peak_process_tree_rss_bytes'])
                rows.append(row)
                pd.DataFrame(rows).to_csv(OUTPUT/'simulated/simulation_metrics.csv',index=False)
                print(f'{dataset.name} completed and scored.',flush=True)
            print('Running all twelve GSE64177 samples, without a benchmark memory cap...',flush=True)
            result = run_one(real,OUTPUT/'real/epykit',float('inf'))
            result['assembly'] = 'hg19'
            write_json(OUTPUT/'real/run_manifest.json',{'epykit':result})
            if result['status']!='ok':
                raise RuntimeError('Real epykit failed; inspect run.log')
            write_json(OUTPUT/'regression_checks.json',checks)
            assert source_inventory(source/'epykit')==current
            assert contract(datasets,real)==frozen
            record['status'] = 'complete'
            record['stage'] = 'complete'
            write_json(OUTPUT/'run.json',record)
            print(f'Complete: {OUTPUT}',flush=True)
        except BaseException as exc:
            record.update(status='failed',error=f'{type(exc).__name__}: {exc}')
            write_json(OUTPUT/'run.json',record)
            raise


if __name__=='__main__':
    main()
