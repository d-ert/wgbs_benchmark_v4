"""Default-path epykit rerun against the preserved, partially completed baseline."""
from pathlib import Path
import fcntl
import hashlib
import json
import shutil
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
SOURCE = Path('/scratch/epykit-calibration-fresh_20261007/src')
BASELINE = ROOT / 'results/genome_baseline'
OUTPUT = ROOT / 'results/default_search_by_20261007'
sys.path[:0] = [str(ROOT / 'vendor310'), str(ROOT / 'src'), str(ROOT)]

from benchmark import contract, preflight
from analysis.run_publication_workflow import simulated, run_real
from wgbs_v3.provenance import snapshot, source_inventory, write_json


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    if OUTPUT.exists():
        raise ValueError(f'Preserving existing output: {OUTPUT}')
    datasets = [ROOT / 'data/genome_autosomes_signal', ROOT / 'data/genome_autosomes_null']
    real = ROOT / 'data/real'
    baseline = json.loads((BASELINE / 'run.json').read_text())
    with (ROOT / 'results/.workflow.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print('Checking original input and harness hashes...', flush=True)
        frozen = contract(datasets, real)
        if frozen != baseline['contract']:
            raise ValueError('Original benchmark inputs or harness changed')
        preflight(SOURCE)
        OUTPUT.mkdir()
        provenance = snapshot(ROOT, SOURCE, real, OUTPUT)
        source = OUTPUT / 'source_snapshot'
        changed = []
        original = source_inventory(BASELINE / 'source_snapshot/epykit')
        current = source_inventory(source / 'epykit')
        for name in set(original) | set(current):
            if original.get(name) != current.get(name):
                changed.append(name)
        if set(changed) != {'dmr.py', 'tl.py', '_region_search.py'}:
            raise ValueError(f'Unexpected package changes: {changed}')
        shutil.copy2(Path(__file__), source / 'run_candidate.py')
        record = dict(status='running', kind='epykit', baseline=str(BASELINE),
                      baseline_status=baseline['status'], contract=frozen,
                      epykit_source=str(source), epykit_hashes=current,
                      changed_package_files=sorted(changed),
                      region_correction='complete-interval-family-by-v1',
                      note='Default runner settings; frozen competitors reused; real baseline only epykit/DSS completed.')
        write_json(OUTPUT / 'run.json', record)
        try:
            print('Running epykit on both full simulations...', flush=True)
            simulated(datasets, OUTPUT / 'simulated', ('epykit',), source)
            print('Running epykit on all twelve GSE64177 samples...', flush=True)
            run_real(real, OUTPUT / 'real', ('epykit',), source, [])
            checks = []
            for relative in ['simulated/genome_autosomes_signal/epykit',
                             'simulated/genome_autosomes_null/epykit', 'real/epykit']:
                old, new = BASELINE / relative, OUTPUT / relative
                a = json.loads((old / 'runner_config.json').read_text())
                b = json.loads((new / 'runner_config.json').read_text())
                checks.append(dict(dataset=relative,
                                   dml_identical=digest(old / 'dml.tsv') == digest(new / 'dml.tsv'),
                                   dmc_settings_identical=a['dmc'] == b['dmc'],
                                   dmr_settings_identical=a['dmr'] == b['dmr'],
                                   empirical_fdr=b['dmr']['empirical_fdr']))
            write_json(OUTPUT / 'regression_checks.json', checks)
            assert all(c['dml_identical'] and c['dmc_settings_identical']
                       and c['dmr_settings_identical'] and not c['empirical_fdr'] for c in checks)
            assert source_inventory(source / 'epykit') == current
            assert contract(datasets, real) == frozen
            record['status'] = 'complete'
            write_json(OUTPUT / 'run.json', record)
            print(f'Complete: {OUTPUT}', flush=True)
        except BaseException as exc:
            record.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            write_json(OUTPUT / 'run.json', record)
            raise


if __name__ == '__main__':
    main()
