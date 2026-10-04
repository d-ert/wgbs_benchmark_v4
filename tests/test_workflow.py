"""Exercise orchestration and isolation without running expensive statistical callers."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('workflow', ROOT / 'benchmark.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.setattr(w, 'ROOT', tmp_path)
    (tmp_path / 'results').mkdir()
    for name in w.SIMS:
        dataset = tmp_path / 'data' / name
        dataset.mkdir(parents=True)
        (dataset / 'manifest.json').write_text('{}')
        count = dataset / 'counts.tsv'
        count.write_text('test')
        pd.DataFrame([dict(sample_id='x', group='control', path=str(count))]).to_csv(dataset / 'samples.tsv', sep='\t', index=False)
    calls = []
    monkeypatch.setattr(w, 'preflight', lambda p: None)
    monkeypatch.setattr(w, 'prepare_real', lambda *a: None)
    monkeypatch.setattr(w, 'contract', lambda *a: {'frozen': 'same'})
    monkeypatch.setattr(w, 'source_inventory', lambda *a: {'epykit.py': 'v1'})
    monkeypatch.setattr(w, 'snapshot', lambda *a: None)
    monkeypatch.setattr(w.subprocess, 'run', lambda *a, **k: None)

    def simulated(datasets, out, tools, ep_source):
        calls.append(('sim', tools))
        out.mkdir(parents=True)
        pd.DataFrame([dict(dataset=d.name, tool=t, status='ok', score=1 if t=='epykit' else 2,
                           wall_seconds=10, peak_rss_bytes=100)
                      for d in datasets for t in tools]).to_csv(out / 'simulation_metrics.csv', index=False)

    def real(dataset, out, tools, ep_source, extras):
        calls.append(('real', tools))
        records = {}
        for tool in tools:
            (out / tool).mkdir(parents=True)
            (out / tool / 'dmr.tsv').write_text('chrom\tstart\tend\nchr1\t1\t10\n')
            records[tool] = dict(status='ok', monitoring=dict(wall_seconds=10, peak_process_tree_rss_bytes=100))
        w.write_json(out / 'run_manifest.json', records)
    monkeypatch.setattr(w, 'simulated', simulated)
    monkeypatch.setattr(w, 'run_real', real)
    return tmp_path, calls


def args(command, name, **kwargs):
    return argparse.Namespace(command=command, name=name, epykit_source=Path('/unused'),
                              baseline='baseline', smoke=False, dry_run=False, **kwargs)


def test_baseline_then_only_epykit_and_reuse(workspace):
    root, calls = workspace
    w.execute(args('baseline', 'baseline'))
    baseline_bytes = (root / 'results/baseline/simulated/simulation_metrics.csv').read_bytes()
    w.execute(args('epykit', 'candidate'))
    assert calls == [('sim', w.TOOLS), ('real', w.TOOLS), ('sim', ('epykit',)), ('real', ('epykit',))]
    comparison = pd.read_csv(root / 'results/candidate/comparison/simulation_metrics.csv')
    assert len(comparison) == 12
    assert set(comparison[comparison.tool == 'epykit'].source_run) == {'candidate'}
    assert comparison[comparison.tool != 'epykit'].reused.all()
    assert (root / 'results/baseline/simulated/simulation_metrics.csv').read_bytes() == baseline_bytes


def test_changed_inputs_block_rerun(workspace, monkeypatch):
    root, calls = workspace
    w.execute(args('baseline', 'baseline'))
    monkeypatch.setattr(w, 'contract', lambda *a: {'frozen': 'changed'})
    with pytest.raises(ValueError, match='changed'):
        w.execute(args('epykit', 'candidate'))
    assert len(calls) == 2
    assert not (root / 'results/candidate').exists()


def test_failed_baseline_cannot_be_reused(workspace, monkeypatch):
    root, calls = workspace
    def fail(*args):
        raise RuntimeError('caller failure')
    monkeypatch.setattr(w, 'run_real', fail)
    with pytest.raises(RuntimeError):
        w.execute(args('baseline', 'baseline'))
    assert w.read_json(root / 'results/baseline/run.json')['status'] == 'failed'
    with pytest.raises(ValueError, match='not complete'):
        w.execute(args('epykit', 'candidate'))


def test_refuse_overwrite(workspace):
    w.execute(args('baseline', 'baseline'))
    with pytest.raises(ValueError, match='already exists'):
        w.execute(args('baseline', 'baseline'))


def test_source_edit_during_run_rejects_completion(workspace, monkeypatch):
    root, calls = workspace
    inventories = iter([{'code': 'before'}, {'code': 'after'}])
    monkeypatch.setattr(w, 'source_inventory', lambda *a: next(inventories))
    with pytest.raises(ValueError, match='changed during'):
        w.execute(args('baseline', 'baseline'))
    assert w.read_json(root / 'results/baseline/run.json')['status'] == 'failed'
