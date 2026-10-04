#!/usr/bin/env python3
"""Frozen baseline, epykit-only reruns, and comparisons against cached competitors."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'vendor310'), str(ROOT / 'src'), str(ROOT)]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'data/mplcache'))
import pandas as pd
from analysis.run_publication_workflow import environment, prepare_real, real_samples, run_real, simulated
from wgbs_v3.benchmark import TOOLS
from wgbs_v3.provenance import dataset_inventory, snapshot, source_inventory, write_json

SIMS = ('chr1_3_full_baseline', 'calibration_v5_full_null_870001')
GENOME = ('genome_autosomes_signal', 'genome_autosomes_null')
PILOTS = ('chr1_3_206_pilot', 'chr1_3_206_null_pilot')


def read_json(path):
    return json.loads(Path(path).read_text())


def run_path(name):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
        raise ValueError('Run names must contain only letters, numbers, _, . or -')
    return ROOT / 'results' / name


def harness_inventory():
    files = [ROOT / 'benchmark.py', ROOT / 'run_epykit.py', ROOT / 'run_r_tool.R']
    files += list((ROOT / 'analysis').glob('*.*')) + list((ROOT / 'src/wgbs_v3').glob('*.py'))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files) if p.is_file()}


def contract(datasets, real_input):
    return dict(simulated={p.name: dataset_inventory(p) for p in datasets},
                real=dataset_inventory(real_input), harness=harness_inventory())


def check_complete(path):
    record = read_json(path / 'run.json')
    if record['status'] != 'complete':
        raise ValueError(f'Run is not complete: {path}')
    return record


def validate_baseline(path):
    record = check_complete(path)
    if record['kind'] != 'baseline' or record['tools'] != list(TOOLS):
        raise ValueError('Expected a complete six-tool baseline')
    return record


def preflight(ep_source):
    environment(ep_source)
    if not (ep_source / 'epykit/__init__.py').is_file():
        raise ValueError(f'No epykit package at {ep_source}')
    if sys.version_info[:2] != (3, 10):
        raise ValueError('Copied vendor dependencies require Python 3.10')
    env = environment(ep_source)
    subprocess.run([sys.executable, '-c',
        'import epykit, polars, pandas, numpy, scipy, pyarrow, psutil, sklearn, networkx'],
        env=env, check=True)
    subprocess.run(['Rscript', '-e',
        'p <- c("data.table","bsseq","DSS","methylKit","dmrseq","DMRcate",'
        '"jsonlite","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db"); '
        'stopifnot(all(vapply(p, requireNamespace, logical(1), quietly=TRUE))); '
        'cat("R dependencies OK\\n")'], env=env, check=True)


def compare(baseline, candidate, out):
    base = validate_baseline(baseline)
    new = check_complete(candidate)
    if new['kind'] != 'epykit' or Path(new['baseline']) != baseline:
        raise ValueError('Candidate does not belong to this baseline')
    if new['contract'] != base['contract']:
        raise ValueError('Input or harness contracts differ')
    out.mkdir(parents=True, exist_ok=True)
    before = pd.read_csv(baseline / 'simulated/simulation_metrics.csv')
    after = pd.read_csv(candidate / 'simulated/simulation_metrics.csv')
    if set(after.tool) != {'epykit'} or set(after.dataset) != set(before.dataset):
        raise ValueError('Candidate simulation coverage differs')
    before['source_run'] = baseline.name
    after['source_run'] = candidate.name
    before['reused'] = True
    after['reused'] = False
    combined = pd.concat([before[before.tool != 'epykit'], after], ignore_index=True)
    combined.to_csv(out / 'simulation_metrics.csv', index=False)
    old = before[before.tool == 'epykit'].set_index('dataset')
    new_metrics = after.set_index('dataset')
    changes = []
    for dataset in old.index:
        for col in old.select_dtypes(include='number').columns.intersection(new_metrics.select_dtypes(include='number').columns):
            a, b = old.at[dataset, col], new_metrics.at[dataset, col]
            changes.append(dict(dataset=dataset, metric=col, baseline=a, candidate=b, delta=b-a))
    pd.DataFrame(changes).to_csv(out / 'epykit_changes.csv', index=False)
    real = out / 'real'
    real.mkdir(exist_ok=True)
    sources = {}
    performance = []
    for tool in TOOLS:
        origin = candidate if tool == 'epykit' else baseline
        dest = real / tool
        dest.mkdir(exist_ok=True)
        shutil.copy2(origin / 'real' / tool / 'dmr.tsv', dest / 'dmr.tsv')
        sources[tool] = dict(source_run=origin.name, reused=tool != 'epykit')
        record = read_json(origin / 'real/run_manifest.json')[tool]
        performance.append(dict(tool=tool, **sources[tool], status=record['status'],
            wall_seconds=record['monitoring']['wall_seconds'],
            peak_rss_bytes=record['monitoring']['peak_process_tree_rss_bytes']))
    pd.DataFrame(performance).to_csv(out / 'real_performance.csv', index=False)
    subprocess.run(['Rscript', str(ROOT / 'analysis/annotate_real_dmrs.R'),
                    str(real), str(real / 'concordance'), ','.join(TOOLS)], check=True)
    write_json(out / 'sources.json', dict(baseline=str(baseline), candidate=str(candidate),
        real=sources, note='Competitor calls and timings reused from baseline; real concordance is descriptive.'))


def execute(args):
    ep_source = args.epykit_source.resolve()
    preflight(ep_source)
    out = run_path(args.name)
    if out.exists():
        raise ValueError(f'Output already exists: {out}. Use a new --name; previous runs are preserved.')
    base = None
    if args.command == 'epykit':
        base_path = run_path(args.baseline)
        base = validate_baseline(base_path)
        datasets = [Path(p) for p in base['datasets']]
        real_input = Path(base['real_input'])
        tools = ('epykit',)
    else:
        if args.smoke and getattr(args, 'whole_genome', False):
            raise ValueError('--smoke and --whole-genome are mutually exclusive')
        selection = GENOME if getattr(args, 'whole_genome', False) else SIMS
        datasets = [ROOT / 'data' / name for name in (PILOTS if args.smoke else selection)]
        real_input = ROOT / 'data' / ('real_smoke' if args.smoke else 'real')
        tools = TOOLS
    for dataset in datasets:
        if not (dataset / 'manifest.json').is_file():
            raise FileNotFoundError(dataset)
        if dataset.name in GENOME:
            validation = read_json(dataset / 'validation.json')
            actual_hash = hashlib.sha256((dataset / 'manifest.json').read_bytes()).hexdigest()
            if not validation.get('passed') or validation.get('dataset_manifest_sha256') != actual_hash:
                raise ValueError(f'Whole-genome input has not passed validation: {dataset}')
        sheet = pd.read_csv(dataset / 'samples.tsv', sep='\t')
        if not all(Path(p).is_file() for p in sheet.path):
            raise ValueError(f'Missing count files in {dataset}')
    plan = dict(kind=args.command, tools=list(tools), datasets=list(map(str, datasets)),
                real_input=str(real_input), output=str(out), epykit_source=str(ep_source))
    if args.dry_run:
        if base is None:
            real_samples(ROOT / "data/real_cov")
        print(json.dumps(plan, indent=2))
        return
    # Lock the entire workflow, including real data and normalization.
    with (ROOT / 'results/.workflow.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if base is None:
            prepare_real(ROOT / 'data/real_cov', real_input, 10000 if args.smoke else None)
        print('Hashing frozen inputs and benchmark code...', flush=True)
        frozen = contract(datasets, real_input)
        if base is not None and frozen != base['contract']:
            raise ValueError('Baseline input or benchmark code changed. Restore it or create a new baseline.')
        ep_hashes = source_inventory(ep_source / 'epykit')
        out.mkdir()
        record = dict(**plan, name=args.name, status='running', contract=frozen,
                      epykit_hashes=ep_hashes, baseline=str(base_path) if base else None,
                      smoke=base['smoke'] if base else args.smoke)
        write_json(out / 'run.json', record)
        try:
            # Real data also receives an actual source snapshot and input hashes.
            snapshot(ROOT, ep_source, real_input, out)
            print('Running simulations...', flush=True)
            simulated(datasets, out / 'simulated', tools, ep_source)
            print('Running real data...', flush=True)
            run_real(real_input, out / 'real', tools, ep_source, [])
            if source_inventory(ep_source / 'epykit') != ep_hashes:
                raise ValueError('epykit source changed during the run')
            if contract(datasets, real_input) != frozen:
                raise ValueError('Inputs or benchmark code changed during the run')
            record['status'] = 'complete'
            write_json(out / 'run.json', record)
        except BaseException as exc:
            record.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            write_json(out / 'run.json', record)
            raise
        if base is not None:
            compare(base_path, out, out / 'comparison')
        print(f'Complete: {out}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    commands = ap.add_subparsers(dest='command', required=True)
    for name in ('baseline', 'epykit'):
        sub = commands.add_parser(name)
        sub.add_argument('--name', default='baseline' if name == 'baseline' else None, required=name == 'epykit')
        sub.add_argument('--epykit-source', type=Path, default=Path('/scratch/epykit3/src'))
        sub.add_argument('--dry-run', action='store_true', help='Check dependencies and show the plan; run no callers')
        if name == 'baseline':
            sub.add_argument('--smoke', action='store_true', help='Use simulation pilots and 10,000 real rows per sample')
            sub.add_argument('--whole-genome', action='store_true', help='Use prepared chr1–22 signal/null simulations')
        else:
            sub.add_argument('--baseline', default='baseline')
    sub = commands.add_parser('compare', help='Rebuild comparison without rerunning callers')
    sub.add_argument('--baseline', default='baseline')
    sub.add_argument('--run', required=True)
    args = ap.parse_args()
    os.chdir(ROOT)
    if args.command == 'compare':
        candidate = run_path(args.run)
        compare(run_path(args.baseline), candidate, candidate / 'comparison')
    else:
        execute(args)


if __name__ == '__main__':
    main()
