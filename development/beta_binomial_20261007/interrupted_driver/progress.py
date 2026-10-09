"""Read-only progress from atomic chromosome outputs and completed manifests.

This reports output progress; a live exec handle must independently prove
that a still-running process exists.
"""
from pathlib import Path
import datetime
import json
import sys

ROOT = Path('/scratch/wgbs_benchmark_v4')
RUN = ROOT/'results/beta_binomial_20261007'
# Original result caches were reused by earlier diagnostic probes. The
# accepted region-only run retains complete, byte-identical original DMLs.
BASE = ROOT/'results/default_search_by_20261007'
sys.path[:0] = [str(ROOT/'vendor310')]
import pyarrow.parquet as pq

record = json.loads((RUN/'run.json').read_text())
rows = []
for name,relative in [('signal','simulated/genome_autosomes_signal/epykit'),
                      ('null','simulated/genome_autosomes_null/epykit'),('GSE64177','real/epykit')]:
    target = RUN/relative
    if not target.exists():
        rows.append(dict(dataset=name,stage='not_started'))
        continue
    old_manifest = BASE/relative/'store/.cache/dmc/lr/.epykit_dmc_manifest.json'
    old = json.loads(old_manifest.read_text())
    current = target/'store/.cache/dmc/beta_binomial'
    fitted = list(current.glob('chrom=*.parquet'))
    expected_chromosomes = len(old['chroms'])
    expected_rows = sum(entry.get('n_rows',entry.get('n_sites',0)) for entry in old['chroms'])
    fitted_rows = sum(pq.read_metadata(path).num_rows for path in fitted)
    rows.append(dict(dataset=name,converted_samples=len(list(target.rglob('.epykit_raw_manifest.json'))),
        fitted_chromosomes=len(fitted),expected_chromosomes=expected_chromosomes,
        fitted_cpgs=fitted_rows,expected_cpgs=expected_rows,
        fitted_fraction=fitted_rows/expected_rows if expected_rows else None,
        exported_dml=(target/'dml.tsv').exists(),exported_dmr=(target/'dmr.tsv').exists(),
        scored=(target/'score.json').exists()))
print(json.dumps(dict(checked_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    run_status=record['status'],run_stage=record['stage'],datasets=rows),indent=2))
