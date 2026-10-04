"""Interrupted real input normalization must be safely repeatable."""
import hashlib
import json
from pathlib import Path

import pandas as pd
from analysis.run_publication_workflow import prepare_real


def test_resumes_partial_real_preparation(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    for donor in range(1, 7):
        for condition in ('NI', 'MTB'):
            (source / f'GSM1_DC{donor}_{condition}_5mC.cov').write_bytes(
                b'chr1\t10\t10\t33,33\t1\t2\nchr1\t20\t20\t50,00\t2\t2\n')
    target = tmp_path / 'real'
    counts = target / 'counts'
    counts.mkdir(parents=True)
    first = counts / 'GSM1_DC1_MTB_5mC.cov'
    first.write_bytes(b'incomplete previous run')
    (counts / 'GSM1_DC1_NI_5mC.cov.partial').write_bytes(b'incomplete temporary file')
    (target / 'samples.tsv.partial').write_text('incomplete sheet')
    prepare_real(source, target)
    manifest = json.loads((target / 'manifest.json').read_text())
    sheet = pd.read_csv(target / 'samples.tsv', sep='\t')
    assert len(sheet) == len(manifest['samples']) == 12
    assert first.read_bytes() == b'chr1\t10\t10\t33.33\t1\t2\nchr1\t20\t20\t50.00\t2\t2\n'
    assert not list(counts.glob('*.partial'))
    assert not (target / 'samples.tsv.partial').exists()
    for row in manifest['samples']:
        assert hashlib.sha256(Path(row['normalized']).read_bytes()).hexdigest() == row['normalized_sha256']
    assert prepare_real(source, target) == target
