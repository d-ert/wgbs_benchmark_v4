"""Verify that methylKit memory cleanup preserves the native tables."""
import csv
import json
import os
from pathlib import Path
import subprocess

ROOT = Path('/scratch/wgbs_benchmark_v4')
HERE = Path(__file__).resolve().parent
fixture = HERE/'smoke_input'
fixture.mkdir(exist_ok=True)
with (ROOT/'data/real/samples.tsv').open() as handle:
    reader = csv.DictReader(handle,delimiter='\t')
    fields = reader.fieldnames
    rows = list(reader)
for row in rows:
    source = Path(row['path'])
    dest = fixture/source.name
    with source.open() as read,dest.open('w') as write:
        for index,line in enumerate(read):
            if index==5000:break
            write.write(line)
    row['path'] = str(dest)
with (fixture/'samples.tsv').open('w') as handle:
    writer = csv.DictWriter(handle,fields,delimiter='\t')
    writer.writeheader();writer.writerows(rows)
(fixture/'manifest.json').write_text(json.dumps({'assembly':'hg19'}))
env = dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
for name,runner in [('original','run_r_tool_original.R'),('cleanup','run_r_tool_uncapped.R')]:
    out = HERE/('smoke_'+name)
    out.mkdir(exist_ok=True)
    with (out/'run.log').open('w') as log:
        subprocess.run(['Rscript',str(HERE/runner),'methylKit',str(fixture),str(out)],
                       env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for name in ['dml.tsv','dmr.tsv','dmr_candidates.tsv']:
    assert (HERE/'smoke_original'/name).read_bytes()==(HERE/'smoke_cleanup'/name).read_bytes(),name
(HERE/'smoke_verified.json').write_text(json.dumps({'status':'verified',
    'scope':'first 5000 coverage rows per real sample; same 12 sample labels',
    'byte_identical_tables':['dml.tsv','dmr.tsv','dmr_candidates.tsv']},indent=2)+'\n')
print('MethylKit cleanup preserves all three native tables byte-for-byte.',flush=True)
