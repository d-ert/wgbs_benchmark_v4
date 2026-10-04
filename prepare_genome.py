#!/usr/bin/env python3
"""Prepare matched genome-wide signal/null inputs from the 206 BLUEPRINT sources."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'vendor310'), str(ROOT/'src'), str(ROOT)]
from wgbs_v3.provenance import digest, write_json
from simulation.calibrate_parallel import calibrate_genome
from simulation.reference import validate as validate_reference
from simulation.simulate import Scenario, simulate

SOURCE = Path('/scratch/wgbs_benchmark_v3/data/blueprint_206')
MATRIX = Path('/scratch/epykit_wgbs_benchmark/sources/BLUEPRINT/macrophage_counts.tsv.gz')
FASTA = Path('/scratch/mimosa_runs/mimosa_run_e9bcf13e0d08_20260814_152006_ca475c26/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--include-sex-chromosomes', action='store_true')
    args=ap.parse_args()
    chroms=[f'chr{i}' for i in range(1,23)] + (['chrX','chrY'] if args.include_sex_chromosomes else [])
    scope='nuclear' if args.include_sex_chromosomes else 'autosomes'
    work=ROOT/'data'/f'genome_{scope}_preparation'
    work.mkdir(parents=True,exist_ok=True)
    status=work/'status.json'
    def stage(name, **extra):
        write_json(status,dict(stage=name,chromosomes=chroms,**extra))
        print(name,flush=True)
    try:
        frozen_sources={p.name:digest(p) for p in (ROOT/'simulation').glob('*.py')}
        stage('checking_sources')
        verified=json.loads((SOURCE/'verified_sources.json').read_text())
        samples=json.loads((ROOT/'archive/v3/context/source_inventory.json').read_text())
        if isinstance(samples,dict):
            samples=samples['samples']
        assert len(samples)==206
        # Verify source content, not just presence of an old verification file.
        proof=work/'source_verification.json'
        if not proof.exists():
            for i,row in enumerate(verified['verified'],1):
                for kind in ('coverage','call'):
                    path=SOURCE/f"{row['sample_id']}.{kind}.bw"
                    if digest(path)!=row[f'{kind}_sha256']:
                        raise ValueError(f'Source changed: {path}')
                if i%20==0: print(f'Verified {i}/206 source pairs',flush=True)
            write_json(proof,dict(passed=True,source_manifest_sha256=digest(SOURCE/'verified_sources.json')))
        calibration=work/'calibration'
        if not (calibration/'manifest.json').exists():
            stage('calibrating')
            calibrate_genome(MATRIX,calibration,samples,SOURCE,chroms,FASTA)
        manifest=json.loads((calibration/'manifest.json').read_text())
        if manifest['scope']['chromosomes']!=chroms: raise ValueError('Calibration scope mismatch')
        if not (calibration/'reference_validation.json').exists():
            stage('validating_reference')
            validate_reference(calibration,FASTA)
        for mode,label in [('dmr','signal'),('null','null')]:
            out=ROOT/'data'/f'genome_{scope}_{label}'
            if not (out/'manifest.json').exists():
                stage(f'simulating_{label}')
                simulate(calibration,out,Scenario(seed=20261003,n_per_group=5,coverage=20,
                         delta=.2,mode=mode),chromosomes=chroms)
            if not (out/'validation.json').exists():
                stage(f'validating_{label}')
                from simulation.validate_streaming import validate_dataset
                validate_dataset(out,chroms)
        if frozen_sources!={p.name:digest(p) for p in (ROOT/'simulation').glob('*.py')}:
            raise ValueError('Simulation code changed during preparation')
        shutil.copytree(ROOT/'simulation',work/'source_snapshot',dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__'))
        write_json(work/'source_hashes.json',frozen_sources)
        stage('complete',datasets=[f'genome_{scope}_signal',f'genome_{scope}_null'])
    except BaseException as exc:
        stage('failed',error=f'{type(exc).__name__}: {exc}')
        raise

if __name__=='__main__': main()
