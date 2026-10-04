"""Validate large simulations with chromosome truth and streamed sample counts."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from wgbs_v3.provenance import digest, write_json


def validate_dataset(dataset, chromosomes):
    dataset=Path(dataset)
    meta=json.loads((dataset/'manifest.json').read_text())
    assert meta['chromosomes']==chromosomes
    regions=pd.read_csv(dataset/'dmr_truth.tsv',sep='\t')
    assert not regions.region_id.duplicated().any()
    positions={}; observed={}; cached={}; npositive=0; found_regions=set()
    for chrom in chromosomes:
        truth=pq.read_table(dataset/'truth_parts'/f'{chrom}.parquet').to_pandas()
        assert len(truth)>0 and set(truth.chrom)=={chrom}
        pos=truth.pos.to_numpy()
        assert np.all(np.diff(pos)>0) and np.all(pos>=1)
        assert np.allclose(truth.treatment_mu-truth.control_mu,truth.signed_delta,atol=1e-6)
        assert truth.control_mu.between(.01-1e-6,.99+1e-6).all()
        assert ((truth.treatment_mu>0)&(truth.treatment_mu<1)).all()
        assert np.array_equal(truth.is_dml,truth.region_id!=0)
        assert (truth.loc[truth.region_id==0,'signed_delta']==0).all()
        assert (truth.loc[truth.region_id!=0,'signed_delta']!=0).all()
        for rid,members in truth[truth.region_id>0].groupby('region_id'):
            row=regions.set_index('region_id').loc[rid]
            assert row.chrom==chrom and len(members)==row.n_cpg
            assert (members.pos.min(),members.pos.max())==(row.start,row.end)
            assert len(members)>=5 and np.diff(members.pos).max(initial=0)<=500
            assert np.all(np.sign(members.signed_delta)==row.direction)
            assert np.isclose(members.signed_delta.mean(),row.mean_delta,atol=1e-6)
            assert np.isclose(abs(row.mean_delta),meta['scenario']['delta'],atol=1e-6)
            found_regions.add(rid)
        eligibility=pq.read_table(dataset/'eligibility_parts'/f'{chrom}.parquet').to_pandas()
        assert np.array_equal(eligibility.pos,pos) and set(eligibility.chrom)=={chrom}
        positions[chrom]=pos
        cached[chrom]=eligibility.eligible.to_numpy()
        observed[chrom]=np.zeros(len(pos),dtype=np.uint16)
        npositive+=int(truth.is_dml.sum())
    assert found_regions==set(regions.region_id)
    nsites=sum(map(len,positions.values()))
    assert nsites==meta['n_sites'] and npositive==meta['n_positive_cpg']
    assert len(regions)==meta['n_regions']
    if meta['scenario']['mode']=='null': assert npositive==0 and len(regions)==0
    sheet=pd.read_csv(dataset/'samples.tsv',sep='\t')
    n=meta['scenario']['n_per_group']
    assert sheet.groupby('group').size().to_dict()=={'control':n,'treatment':n}
    assert not sheet.sample_id.duplicated().any() and not sheet.path.duplicated().any()
    per_sample={}; hashes={}; nrows=0
    for sample in sheet.itertuples(index=False):
        total_depth=0; last_chrom=-1;last_pos=0
        for chunk in pd.read_csv(sample.path,sep='\t',header=None,
                    names=['chrom','pos','end','percent','m','u'],chunksize=250000):
            values=chunk[['pos','end','m','u']].to_numpy(dtype=float)
            assert np.isfinite(values).all() and np.all(values==np.floor(values))
            assert (chunk.end==chunk.pos+1).all() and (chunk[['m','u']]>=0).all().all()
            depth=chunk.m+chunk.u
            assert (depth>0).all()
            assert np.allclose(chunk.percent,100*chunk.m/depth,atol=1e-5)
            for chrom,part in chunk.groupby('chrom',sort=False):
                rank=chromosomes.index(chrom)
                pos=part.pos.to_numpy()
                assert rank>=last_chrom and np.all(np.diff(pos)>0)
                if rank==last_chrom: assert pos[0]>last_pos
                last_chrom=rank;last_pos=pos[-1]
                indices=np.searchsorted(positions[chrom],pos)
                assert (indices<len(positions[chrom])).all()
                assert np.array_equal(positions[chrom][indices],pos)
                observed[chrom][indices]+=1
            total_depth+=int(depth.sum());nrows+=len(chunk)
        per_sample[sample.sample_id]=total_depth/nsites
        hashes[sample.sample_id]=digest(Path(sample.path))
        print(f'Validated {dataset.name}: {sample.sample_id}',flush=True)
    eligible=0
    for chrom in chromosomes:
        mask=observed[chrom]==len(sheet)
        assert np.array_equal(mask,cached[chrom])
        eligible+=int(mask.sum())
    assert eligible==meta['n_common_eligible']
    result=dict(passed=True,chromosomes=chromosomes,n_sites=nsites,n_regions=len(regions),
                n_positive_cpg=npositive,n_common_eligible=eligible,n_coverage_rows=nrows,
                per_sample_mean_depth=per_sample,count_file_sha256=hashes,
                dataset_manifest_sha256=digest(dataset/'manifest.json'))
    write_json(dataset/'validation.json',result)
    return result
