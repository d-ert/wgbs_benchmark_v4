import pandas as pd
from wgbs_real.regions import clump_sites, paper_metrics


def test_raw_p_region_geometry():
    sites=pd.DataFrame(dict(chrom=['chr1']*6,pos=[10,30,50,70,90,250],
        pvalue=[1e-6,1e-6,.2,1e-6,1e-6,1e-6],meth_diff=[.1]*6))
    regions=clump_sites(sites)
    assert len(regions)==1
    assert regions.iloc[0].start==10 and regions.iloc[0].end==90
    assert regions.iloc[0].n_cpgs==5 and regions.iloc[0].direction==1


def test_missing_p_and_wrong_direction_do_not_match():
    calls=pd.DataFrame(dict(chrom=['chr1']*2,start=[10,10],end=[90,90],direction=[1,-1]))
    paper=pd.DataFrame(dict(chrom=['chr1'],start=[10],end=[90],direction=['hyper']))
    result,pairs=paper_metrics(calls,paper)
    assert result[0]['one_to_one_matches']==1
    assert result[0]['one_to_one_precision']==.5
    assert result[0]['one_to_one_recall']==1
    assert result[0]['discordant_pairs']==1


def test_empty_calls_have_undefined_precision():
    calls=pd.DataFrame(columns=['chrom','start','end','direction'])
    paper=pd.DataFrame(dict(chrom=['chr1'],start=[10],end=[90],direction=['hyper']))
    result,_=paper_metrics(calls,paper)
    assert result[0]['one_to_one_precision'] is None
    assert result[0]['one_to_one_recall']==0
