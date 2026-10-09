import json
from pathlib import Path

import pandas as pd
import pytest

from wgbs_real.dataset import build_samples, load_paper, validate_cov, check_bed_cov


def test_stale_paths_and_contrast(tmp_path):
    root=tmp_path/'covs';root.mkdir()
    rows=[]
    for i in range(6):
        name=f'GSM{i}.cov';(root/name).write_text('chr1\t11\t11\t40\t2\t3\n')
        rows.append(dict(sample_id=f's{i}',path='/obsolete/'+name,
                         group='Het_AKAP11_KO' if i<3 else 'WT'))
    pd.DataFrame(rows).to_csv(root/'samplesheet.csv',index=False)
    result=build_samples(tmp_path)
    assert result.group.tolist()==['treatment']*3+['control']*3
    assert all(Path(p).is_file() for p in result.path)


def test_wrong_group_balance_refused(tmp_path):
    (tmp_path/'covs').mkdir()
    pd.DataFrame([dict(sample_id='x',path='x.cov',group='WT')]).to_csv(
        tmp_path/'covs/samplesheet.csv',index=False)
    with pytest.raises(ValueError):build_samples(tmp_path)


def test_locale_and_paper_direction(tmp_path):
    p=tmp_path/'paper.csv'
    p.write_text('chr;start;end;length;nCG;diff.meth_mean;Gene.Name\n'
                 'chr1;11;64;54;4;0,3;X\nchr2;20;81;62;5;-0,25;Y\n')
    result=load_paper(p)
    assert result.direction.tolist()==['hyper','hypo']
    assert result.meth_diff.tolist()==[.3,-.25]


@pytest.mark.parametrize('body',[
    'chr1\t11\t11\t40\t2\t3\nchr1\t11\t11\t40\t2\t3\n',
    'chr1\t11\t12\t40\t2\t3\n',
    'chr1\t11\t11\t50\t2\t2\n',
    'chr1\t11\t11\t40\t-1\t6\n'])
def test_malformed_cov_refused(tmp_path,body):
    p=tmp_path/'sample.cov';p.write_text(body)
    with pytest.raises(ValueError):validate_cov(p)


def test_cov_boundary_and_combined_bed_counts(tmp_path):
    p=tmp_path/'sample.cov';p.write_text('chr1\t11\t11\t40\t2\t3\nchr1\t21\t21\t60\t3\t2\n')
    assert validate_cov(p,chunksize=1)['rows']==2
    import gzip
    bed=tmp_path/'sample.bed.gz'
    with gzip.open(bed,'wt') as f:f.write('chr1\t10\t11\t1\t3\t33.33\t1\t2\t50\t2\t5\t40\n')
    assert check_bed_cov(bed,p,limit=1)==1
    p.write_text('chr1\t10\t10\t40\t2\t3\n')
    with pytest.raises(ValueError):check_bed_cov(bed,p,limit=1)
