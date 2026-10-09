import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def build_samples(raw):
    raw=Path(raw)
    sheet=pd.read_csv(raw/'covs/samplesheet.csv')
    mapping={'WT':'control','Het_AKAP11_KO':'treatment'}
    if len(sheet)!=6 or set(sheet.group)!=set(mapping) or any(sheet.group.value_counts()!=3):
        raise ValueError('GSE263850 requires three WT and three heterozygous KO samples')
    if sheet.sample_id.duplicated().any():raise ValueError('Duplicate sample ID')
    sheet['original_group']=sheet.group
    sheet['group']=sheet.group.map(mapping)
    sheet['path']=[str((raw/'covs'/Path(p).name).resolve()) for p in sheet.path]
    if not all(Path(p).is_file() for p in sheet.path):raise ValueError('Missing coverage file')
    return sheet


def load_paper(path):
    p=pd.read_csv(path,sep=';',decimal=',')
    required={'chr','start','end','length','nCG','diff.meth_mean'}
    if not required<=set(p):raise ValueError('Missing published DMR columns')
    if p.duplicated(['chr','start','end']).any():raise ValueError('Duplicate published DMR')
    if not ((p.start>=1)&(p.end>=p.start)&(p.end-p.start+1==p.length)).all():
        raise ValueError('Published coordinates are not 1-based inclusive')
    if p['diff.meth_mean'].isna().any() or (p['diff.meth_mean']==0).any():
        raise ValueError('Missing/zero published methylation direction')
    result=pd.DataFrame(dict(chrom=p.chr,start=p.start.astype(int),end=p.end.astype(int),
        direction=np.where(p['diff.meth_mean']>0,'hyper','hypo'),meth_diff=p['diff.meth_mean'],
        n_cpgs=p.nCG.astype(int),gene=p.get('Gene.Name','')))
    return result


def validate_cov(path,chunksize=500000):
    rows=0;last_chrom=None;last_pos=0;finished=set();chrom_counts={};minimum=None
    names=['chrom','start','end','percent','m','u']
    for d in pd.read_csv(path,sep='\t',header=None,names=names,chunksize=chunksize,
                         dtype={'chrom':str,'start':np.int64,'end':np.int64,'m':np.int64,'u':np.int64}):
        n=d.m+d.u
        if d.isna().any().any() or not ((d.start==d.end)&(d.start>=1)&(d.m>=0)&(d.u>=0)&(n>=5)).all():
            raise ValueError(f'Invalid integer counts/coordinates or coverage below five in {path}')
        if not np.all(np.abs(d.percent.to_numpy()-100*d.m.to_numpy()/n.to_numpy())<=.011):
            raise ValueError(f'Count/percentage mismatch in {path}')
        for chrom,g in d.groupby('chrom',sort=False):
            pos=g.start.to_numpy()
            if chrom!=last_chrom:
                if chrom in finished:raise ValueError('Chromosome reappears in coverage file')
                if last_chrom is not None:finished.add(last_chrom)
                last_pos=0
            if pos[0]<=last_pos or np.any(np.diff(pos)<=0):raise ValueError('Unsorted/duplicate CpG')
            last_chrom,last_pos=chrom,int(pos[-1])
            chrom_counts[chrom]=chrom_counts.get(chrom,0)+len(g)
        rows+=len(d);minimum=int(n.min()) if minimum is None else min(minimum,int(n.min()))
    if not rows:raise ValueError('Empty coverage file')
    return dict(rows=rows,min_coverage=minimum,chromosomes=chrom_counts,sha256=digest(path))


def check_bed_cov(bed,cov,limit=10000):
    checked=0
    with gzip.open(bed,'rt') as a,Path(cov).open() as b:
        for raw,converted in zip(a,b):
            x,y=raw.split(),converted.split()
            if len(x)!=12 or len(y)!=6:raise ValueError('Malformed processed BED/COV')
            expected=(x[0],int(x[1])+1,int(x[9]),int(x[10])-int(x[9]))
            observed=(y[0],int(y[1]),int(y[4]),int(y[5]))
            if expected!=observed or int(x[2])!=int(x[1])+1 or y[1]!=y[2]:
                raise ValueError('BED/COV offset or combined counts mismatch')
            if int(x[9])!=int(x[3])+int(x[6]) or int(x[10])!=int(x[4])+int(x[7]):
                raise ValueError('Combined BED counts differ from strand counts')
            checked+=1
            if checked>=limit:break
    if not checked:raise ValueError('Empty BED/COV comparison')
    return checked


def prepare(raw,out,max_sites=None):
    raw,out=Path(raw).resolve(),Path(out).resolve()
    if out.exists():raise ValueError(f'Preserving existing prepared dataset: {out}')
    sheet=build_samples(raw)
    paper=load_paper(raw/'Paper data/supp table 5_DMR List.csv')
    if len(paper)!=813:raise ValueError('Expected 813 published DMRs')
    out.mkdir(parents=True)
    audit=[]
    for index,row in sheet.iterrows():
        path=Path(row.path)
        bed=raw/(path.stem.removesuffix('.cov')+'.bed.gz')
        checked=check_bed_cov(bed,path)
        if max_sites:
            dest=out/path.name
            with path.open() as src,dest.open('w') as target:
                for i,line in enumerate(src):
                    if i==max_sites:break
                    target.write(line)
            path=dest;sheet.loc[index,'path']=str(path)
        stats=validate_cov(path)
        audit.append(dict(sample_id=row.sample_id,path=str(path),source=str(Path(row.path)),
                          bed_cov_rows_checked=checked,**stats))
        print(f'Validated {row.sample_id}: {stats["rows"]:,} CpGs',flush=True)
    sheet.to_csv(out/'samples.tsv',sep='\t',index=False)
    paper.to_csv(out/'paper_dmrs.tsv',sep='\t',index=False)
    ref=raw/'covs/refseq/refGene.txt.gz'
    genes=pd.read_csv(ref,sep='\t',header=None,usecols=[2,3,4,5,12])
    genes.columns=['chrom','strand','start','end','symbol']
    genes=genes.groupby(['symbol','chrom','strand'],as_index=False).agg(start=('start','min'),end=('end','max'))
    genes['start']+=1
    genes['gene_id']=genes.symbol+':'+genes.chrom+':'+genes.strand
    genes.to_csv(out/'genes.tsv',sep='\t',index=False)
    manifest=dict(kind='real',study='GSE263850',assembly='hg38',strand_collapsed=True,
        coordinate_system='1-based inclusive',min_coverage=5,
        contrast='Het_AKAP11_KO minus WT',paired=False,common_parental_line='SBP009',
        source_dir=str(raw),max_sites=max_sites,samples=audit,paper_regions=813,
        paper_source_sha256=digest(raw/'Paper data/supp table 5_DMR List.csv'),
        refgene_source_sha256=digest(ref))
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return out
