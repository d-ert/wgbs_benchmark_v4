"""Use the existing full-data scorers without modifying truth or matching."""
import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument('kind',choices=['simulation','real'])
    p.add_argument('dataset',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--paper',type=Path)
    a=p.parse_args()
    if a.kind=='simulation':
        from wgbs_v3.score import score
        r=score(a.dataset,a.output/'dml.tsv',a.output/'dmr.tsv',a.output/'score.json')
        print(json.dumps(dict(cpg=r.get('dml'),regions=r.get('dmr'))),flush=True)
    else:
        from wgbs_real.regions import paper_metrics
        paper=pd.read_csv(a.paper,sep='\t')
        calls=pd.read_csv(a.output/'dmr.tsv',sep='\t')
        if calls.duplicated(['chrom','start','end']).any():
            raise ValueError('Duplicate region outputs')
        metrics,pairs=paper_metrics(calls,paper)
        (a.output/'paper_score.json').write_text(json.dumps(dict(metrics=metrics,
            interpretation='Descriptive same-direction published-list concordance; not biological precision/recall'),indent=2)+'\n')
        pairs.to_csv(a.output/'paper_overlap_pairs.csv',index=False)
        print(json.dumps(metrics),flush=True)


if __name__=='__main__':
    main()
