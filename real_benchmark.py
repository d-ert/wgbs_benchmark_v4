#!/usr/bin/env python3
"""Additional real-study pipeline; recorded baseline harness remains frozen."""
import argparse
import json
import re
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'vendor310'),str(ROOT/'src'),str(ROOT)]
from wgbs_real.pipeline import run,TOOLS
from wgbs_real.dataset import prepare
from wgbs_real.comparison import compare

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare');p.add_argument('--source',type=Path,default=ROOT/'data/GSE263850');p.add_argument('--out',type=Path,required=True)
    p=sub.add_parser('run');p.add_argument('--source',type=Path,default=ROOT/'data/GSE263850')
    p.add_argument('--name',required=True);p.add_argument('--profile',choices=['paper_aligned','native'],default='paper_aligned')
    p.add_argument('--tools',default=','.join(TOOLS));p.add_argument('--smoke-sites',type=int)
    p.add_argument('--epykit-source',type=Path,default=Path('/scratch/epykit-calibration-fresh_20261007/src'))
    for name in ['status','compare']:
        p=sub.add_parser(name);p.add_argument('--name',required=True)
    args=parser.parse_args()
    if args.command=='prepare':prepare(args.source,args.out);return
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',args.name):raise ValueError('Invalid run name')
    out=ROOT/'results'/args.name
    if args.command=='run':run(args.source,out,args.epykit_source,args.profile,args.tools.split(','),args.smoke_sites)
    elif args.command=='compare':compare(out)
    else:
        print((out/'run.json').read_text())
        path=out/'real/run_manifest.json'
        if path.exists():print(path.read_text())

if __name__=='__main__':main()
