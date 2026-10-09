"""One caller, including methylKit's explicitly labelled region adapter."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from wgbs_real.regions import clump_export

def main():
    tool,dataset,out,snapshot,profile=sys.argv[1:]
    out,snapshot=Path(out),Path(snapshot)
    env=dict(os.environ,WGBS_REAL_PROFILE=profile)
    if tool=='epykit':
        env['WGBS_EPYKIT_PROFILE']='paper_aligned' if profile=='paper_aligned' else 'baseline'
        command=[sys.executable,str(snapshot/'run_epykit_profile.py'),dataset,str(out)]
    else:command=['Rscript',str(snapshot/'run_r_profile.R'),tool,dataset,str(out)]
    subprocess.run(command,env=env,check=True)
    if tool=='methylKit' and profile=='paper_aligned':
        started=time.monotonic()
        count=clump_export(out/'dml.tsv',out/'dmr.tsv')
        (out/'dmr_candidates.tsv').write_bytes((out/'dmr.tsv').read_bytes())
        (out/'region_adapter.json').write_text(json.dumps(dict(calls=count,seconds=time.monotonic()-started,
            method='native methylKit CpG p-values; raw-p clumping',raw_p_threshold=1e-5,
            merge_bp=100,minlen_strict=50,minCG_strict=3,pct_sig=.5,region_fdr_control=False),indent=2)+'\n')

if __name__=='__main__':main()
