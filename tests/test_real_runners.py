import ast
import subprocess
from pathlib import Path
from wgbs_real.runners import build_r_runner,build_epykit_runner

ROOT=Path(__file__).resolve().parents[1]

def test_native_sources_remain_untouched_and_adapters_parse(tmp_path):
    original=(ROOT/'run_r_tool.R').read_bytes()
    r=tmp_path/'r.R';p=tmp_path/'p.py'
    build_r_runner(ROOT/'run_r_tool.R',r)
    build_epykit_runner(ROOT/'run_epykit.py',p)
    assert (ROOT/'run_r_tool.R').read_bytes()==original
    ast.parse(p.read_text())
    result=subprocess.run(['Rscript','-e',f'parse(file="{r}");cat("Parsed\\n")'],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
