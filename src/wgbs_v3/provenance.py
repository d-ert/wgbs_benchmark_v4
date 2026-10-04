"""Content hashes and immutable source snapshots for reproducible experiments."""
from __future__ import annotations

import importlib.metadata
import json
import platform
import shutil
from pathlib import Path

import hashlib


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def source_inventory(root: Path) -> dict:
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix in {".py", ".R", ".toml"} and "__pycache__" not in p.parts}


def dataset_inventory(dataset: Path) -> dict:
    # Hash content, not just the manifest that names it.
    paths = [dataset / "manifest.json", dataset / "samples.tsv", dataset / "dmr_truth.tsv"]
    for directory in ("truth_parts", "eligibility_parts", "counts"):
        paths.extend(sorted((dataset / directory).glob("*")))
    return {str(p.relative_to(dataset)): digest(p) for p in paths if p.is_file()}


def snapshot(root: Path, ep_source: Path, dataset: Path, out: Path) -> dict:
    target = out / "source_snapshot"
    target.mkdir(exist_ok=True)
    files = list((root / "src/wgbs_v3").glob("*.py")) + [root / "run_epykit.py", root / "run_r_tool.R", root / "benchmark.py"] + list((root / "analysis").glob("*.py")) + list((root / "analysis").glob("*.R"))
    for p in files:
        relative = p.relative_to(root)
        destination = target / "benchmark" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, destination)
    shutil.copytree(ep_source / "epykit", target / "epykit", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    versions = {}
    for name in ("numpy", "scipy", "pandas", "polars", "pyarrow", "statsmodels", "scikit-learn"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "unavailable"
    result = {"python": platform.python_version(), "platform": platform.platform(),
              "versions": versions, "epykit_source": str(ep_source.resolve()),
              "source_hashes": source_inventory(target), "dataset_hashes": dataset_inventory(dataset)}
    write_json(out / "provenance.json", result)
    return result
