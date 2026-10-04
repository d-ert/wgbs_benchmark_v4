"""Check that every calibrated coordinate is a forward-strand GRCh38 CpG."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from wgbs_v3.provenance import digest


def fasta_index(path: Path) -> dict:
    index = {}
    with path.open("rb") as stream:
        chrom = None
        while True:
            offset = stream.tell()
            line = stream.readline()
            if not line:
                break
            if line.startswith(b">"):
                chrom = line[1:].split()[0].decode()
                index[chrom] = {"offset": None, "line_bases": None, "line_bytes": None,
                                "length": 0}
            elif chrom is not None:
                entry = index[chrom]
                if entry["offset"] is None:
                    entry["offset"] = offset
                    entry["line_bases"] = len(line.rstrip(b"\r\n"))
                    entry["line_bytes"] = len(line)
                entry["length"] += len(line.rstrip(b"\r\n"))
    return index


class ReferenceGenome:
    def __init__(self, path: Path):
        self.path = path
        self.index = fasta_index(path)
        self.bytes = np.memmap(path, mode="r", dtype=np.uint8)

    def is_cpg(self, chrom: str, pos: np.ndarray) -> np.ndarray:
        name = chrom.removeprefix("chr")
        entry = self.index.get(name)
        if entry is None:
            raise ValueError(f"Reference missing chromosome {chrom}")
        valid = (pos >= 1) & (pos < entry["length"])
        result = np.zeros(len(pos), dtype=bool)
        zero = pos[valid] - 1
        offset = entry["offset"] + (zero // entry["line_bases"]) * entry["line_bytes"] + zero % entry["line_bases"]
        next_pos = zero + 1
        offset_g = entry["offset"] + (next_pos // entry["line_bases"]) * entry["line_bytes"] + next_pos % entry["line_bases"]
        result[valid] = (self.bytes[offset] == ord("C")) & (self.bytes[offset_g] == ord("G"))
        return result


def validate(calibration: Path, fasta: Path) -> dict:
    manifest = json.loads((calibration / "manifest.json").read_text())
    reference = ReferenceGenome(fasta)
    bad = 0
    total = 0
    examples = []
    for filename in manifest["parquet_files"]:
        table = pq.read_table(filename, columns=["chrom", "pos"])
        if table.num_rows == 0:
            continue
        chrom = table["chrom"][0].as_py()
        pos = table["pos"].to_numpy().astype(np.int64)
        mismatch = ~reference.is_cpg(chrom, pos)
        bad += int(mismatch.sum())
        total += len(pos)
        for p in pos[np.flatnonzero(mismatch)[:max(0, 10 - len(examples))]]:
            examples.append(f"{chrom}:{p}")
    result = {"reference_sha256": digest(fasta), "calibration_sha256": digest(calibration / "manifest.json"),
              "checked": total, "non_cpg": bad, "examples": examples,
              "passed": total > 0 and bad == 0}
    (calibration / "reference_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["passed"]:
        raise ValueError(f"Reference check failed: {bad}/{total} non-CpGs; {examples}")
    return result
