"""Bounded count blocks aligned to a predeclared coordinate universe."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

from ._beta_binomial import _counts

ALGORITHM_REVISION = "sorted-parquet-count-blocks-v1"
MAX_BLOCK_SIZE = 32768


def sample_parts(store, chrom, sample):
    return sorted((Path(store) / f"sample={sample}" / f"chrom={chrom}").glob("*.parquet"))


def canonical_sites(store, chrom, samples, *, intersect=False):
    """Coordinate-only union/intersection; counts and group effects are unread.

    Strand is descriptive, resolved to the first known value across samples.
    Counts are validated separately; duplicate input coordinates are fatal.
    """
    frames = []
    for sample in samples:
        parts = sample_parts(store, chrom, sample)
        if not parts:
            if intersect:
                return pl.DataFrame(schema={"pos": pl.Int32, "strand": pl.String})
            continue
        frames.append(pl.scan_parquet(parts).select("pos", "strand").unique(subset=["pos"]))
    if not frames:
        return pl.DataFrame(schema={"pos": pl.Int32, "strand": pl.String})
    combined = (
        pl.concat(frames)
        .group_by("pos")
        .agg(
            pl.len().alias("_support"),
            pl.col("strand").filter(pl.col("strand") != "*").first().fill_null("*").alias("strand"),
        )
    )
    if intersect:
        combined = combined.filter(pl.col("_support") == len(samples))
    return combined.select("pos", "strand").sort("pos").collect(engine="streaming")


def _sample_batches(store, chrom, sample, block_size):
    previous = None
    for path in sample_parts(store, chrom, sample):
        for batch in pq.ParquetFile(path).iter_batches(
            batch_size=block_size, columns=["pos", "N_meth", "coverage"], use_threads=False
        ):
            pos = batch.column(0).to_numpy(zero_copy_only=False)
            m, n = _counts(
                batch.column(1).to_numpy(zero_copy_only=False),
                batch.column(2).to_numpy(zero_copy_only=False),
                1,
            )
            if (
                not np.issubdtype(pos.dtype, np.integer)
                or np.any(np.diff(pos) <= 0)
                or (len(pos) and previous is not None and pos[0] <= previous)
            ):
                raise ValueError(
                    f"{sample}/{chrom} coordinates must be strictly increasing; "
                    "duplicate or unsorted records are not supported"
                )
            if len(pos):
                previous = pos[-1]
                yield pos, m, n


class _Cursor:
    def __init__(self, batches):
        self.batches = iter(batches)
        self.current = next(self.batches, None)
        self.offset = 0

    def fill(self, positions, m, n):
        end = positions[-1]
        while self.current is not None:
            pos, source_m, source_n = self.current
            stop = np.searchsorted(pos, end, side="right")
            part = pos[self.offset : stop]
            if len(part):
                index = np.searchsorted(positions, part)
                match = index < len(positions)
                match[match] &= positions[index[match]] == part[match]
                m[index[match]] = source_m[self.offset : stop][match]
                n[index[match]] = source_n[self.offset : stop][match]
            self.offset = stop
            if stop < len(pos):
                break
            self.current, self.offset = next(self.batches, None), 0

    def validate_remaining(self):
        for _ in self.batches:
            pass


def iter_count_blocks(store, chrom, samples, canonical_positions, *, block_size=32768):
    """Yield (coordinates, M, N) with <=32768 rows and missing counts at N=0.

    Inputs follow sorted single-coordinate CpG stores. Arrow input batches and
    output count matrices are bounded by block_size; no whole-chromosome
    sample matrix is built. Decoding also depends on source Parquet row groups.
    Consume the iterator fully to validate records beyond the declared family.
    """
    if not isinstance(block_size, (int, np.integer)) or not 1 <= block_size <= MAX_BLOCK_SIZE:
        raise ValueError("score block_size must be in [1,32768]")
    positions = np.asarray(canonical_positions)
    if (
        positions.ndim != 1
        or not np.issubdtype(positions.dtype, np.integer)
        or np.any(np.diff(positions) <= 0)
    ):
        raise ValueError("Canonical positions must be a strictly increasing integer vector")
    if len(set(samples)) != len(samples):
        raise ValueError("Samples must be unique")
    cursors = [_Cursor(_sample_batches(store, chrom, s, block_size)) for s in samples]
    for start in range(0, len(positions), block_size):
        block = positions[start : start + block_size]
        m, n = (
            np.zeros((len(block), len(samples)), dtype=np.int64),
            np.zeros((len(block), len(samples)), dtype=np.int64),
        )
        for j, cursor in enumerate(cursors):
            cursor.fill(block, m[:, j], n[:, j])
        yield block, m, n
    for cursor in cursors:
        cursor.validate_remaining()
