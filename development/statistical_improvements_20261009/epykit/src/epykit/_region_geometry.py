"""Coordinate-only nested multiscale families; no effect/seed selection.

Site indices and genomic bounds are half-open. Genomic coordinates retain the
input point origin (the package's CpG stores normally use 1-based positions).
Compact numeric IDs are stable genomic-bound pairs, namespaced by chromosome
and parent/window hypothesis type when exported. They are not merged calls.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

ALGORITHM_REVISION = "coordinate-parent-nested-halfstep-v1"
DEFAULT_SCALES = (3, 5, 10, 20, 40)
INTERVAL_DTYPE = np.dtype(
    [
        ("left", "i8"),
        ("right", "i8"),
        ("start", "i8"),
        ("end", "i8"),
        ("scale", "i4"),
        ("parent_index", "i8"),
        ("parent_id", "u8"),
        ("test_id", "u8"),
    ]
)
PARENT_DTYPE = np.dtype(
    [
        ("left", "i8"),
        ("right", "i8"),
        ("start", "i8"),
        ("end", "i8"),
        ("n_cpgs", "i4"),
        ("parent_id", "u8"),
    ]
)


@dataclass(frozen=True)
class RegionFamily:
    parents: np.ndarray
    intervals: np.ndarray

    @property
    def hypothesis_counts(self):
        return dict(intervals=len(self.intervals), parents=len(self.parents))


def _positions(positions):
    p = np.asarray(positions)
    if p.ndim != 1:
        raise ValueError("Positions must be a strictly increasing vector")
    if len(p) and (
        not np.issubdtype(p.dtype, np.integer)
        or np.any(p < 0)
        or np.any(p >= np.iinfo(np.int32).max)
        or np.any(np.diff(p) <= 0)
    ):
        raise ValueError(
            "Positions must be nonnegative strictly increasing integers "
            "within the canonical Int32 chromosome domain"
        )
    return np.ascontiguousarray(p, dtype=np.int64)


def _parameters(scales, max_gap, max_width):
    values = tuple(scales)
    if not values or any(not isinstance(k, (int, np.integer)) or k < 1 for k in values):
        raise ValueError("Scales must be positive integer CpG counts")
    if (
        not isinstance(max_gap, (int, np.integer))
        or max_gap < 0
        or not isinstance(max_width, (int, np.integer))
        or max_width < 1
    ):
        raise ValueError("Require nonnegative max_gap and positive max_width")
    return np.array(sorted(set(values)), dtype=np.int64)


def _pair_ids(start, end):
    # Cantor pairing is injective. The Int32 coordinate domain guarantees the
    # product fits uint64 even before division; no probabilistic hash collision.
    a, b = np.asarray(start, dtype=np.uint64), np.asarray(end, dtype=np.uint64)
    total = a + b
    return total * (total + np.uint64(1)) // np.uint64(2) + b


@njit(cache=True)
def _parent_bounds(p, cap, gap, width):
    bounds = np.empty((len(p), 2), dtype=np.int64)
    left, count = 0, 0
    while left < len(p):
        right = left + 1
        while right < len(p):
            if (
                right - left >= cap
                or p[right] - p[right - 1] > gap
                or p[right] - p[left] + 1 > width
            ):
                break
            right += 1
        bounds[count, 0], bounds[count, 1] = left, right
        count += 1
        left = right
    return bounds[:count]


@njit(cache=True)
def _child_bounds(bounds, scales):
    count = 0
    for j in range(len(bounds)):
        left, right = bounds[j, 0], bounds[j, 1]
        for k in scales:
            if right - left < k:
                continue
            step = max(1, k // 2)
            count += (right - k - left) // step + 1
            if (right - k - left) % step:
                count += 1
    result = np.empty((count, 4), dtype=np.int64)
    offset = 0
    for j in range(len(bounds)):
        left, right = bounds[j, 0], bounds[j, 1]
        for k in scales:
            if right - left < k:
                continue
            step = max(1, k // 2)
            for start in range(left, right - k + 1, step):
                result[offset, 0], result[offset, 1] = start, start + k
                result[offset, 2], result[offset, 3] = k, j
                offset += 1
            if (right - k - left) % step:
                result[offset, 0], result[offset, 1] = right - k, right
                result[offset, 2], result[offset, 3] = k, j
                offset += 1
    return result


def parent_regions(positions, *, max_cpgs=40, max_gap=500, max_width=2000, scales=DEFAULT_SCALES):
    """Disjoint greedy coordinate parents and their complete child family.

    A parent ends at a size, gap or width limit; small tails remain parents
    even when no child fits. Children never cross parent cuts. Each scale
    anchors every max(1,k//2) sites, with a final right-aligned tail anchor.
    Counts report unavailable parents too; inference treats no-child parents
    as p=1. Scales larger than the parent cap produce no children.
    """
    p = _positions(positions)
    ks = _parameters(scales, max_gap, max_width)
    if not isinstance(max_cpgs, (int, np.integer)) or max_cpgs < 1:
        raise ValueError("max_cpgs must be a positive integer")
    bounds = _parent_bounds(p, max_cpgs, max_gap, max_width)
    parents = np.empty(len(bounds), dtype=PARENT_DTYPE)
    if len(bounds):
        parents["left"], parents["right"] = bounds[:, 0], bounds[:, 1]
        parents["start"], parents["end"] = p[bounds[:, 0]], p[bounds[:, 1] - 1] + 1
        parents["n_cpgs"] = bounds[:, 1] - bounds[:, 0]
        parents["parent_id"] = _pair_ids(parents["start"], parents["end"])
    child = _child_bounds(bounds, ks)
    intervals = np.empty(len(child), dtype=INTERVAL_DTYPE)
    if len(child):
        intervals["left"], intervals["right"] = child[:, 0], child[:, 1]
        intervals["start"], intervals["end"] = p[child[:, 0]], p[child[:, 1] - 1] + 1
        intervals["scale"], intervals["parent_index"] = child[:, 2], child[:, 3]
        intervals["parent_id"] = parents["parent_id"][child[:, 3]]
        intervals["test_id"] = _pair_ids(intervals["start"], intervals["end"])
    return RegionFamily(parents=parents, intervals=intervals)


def multiscale_intervals(positions, *, scales=DEFAULT_SCALES, max_gap=500, max_width=2000):
    """Return the fixed nested interval family with stable IDs.

    The default parent cap is 40 CpGs, expanding if a larger scale is explicitly
    requested. Use parent_regions to declare a different cap and retain both
    parent and interval hypotheses. No significant CpG seed is inspected.
    """
    ks = _parameters(scales, max_gap, max_width)
    return parent_regions(
        positions, max_cpgs=max(40, int(ks.max())), scales=ks, max_gap=max_gap, max_width=max_width
    ).intervals
