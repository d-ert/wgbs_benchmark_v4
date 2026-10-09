"""Brute-force coordinate-only family and processing-boundary oracles."""

import importlib
import importlib.util

import numpy as np
import pytest


def geometry():
    assert importlib.util.find_spec("epykit._region_geometry") is not None, (
        "Region geometry is missing"
    )
    return importlib.import_module("epykit._region_geometry")


def brute(positions, scales=(3, 5, 10, 20, 40), cap=40, gap=500, width=2000):
    parents, children = [], []
    left = 0
    while left < len(positions):
        right = left + 1
        while right < len(positions):
            if (
                right - left >= cap
                or positions[right] - positions[right - 1] > gap
                or positions[right] - positions[left] + 1 > width
            ):
                break
            right += 1
        parent = len(parents)
        parents.append((left, right, int(positions[left]), int(positions[right - 1]) + 1))
        for k in sorted(set(scales)):
            if right - left < k:
                continue
            anchors = list(range(left, right - k + 1, max(1, k // 2)))
            if anchors[-1] != right - k:
                anchors.append(right - k)
            for start in anchors:
                children.append(
                    (
                        start,
                        start + k,
                        int(positions[start]),
                        int(positions[start + k - 1]) + 1,
                        k,
                        parent,
                    )
                )
        left = right
    return parents, children


@pytest.mark.parametrize(
    "positions",
    [
        [],
        [10, 20],
        [10, 20, 30],
        list(range(1, 423, 7)),
        [100, 600, 1100],
        [100, 601, 1101],
        list(range(0, 5000, 103)),
    ],
)
def test_complete_family_matches_brute_force_for_small_and_edge_chromosomes(positions):
    g = geometry()
    p = np.array(positions, dtype=np.int64)
    expected_parents, expected_children = brute(p)
    family = g.parent_regions(p)
    actual_parents = [
        (int(r["left"]), int(r["right"]), int(r["start"]), int(r["end"])) for r in family.parents
    ]
    actual_children = [
        (
            int(r["left"]),
            int(r["right"]),
            int(r["start"]),
            int(r["end"]),
            int(r["scale"]),
            int(r["parent_index"]),
        )
        for r in family.intervals
    ]
    assert actual_parents == expected_parents
    assert actual_children == expected_children
    np.testing.assert_array_equal(g.multiscale_intervals(p), family.intervals)
    assert family.hypothesis_counts == {
        "intervals": len(expected_children),
        "parents": len(expected_parents),
    }


def test_exact_width_and_gap_limits_are_inclusive():
    g = geometry()
    a = g.parent_regions([0, 999, 1999], max_gap=2000, max_width=2000)
    b = g.parent_regions([0, 1000, 2000], max_gap=2000, max_width=2000)
    assert len(a.intervals) == 1 and a.intervals["end"][0] - a.intervals["start"][0] == 2000
    assert not len(b.intervals)
    assert len(g.multiscale_intervals([100, 600, 1100])) == 1
    assert not len(g.multiscale_intervals([100, 601, 1101]))


def test_half_step_tail_and_scale_ids_are_unique_and_order_invariant():
    g = geometry()
    p = np.arange(40, dtype=np.int64) * 10
    a = g.multiscale_intervals(p, scales=(5, 3, 10, 20, 40, 5))
    b = g.multiscale_intervals(p, scales=(3, 5, 10, 20, 40))
    np.testing.assert_array_equal(a, b)
    starts = a["left"][a["scale"] == 5].tolist()
    assert starts == list(range(0, 36, 2)) + [35]
    assert len(np.unique(a["test_id"])) == len(a)
    assert np.issubdtype(a["test_id"].dtype, np.integer)


def test_every_child_has_one_fixed_parent_and_forbidden_gaps_are_never_crossed():
    g = geometry()
    rng = np.random.default_rng(401)
    p = np.cumsum(rng.integers(1, 602, size=400))
    family = g.parent_regions(p)
    for child in family.intervals:
        parent = family.parents[child["parent_index"]]
        assert parent["left"] <= child["left"] < child["right"] <= parent["right"]
        assert parent["parent_id"] == child["parent_id"]
        assert np.diff(p[child["left"] : child["right"]]).max(initial=0) <= 500
        assert child["end"] - child["start"] <= 2000


def test_parent_aligned_processing_chunks_preserve_test_ids_and_site_indices():
    g = geometry()
    p = np.arange(153, dtype=np.int64) * 13
    full = g.parent_regions(p)
    parts = []
    for parent in full.parents:
        left, right = int(parent["left"]), int(parent["right"])
        intervals = g.multiscale_intervals(p[left:right]).copy()
        intervals["left"] += left
        intervals["right"] += left
        intervals["parent_index"] += int(np.searchsorted(full.parents["left"], left))
        parts.append(intervals)
    np.testing.assert_array_equal(np.concatenate(parts), full.intervals)


@pytest.mark.parametrize("positions", [[1, 1, 2], [2, 1], [-1, 0, 2], [1.0, 2.0, 3.0], [[1, 2]]])
def test_invalid_coordinates_are_refused(positions):
    with pytest.raises(ValueError):
        geometry().multiscale_intervals(positions)


@pytest.mark.parametrize(
    "kwargs",
    [{"scales": ()}, {"scales": (0, 3)}, {"scales": (2.5,)}, {"max_gap": -1}, {"max_width": 0}],
)
def test_invalid_geometry_parameters_are_refused(kwargs):
    with pytest.raises(ValueError):
        geometry().multiscale_intervals([1, 2, 3], **kwargs)
