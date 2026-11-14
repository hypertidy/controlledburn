# test_options.py -- burn() options from the Rust crate (0.3): merge_parts,
# coalesce_runs and id_offset, with the crate's defaults.
#
# Copyright (c) 2026 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np
import pytest
import shapely

import controlledburn as cb

G3 = dict(extent=(0, 3, 0, 3), shape=(3, 3))
HALVES = shapely.MultiPolygon([shapely.box(0, 0, 1.5, 3), shapely.box(1.5, 0, 3, 3)])


def _same(a, b):
    for k in ("runs", "edges", "lines", "points"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k


def test_merge_parts_makes_two_halves_a_full_cell():
    r = cb.burn(HALVES, **G3)
    assert len(r.edges) == 0
    assert len(r.runs) == 3  # one coalesced run per row
    assert (r.runs["col_end"] - r.runs["col_start"]).sum() == 9
    old = cb.burn(HALVES, **G3, merge_parts=False, coalesce_runs=False)
    assert len(old.edges) == 6  # 0.5 from each part in the middle column
    dense = dict(values=[1.0], fn="sum", background=0.0, edge_policy="fraction")
    full = cb.materialize(r, **dense)
    np.testing.assert_array_equal(full, np.ones((3, 3)))
    np.testing.assert_array_equal(cb.materialize(old, **dense), full)


def test_coalesce_runs_off_gives_one_cell_runs():
    sq = shapely.box(1, 0.5, 9, 9.5)  # sides on column walls
    g = dict(extent=(0, 10, 0, 10), shape=(10, 10))
    a = cb.burn(sq, **g)
    b = cb.burn(sq, **g, coalesce_runs=False)
    assert len(b.runs) > len(a.runs)
    assert (a.runs["col_end"] - a.runs["col_start"]).sum() == \
        (b.runs["col_end"] - b.runs["col_start"]).sum()


def test_id_offset_concatenates():
    boxes = [shapely.box(i, 0.2, i + 0.8, 2.7) for i in range(8)]
    g = dict(extent=(0, 8, 0, 3), shape=(6, 16))
    whole = cb.burn(boxes, **g)
    a = cb.burn(boxes[:5], **g)
    b = cb.burn(boxes[5:], **g, id_offset=5)
    for k in ("runs", "edges"):
        assert np.array_equal(np.concatenate([getattr(a, k), getattr(b, k)]),
                              getattr(whole, k))
    with pytest.raises(ValueError):
        cb.burn(boxes, **g, id_offset=-1)


def test_arrow_path_takes_the_options():
    pa = pytest.importorskip("pyarrow")
    pytest.importorskip("nanoarrow")
    ws = [shapely.to_wkb(HALVES), shapely.to_wkb(shapely.box(0.2, 0.2, 2.2, 2.2))]
    arr = pa.chunked_array([pa.array(ws[:1], type=pa.binary()),
                            pa.array(ws[1:], type=pa.binary())])
    for opts in (dict(), dict(merge_parts=False, coalesce_runs=False, id_offset=4)):
        _same(cb.burn(ws, **G3, **opts), cb.burn(arr, **G3, **opts))
