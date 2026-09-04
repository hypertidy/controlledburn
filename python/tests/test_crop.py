# test_crop.py -- tests for BurnResult.crop(): the sparse sub-window
# crop that mirrors the R package's crop_burn(). Covers overlap,
# edge-aligned windows, non-overlap, run clipping/re-basing, the
# crop().materialize() tile chain, and argument validation.
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np
import pytest
import shapely

import controlledburn as cb

# extent is (xmin, xmax, ymin, ymax); shape is (nrow, ncol). Cell size 1.
G10 = dict(extent=(0, 10, 0, 10), shape=(10, 10))
BOX = shapely.box(2.5, 4.5, 6.5, 8.5)


def wkb(geom):
    return shapely.to_wkb(geom)


def dense(r, **kw):
    return r.materialize(fn="sum", edge_policy="fraction", background=0.0, **kw)


# ---- overlap: shape, extent, re-basing ------------------------------

def test_crop_overlap_shape_and_extent():
    r = cb.burn([wkb(BOX)], **G10)
    sub = r.crop((3, 7, 5, 9))
    assert sub.shape == (4, 4)
    assert sub.extent == (3.0, 7.0, 5.0, 9.0)


def test_crop_rebases_indices_to_zero():
    r = cb.burn([wkb(BOX)], **G10)
    sub = r.crop((3, 7, 5, 9))
    nrow, ncol = sub.shape
    if len(sub.runs):
        assert sub.runs["row"].min() >= 0 and sub.runs["row"].max() < nrow
        assert sub.runs["col_start"].min() >= 0
        assert sub.runs["col_end"].max() <= ncol
    if len(sub.edges):
        assert sub.edges["row"].min() >= 0 and sub.edges["row"].max() < nrow
        assert sub.edges["col"].min() >= 0 and sub.edges["col"].max() < ncol


# ---- crop().materialize() equals slicing the full dense array --------

def test_crop_matches_full_slice():
    r = cb.burn([wkb(BOX)], **G10)
    full = dense(r)
    # window x in [3, 7], y in [5, 9] -> rows 1..4, cols 3..6 (0-based)
    sub = dense(r.crop((3, 7, 5, 9)))
    np.testing.assert_allclose(sub, full[1:5, 3:7])


def test_crop_chain_reads_as_tile_workflow():
    r = cb.burn([wkb(BOX)], **G10)
    tile = r.crop((3, 7, 5, 9)).materialize(fn="sum", edge_policy="fraction",
                                            background=0.0)
    assert tile.shape == (4, 4)


# ---- edge-aligned windows (snap-to-boundary / eps behaviour) ---------

def test_crop_edge_aligned_window():
    # window edges fall exactly on cell boundaries (x=2,6 ; y=4,8)
    r = cb.burn([wkb(BOX)], **G10)
    sub = r.crop((2, 6, 4, 8))
    assert sub.shape == (4, 4)
    assert sub.extent == (2.0, 6.0, 4.0, 8.0)
    np.testing.assert_allclose(dense(sub), dense(r)[2:6, 2:6])


def test_crop_full_window_is_identity():
    r = cb.burn([wkb(BOX)], **G10)
    sub = r.crop((0, 10, 0, 10))
    assert sub.shape == (10, 10)
    assert sub.extent == (0.0, 10.0, 0.0, 10.0)
    np.testing.assert_array_equal(dense(sub), dense(r))


def test_crop_window_larger_than_grid_clamps():
    r = cb.burn([wkb(BOX)], **G10)
    sub = r.crop((-5, 15, -5, 15))
    assert sub.shape == (10, 10)
    assert sub.extent == (0.0, 10.0, 0.0, 10.0)
    np.testing.assert_array_equal(dense(sub), dense(r))


# ---- run clipping and re-basing on a full-grid polygon --------------

def test_crop_clips_and_rebases_runs():
    # a polygon covering the whole grid: every row is one run, col 0..9
    r = cb.burn([wkb(shapely.box(0, 0, 10, 10))], **G10)
    assert len(r.edges) == 0            # grid-aligned -> pure interior
    sub = r.crop((3, 7, 3, 7))
    assert sub.shape == (4, 4)
    assert sub.runs["col_start"].min() >= 0
    assert sub.runs["col_end"].max() <= 4
    assert sub.runs["row"].min() >= 0 and sub.runs["row"].max() < 4
    np.testing.assert_allclose(dense(sub), dense(r)[3:7, 3:7])


# ---- non-overlapping window: warn + empty tables --------------------

def test_crop_non_overlap_warns_and_returns_empty():
    r = cb.burn([wkb(BOX)], **G10)
    with pytest.warns(UserWarning, match="does not overlap"):
        sub = r.crop((100, 200, 100, 200))
    assert sub.shape == (0, 0)
    assert sub.extent == (100.0, 200.0, 100.0, 200.0)
    assert len(sub.runs) == 0 and len(sub.edges) == 0
    assert len(sub.lines) == 0 and len(sub.points) == 0


# ---- non-polygon tables are cropped and re-based too ----------------

def test_crop_line_table_cropped_and_rebased():
    line = shapely.LineString([(0.5, 0.5), (9.5, 9.5)])
    r = cb.burn([wkb(line)], **G10)
    assert len(r.lines) > 0
    sub = r.crop((3, 7, 3, 7))
    assert sub.shape == (4, 4)
    assert len(sub.lines) > 0
    assert len(sub.lines) <= len(r.lines)
    assert sub.lines["row"].min() >= 0 and sub.lines["row"].max() < 4
    assert sub.lines["col"].min() >= 0 and sub.lines["col"].max() < 4


def test_crop_point_table_cropped_and_rebased():
    pts = [shapely.Point(1.5, 1.5), shapely.Point(5.5, 5.5)]
    r = cb.burn([wkb(p) for p in pts], **G10)
    sub = r.crop((4, 8, 4, 8))
    # only the (5.5, 5.5) point falls inside the window
    assert len(sub.points) == 1
    assert sub.points["row"].min() >= 0 and sub.points["row"].max() < 4
    assert sub.points["col"].min() >= 0 and sub.points["col"].max() < 4


# ---- argument validation -------------------------------------------

def test_crop_requires_extent_and_shape():
    r = cb.burn([wkb(BOX)], **G10)
    bare = cb.BurnResult(r.runs, r.edges, r.lines, r.points)  # extent/shape None
    with pytest.raises(ValueError, match="extent and shape"):
        bare.crop((3, 7, 5, 9))


def test_crop_bad_target_length_raises():
    r = cb.burn([wkb(BOX)], **G10)
    with pytest.raises(ValueError, match="xmin, xmax, ymin, ymax"):
        r.crop((3, 7, 5))
