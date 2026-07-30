# controlledburn -- scanline polygon/line/point rasterization with exact
# coverage fractions, O(perimeter) sparse output.
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from . import _core

__all__ = ["burn", "materialize", "BurnResult"]
__version__ = "0.1.0"

_RUNS_DTYPE = np.dtype(
    [("row", "i4"), ("col_start", "i4"), ("col_end", "i4"), ("id", "i4")]
)
_EDGES_DTYPE = np.dtype(
    [("row", "i4"), ("col", "i4"), ("fraction", "f4"), ("id", "i4")]
)
_LINES_DTYPE = np.dtype(
    [("row", "i4"), ("col", "i4"), ("length", "f4"), ("id", "i4")]
)
_POINTS_DTYPE = np.dtype([("row", "i4"), ("col", "i4"), ("id", "i4")])


def _structured(cols: dict, dtype: np.dtype) -> np.ndarray:
    n = len(cols[dtype.names[0]])
    out = np.empty(n, dtype=dtype)
    for name in dtype.names:
        out[name] = cols[name]
    return out


@dataclass
class BurnResult:
    """Sparse rasterization output: the four-table contract.

    Tables are numpy structured arrays. Row/col indices are 1-BASED with
    row 1 at the TOP of the grid -- identical values to the R package,
    so cross-language fixtures compare directly. ``pandas.DataFrame(r.runs)``
    works as-is.

    runs   : interior RLE          (row, col_start, col_end, id)
    edges  : polygon boundary cells (row, col, fraction, id), fraction in [0, 1]
    lines  : line cells            (row, col, length, id), length in CRS units
    points : point cells           (row, col, id)
    """

    runs: np.ndarray
    edges: np.ndarray
    lines: np.ndarray
    points: np.ndarray
    bounds: tuple = field(default=None)
    shape: tuple = field(default=None)


def burn(
    wkb: Sequence,
    bounds: Sequence[float],
    shape: Sequence[int],
    mode: str = "coverage",
) -> BurnResult:
    """Burn WKB geometries onto a regular grid, returning sparse tables.

    Parameters
    ----------
    wkb : sequence of bytes-like
        One WKB blob per geometry (ISO WKB or EWKB; both byte orders;
        Z/M ordinates are skipped). ``shapely.to_wkb(geom)`` output works
        directly. Geometry k (0-based position) is assigned id k + 1 in
        the output tables. ``None`` entries are skipped.
    bounds : (xmin, ymin, xmax, ymax)
        Grid extent, rasterio-style bounds ordering. (Note this differs
        from the R package's ``extent = c(xmin, xmax, ymin, ymax)``.)
    shape : (nrow, ncol)
        Grid dimensions, numpy-style ordering. (The R package's
        ``dimension`` is ``c(ncol, nrow)``.)
    mode : str
        ``"coverage"`` (default) computes exact coverage fractions for
        polygon boundary cells. ``"approx"`` uses the cell-center rule
        (fasterize semantics): boundary cells are included as full run
        cells iff the cell center is inside the polygon. No edges are
        produced for polygons in approx mode. Lines and points are
        unaffected by mode.

    Non-fatal problems (unparseable WKB, GeometryCollection input) are
    raised as warnings and the offending geometry is skipped;
    GeometryCollections must be split into homogeneous groups upstream.
    """
    if mode not in ("coverage", "approx"):
        raise ValueError("mode must be 'coverage' or 'approx'")
    xmin, ymin, xmax, ymax = (float(v) for v in bounds)
    nrow, ncol = (int(v) for v in shape)

    raw = _core.burn_wkb(list(wkb), xmin, ymin, xmax, ymax, ncol, nrow, mode)

    for geom_index, message in raw["notes"]:
        warnings.warn(f"geometry {geom_index}: {message}", stacklevel=2)

    return BurnResult(
        runs=_structured(raw["runs"], _RUNS_DTYPE),
        edges=_structured(raw["edges"], _EDGES_DTYPE),
        lines=_structured(raw["lines"], _LINES_DTYPE),
        points=_structured(raw["points"], _POINTS_DTYPE),
        bounds=(xmin, ymin, xmax, ymax),
        shape=(nrow, ncol),
    )


def materialize(
    result: BurnResult,
    shape: Optional[Sequence[int]] = None,
    values: Optional[Sequence[float]] = None,
    fn: str = "last",
    edge_policy: str = "threshold",
    threshold: float = 0.5,
    background: float = np.nan,
) -> np.ndarray:
    """Materialize polygon output (runs + edges) into a dense 2D array.

    This is the fasterize-style consumer of the sparse tables. Row 1 of
    the tables maps to array row 0 (top of grid).

    Parameters
    ----------
    result : BurnResult
    shape : (nrow, ncol), optional
        Defaults to ``result.shape``.
    values : sequence of float, optional
        Burn value per geometry id (value for id k is ``values[k - 1]``).
        Defaults to burning the id itself.
    fn : str
        Per-pixel reduction when several geometries touch the same cell:
        one of "first", "last", "sum", "min", "max", "count", "any".
    edge_policy : str
        Boundary-cell handling. "threshold": include the cell iff its
        coverage fraction >= ``threshold`` (approximates -- but does not
        equal -- the fasterize cell-center rule). "fraction": write
        value * fraction (exact, area-conserving).
    threshold : float
    background : float
        Fill value for untouched cells (default NaN).
    """
    if shape is None:
        shape = result.shape
    if shape is None:
        raise ValueError("shape is required when result.shape is not set")
    nrow, ncol = (int(v) for v in shape)

    return _core.materialize(
        np.ascontiguousarray(result.runs["row"]),
        np.ascontiguousarray(result.runs["col_start"]),
        np.ascontiguousarray(result.runs["col_end"]),
        np.ascontiguousarray(result.runs["id"]),
        np.ascontiguousarray(result.edges["row"]),
        np.ascontiguousarray(result.edges["col"]),
        np.ascontiguousarray(result.edges["fraction"]),
        np.ascontiguousarray(result.edges["id"]),
        ncol,
        nrow,
        [] if values is None else [float(v) for v in values],
        fn,
        edge_policy,
        float(threshold),
        float(background),
    )
