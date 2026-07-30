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

__all__ = ["burn", "materialize", "as_wkb_list", "resolve_grid", "BurnResult"]
__version__ = "0.2.0"

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


# --------------------------------------------------------------------
# input coercion -- the Python analogue of the R package's
# as_wkb_list() S3 generic
# --------------------------------------------------------------------

def as_wkb_list(x) -> list:
    """Coerce assorted geometry containers to a list of WKB blobs.

    Accepts: raw WKB bytes, a single shapely geometry, a shapely
    geometry array, a geopandas GeoSeries or GeoDataFrame, or any
    sequence mixing WKB bytes, shapely geometries and ``None``.
    ``None`` entries are preserved (burn() skips them but they still
    consume an id).

    shapely is imported only when shapely geometries are actually
    present in the input; WKB-only paths have no shapely dependency.
    """
    if isinstance(x, (bytes, bytearray, memoryview)):
        return [x]

    if hasattr(x, "geom_type"):                            # one shapely geom
        import shapely
        return [shapely.to_wkb(x)]

    if hasattr(x, "geometry") and hasattr(x, "columns"):   # GeoDataFrame
        x = x.geometry

    if hasattr(x, "to_wkb"):                               # GeoSeries etc.
        return list(x.to_wkb())

    out = []
    for g in x:
        if g is None or isinstance(g, (bytes, bytearray, memoryview)):
            out.append(g)
        elif hasattr(g, "geom_type"):
            import shapely
            out.append(shapely.to_wkb(g))
        else:
            raise TypeError(
                f"cannot interpret {type(g).__name__} as geometry; "
                "pass WKB bytes, a shapely geometry, or None"
            )
    return out


# --------------------------------------------------------------------
# grid parameter resolution -- the Python analogue of the R package's
# .resolve_grid_params()
# --------------------------------------------------------------------

def resolve_grid(
    geoms,
    bounds: Optional[Sequence[float]] = None,
    shape: Optional[Sequence[int]] = None,
    resolution=None,
    size: int = 256,
):
    """Fill in ``bounds`` and ``shape`` the way the R package does.

    ``bounds`` defaults to the bounding box of ``geoms`` (requires
    shapely). ``shape`` defaults to a grid of at most ``size`` cells
    along the longer axis, preserving aspect ratio. ``resolution``
    (scalar or ``(dx, dy)``) computes ``shape`` from the extent and is
    mutually exclusive with ``shape``.

    Returns ``(bounds, shape)`` with bounds as
    ``(xmin, ymin, xmax, ymax)`` and shape as ``(nrow, ncol)``.
    """
    if bounds is None:
        try:
            import shapely
        except ImportError:
            raise ImportError(
                "deriving bounds from geometry requires shapely; "
                "pass bounds=(xmin, ymin, xmax, ymax) explicitly"
            ) from None
        wkb = [w for w in as_wkb_list(geoms) if w is not None]
        if not wkb:
            raise ValueError("cannot derive bounds from empty geometry input")
        b = np.atleast_2d(shapely.bounds(shapely.from_wkb(wkb)))
        bounds = (float(np.nanmin(b[:, 0])), float(np.nanmin(b[:, 1])),
                  float(np.nanmax(b[:, 2])), float(np.nanmax(b[:, 3])))

    xmin, ymin, xmax, ymax = (float(v) for v in bounds)
    if not (xmax > xmin and ymax > ymin):
        raise ValueError(
            "invalid extent: xmax must be > xmin, ymax must be > ymin"
        )

    if shape is not None and resolution is not None:
        raise ValueError("specify 'shape' or 'resolution', not both")

    if resolution is not None:
        res = np.broadcast_to(np.asarray(resolution, dtype=float), (2,))
        if not np.all(res > 0):
            raise ValueError("resolution must be positive")
        shape = (int(np.ceil((ymax - ymin) / res[1])),
                 int(np.ceil((xmax - xmin) / res[0])))

    if shape is None:
        w, h = xmax - xmin, ymax - ymin
        m = max(w, h)
        shape = (int(np.ceil(size * h / m)), int(np.ceil(size * w / m)))

    nrow, ncol = int(shape[0]), int(shape[1])
    if nrow <= 0 or ncol <= 0:
        raise ValueError("ncol and nrow must be positive")

    return (xmin, ymin, xmax, ymax), (nrow, ncol)


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

    def __repr__(self) -> str:
        # Summary, not a table dump: results can hold millions of rows.
        # Mirrors the R package's print.controlledburn.
        if self.shape is not None:
            nrow, ncol = self.shape
            grid = f"{ncol} x {nrow} grid"
            total = float(nrow) * float(ncol)
        else:
            grid = "grid unknown"
            total = None
        interior = int(
            (self.runs["col_end"] - self.runs["col_start"] + 1).sum()
        ) if len(self.runs) else 0
        touched = interior + len(self.edges) + len(self.lines) + len(self.points)
        ids = np.unique(np.concatenate([
            self.runs["id"], self.edges["id"],
            self.lines["id"], self.points["id"],
        ]))
        n = len(ids)
        out = [f"<BurnResult> {grid}, {n} "
               f"{'geometry' if n == 1 else 'geometries'}",
               f"  runs:   {len(self.runs)} ({interior} interior cells)",
               f"  edges:  {len(self.edges)} polygon boundary cells"]
        if len(self.lines):
            out.append(f"  lines:  {len(self.lines)} cells")
        if len(self.points):
            out.append(f"  points: {len(self.points)} cells")
        if total:
            out.append(f"  sparsity: {100.0 * (1.0 - touched / total):.1f}% empty")
        return "\n".join(out)


def burn(
    geoms,
    bounds: Optional[Sequence[float]] = None,
    shape: Optional[Sequence[int]] = None,
    resolution=None,
    mode: str = "coverage",
) -> BurnResult:
    """Burn geometries onto a regular grid, returning sparse tables.

    Parameters
    ----------
    geoms : geometry input
        WKB bytes (a single blob or a sequence), shapely geometries (a
        single geometry or a sequence), a geopandas GeoSeries or
        GeoDataFrame, or a sequence mixing WKB, shapely and ``None``.
        See :func:`as_wkb_list`. WKB may be ISO WKB or EWKB, either
        byte order; Z/M ordinates are skipped. Geometry k (0-based
        position) is assigned id k + 1 in the output tables. ``None``
        entries are skipped.
    bounds : (xmin, ymin, xmax, ymax), optional
        Grid extent, rasterio-style bounds ordering. Defaults to the
        bounding box of the input geometry (requires shapely). (Note
        this differs from the R package's
        ``extent = c(xmin, xmax, ymin, ymax)``.)
    shape : (nrow, ncol), optional
        Grid dimensions, numpy-style ordering. Defaults to a grid of at
        most 256 cells along the longer axis of the extent, preserving
        aspect ratio -- the same rule as the R package. Mutually
        exclusive with ``resolution``. (The R package's ``dimension``
        is ``c(ncol, nrow)``.)
    resolution : float or (dx, dy), optional
        Cell size; ``shape`` is computed as ``ceil(extent / resolution)``.
        Mutually exclusive with ``shape``.
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

    wkb = as_wkb_list(geoms)
    bounds, shape = resolve_grid(wkb, bounds, shape, resolution)
    xmin, ymin, xmax, ymax = bounds
    nrow, ncol = shape

    raw = _core.burn_wkb(list(wkb), xmin, ymin, xmax, ymax, ncol, nrow, mode)

    for geom_index, message in raw["notes"]:
        warnings.warn(f"geometry {geom_index}: {message}", stacklevel=2)

    return BurnResult(
        runs=_structured(raw["runs"], _RUNS_DTYPE),
        edges=_structured(raw["edges"], _EDGES_DTYPE),
        lines=_structured(raw["lines"], _LINES_DTYPE),
        points=_structured(raw["points"], _POINTS_DTYPE),
        bounds=bounds,
        shape=shape,
    )


# --------------------------------------------------------------------
# materialization
# --------------------------------------------------------------------

_FNS = ("first", "last", "sum", "min", "max", "count", "any")


def materialize(
    result: BurnResult,
    shape: Optional[Sequence[int]] = None,
    values: Optional[Sequence[float]] = None,
    fn: str = "last",
    edge_policy: str = "threshold",
    threshold: float = 0.5,
    background: float = np.nan,
) -> np.ndarray:
    """Materialize a BurnResult into a dense 2D array.

    Row 1 of the tables maps to array row 0 (top of grid).

    This is the fasterize-style POLYGON consumer of the sparse tables:
    ``values[k - 1]`` is burned for geometry id k (default: the id
    itself), combined per-pixel by ``fn``, with ``edge_policy``
    controlling boundary cells.

    Line and point input raise NotImplementedError rather than
    returning an empty array. What a dense line or point raster should
    mean (burn the length? a value? which reduction?) is an aggregation
    question this package deliberately does not answer yet; the sparse
    ``lines`` and ``points`` tables are the product, and dense
    consumption is a few lines of numpy over them. See the tracking
    issue for the design discussion. Mixed-kind results raise for the
    same reason, compounded: the measures have different units
    (dimensionless coverage, length in CRS units, counts).

    Parameters
    ----------
    result : BurnResult
    shape : (nrow, ncol), optional
        Defaults to ``result.shape``.
    values : sequence of float, optional
        Per-geometry burn value; entry for id k is ``values[k - 1]``.
        Defaults to burning the id itself.
    fn : str
        Per-pixel reduction when several records touch the same cell:
        one of "first", "last", "sum", "min", "max", "count", "any".
    edge_policy : str
        Polygon boundary-cell handling. "threshold": include the cell
        iff its coverage fraction >= ``threshold`` (approximates -- but
        does not equal -- the fasterize cell-center rule; for true
        parity use ``mode="approx"`` at burn time). "fraction": write
        value * fraction (exact, area-conserving).
    threshold : float
    background : float
        Fill value for untouched cells (default NaN).
    """
    if fn not in _FNS:
        raise ValueError(f"fn must be one of {'/'.join(_FNS)}")

    kinds = []
    if len(result.runs) or len(result.edges):
        kinds.append("polygon")
    if len(result.lines):
        kinds.append("line")
    if len(result.points):
        kinds.append("point")
    if len(kinds) > 1:
        raise ValueError(
            f"mixed geometry kinds present ({', '.join(kinds)}); the "
            "measures have different units (coverage fraction, length "
            "in CRS units, count) and materialize() handles polygon "
            "output only. Burn each kind separately, or consume the "
            "sparse tables directly."
        )

    if shape is None:
        shape = result.shape
    if shape is None:
        raise ValueError("shape is required when result.shape is not set")
    nrow, ncol = (int(v) for v in shape)

    if kinds and kinds[0] != "polygon":
        raise NotImplementedError(
            f"materialize() handles polygon output (runs + edges) only; "
            f"this result contains {kinds[0]} records. Dense {kinds[0]} "
            "rasterization semantics are deliberately unresolved -- "
            "consume the sparse tables directly (numpy indexing on "
            f"result.{kinds[0]}s) or see the tracking issue."
        )

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
