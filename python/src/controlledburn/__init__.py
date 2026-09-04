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
__version__ = "0.3.0"

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

    if hasattr(x, "geometry") and hasattr(x, "columns"):   # GeoDataFrame
        x = x.geometry

    if hasattr(x, "to_wkb") and hasattr(x, "__len__"):     # GeoSeries etc.
        return list(x.to_wkb())

    if hasattr(x, "geom_type"):                            # one shapely geom
        import shapely
        return [shapely.to_wkb(x)]

    if _is_arrow(x):                                       # Arrow C data interface
        return _arrow_to_wkb_list(x)

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
    extent: Optional[Sequence[float]] = None,
    shape: Optional[Sequence[int]] = None,
    resolution=None,
    size: int = 256,
):
    """Fill in ``extent`` and ``shape`` the way the R package does.

    ``extent`` defaults to the bounding box of ``geoms`` (requires
    shapely). ``shape`` defaults to a grid of at most ``size`` cells
    along the longer axis, preserving aspect ratio. ``resolution``
    (scalar or ``(dx, dy)``) computes ``shape`` from the extent and is
    mutually exclusive with ``shape``.

    Returns ``(extent, shape)`` with extent as
    ``(xmin, xmax, ymin, ymax)`` (matching the R package's extent
    ordering) and shape as ``(nrow, ncol)``.
    """
    if extent is None:
        try:
            import shapely
        except ImportError:
            raise ImportError(
                "deriving extent from geometry requires shapely; "
                "pass extent=(xmin, xmax, ymin, ymax) explicitly"
            ) from None
        wkb = [w for w in as_wkb_list(geoms) if w is not None]
        if not wkb:
            raise ValueError("cannot derive extent from empty geometry input")
        b = np.atleast_2d(shapely.bounds(shapely.from_wkb(wkb)))
        extent = (float(np.nanmin(b[:, 0])), float(np.nanmax(b[:, 2])),
                  float(np.nanmin(b[:, 1])), float(np.nanmax(b[:, 3])))

    xmin, xmax, ymin, ymax = (float(v) for v in extent)
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

    return (xmin, xmax, ymin, ymax), (nrow, ncol)


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
    extent: tuple = field(default=None)
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

    def crop(self, target: Sequence[float]) -> "BurnResult":
        """Crop to a target window, returning a new BurnResult.

        The sparse analogue of a raster crop and the Python counterpart
        of the R package's ``crop_burn()``: filter and clip the four
        tables (runs, edges, lines, points) to a sub-window, re-basing
        row/col indices to 1 and snapping the window outward to whole
        cells. No dense allocation -- just structured-array filtering.

        Chains with :meth:`materialize` to cut tiles from a single burn::

            tile = r.crop((150, 250, 150, 250)).materialize(fn="sum",
                                                             edge_policy="fraction")

        Parameters
        ----------
        target : (xmin, xmax, ymin, ymax)
            Target window in CRS units, same ordering as ``extent`` and
            as the R package's ``crop_burn()`` target. The window is
            snapped outward to cell boundaries.

        Returns
        -------
        BurnResult
            A result covering only the target window, with row/col
            indices re-based to 1 and ``extent``/``shape`` updated to
            the snapped window. A window that does not overlap the grid
            warns and returns empty tables with ``shape == (0, 0)``.
        """
        if self.extent is None or self.shape is None:
            raise ValueError("crop requires extent and shape to be set")
        if len(target) != 4:
            raise ValueError("target must be (xmin, xmax, ymin, ymax)")

        xmin, xmax, ymin, ymax = (float(v) for v in self.extent)
        nrow, ncol = int(self.shape[0]), int(self.shape[1])
        txmin, txmax, tymin, tymax = (float(v) for v in target)

        dx = (xmax - xmin) / ncol
        dy = (ymax - ymin) / nrow

        # 1-based inclusive row/col limits, snapped outward. The eps
        # nudge avoids floor/ceil flips when the target aligns exactly
        # with a cell boundary.
        eps = 1e-8
        col_lo = max(1, int(np.floor((txmin - xmin) / dx + eps)) + 1)
        col_hi = min(ncol, int(np.ceil((txmax - xmin) / dx - eps)))
        row_hi = min(nrow, int(np.ceil((ymax - tymin) / dy - eps)))
        row_lo = max(1, int(np.floor((ymax - tymax) / dy + eps)) + 1)

        if col_lo > col_hi or row_lo > row_hi:
            warnings.warn("target extent does not overlap the grid",
                          stacklevel=2)
            return BurnResult(
                runs=self.runs[:0].copy(),
                edges=self.edges[:0].copy(),
                lines=self.lines[:0].copy(),
                points=self.points[:0].copy(),
                extent=(txmin, txmax, tymin, tymax),
                shape=(0, 0),
            )

        new_extent = (
            xmin + (col_lo - 1) * dx,
            xmin + col_hi * dx,
            ymax - row_hi * dy,
            ymax - (row_lo - 1) * dy,
        )
        new_shape = (row_hi - row_lo + 1, col_hi - col_lo + 1)

        def crop_single(a):
            if len(a) == 0:
                return a.copy()
            keep = ((a["row"] >= row_lo) & (a["row"] <= row_hi) &
                    (a["col"] >= col_lo) & (a["col"] <= col_hi))
            out = a[keep].copy()
            out["row"] -= row_lo - 1
            out["col"] -= col_lo - 1
            return out

        def crop_runs(a):
            if len(a) == 0:
                return a.copy()
            keep = ((a["row"] >= row_lo) & (a["row"] <= row_hi) &
                    (a["col_end"] >= col_lo) & (a["col_start"] <= col_hi))
            out = a[keep].copy()
            out["row"] -= row_lo - 1
            out["col_start"] = np.maximum(out["col_start"], col_lo) - (col_lo - 1)
            out["col_end"] = np.minimum(out["col_end"], col_hi) - (col_lo - 1)
            return out

        return BurnResult(
            runs=crop_runs(self.runs),
            edges=crop_single(self.edges),
            lines=crop_single(self.lines),
            points=crop_single(self.points),
            extent=new_extent,
            shape=new_shape,
        )

    def materialize(
        self,
        shape: Optional[Sequence[int]] = None,
        values: Optional[Sequence[float]] = None,
        fn: str = "last",
        edge_policy: str = "threshold",
        threshold: float = 0.5,
        background: float = np.nan,
    ) -> np.ndarray:
        """Dense materialization of this result.

        Method form of the module-level :func:`materialize`, so the
        crop/materialize tile workflow reads as a chain::

            r.crop((150, 250, 150, 250)).materialize(fn="sum")

        See :func:`materialize` for the full parameter documentation.
        """
        return materialize(
            self, shape=shape, values=values, fn=fn,
            edge_policy=edge_policy, threshold=threshold,
            background=background,
        )


def burn(
    geoms,
    extent: Optional[Sequence[float]] = None,
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
    extent : (xmin, xmax, ymin, ymax), optional
        Grid extent, matching the R package's
        ``extent = c(xmin, xmax, ymin, ymax)`` ordering. Defaults to the
        bounding box of the input geometry (requires shapely).
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

    if _is_arrow(geoms) and not _is_shapely_like(geoms):
        return _burn_arrow(geoms, extent, shape, resolution, mode)

    wkb = as_wkb_list(geoms)
    extent, shape = resolve_grid(wkb, extent, shape, resolution)
    xmin, xmax, ymin, ymax = extent
    nrow, ncol = shape

    raw = _core.burn_wkb(list(wkb), xmin, ymin, xmax, ymax, ncol, nrow, mode)

    for geom_index, message in raw["notes"]:
        warnings.warn(f"geometry {geom_index}: {message}", stacklevel=2)

    return BurnResult(
        runs=_structured(raw["runs"], _RUNS_DTYPE),
        edges=_structured(raw["edges"], _EDGES_DTYPE),
        lines=_structured(raw["lines"], _LINES_DTYPE),
        points=_structured(raw["points"], _POINTS_DTYPE),
        extent=extent,
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
            f"dense materialization of {kinds[0]} input is deliberately "
            "unresolved -- consume the sparse tables directly "
            "or see https://github.com/hypertidy/controlledburn/issues/13"
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


# --------------------------------------------------------------------
# Arrow input -- consume WKB directly from the Arrow C data interface
# (GeoArrow "wkb" / plain (large_)binary) via nanoarrow. No shapely, and
# no per-geometry Python objects on the fast path (see burn()).
# --------------------------------------------------------------------

def _is_arrow(x) -> bool:
    return hasattr(x, "__arrow_c_array__") or hasattr(x, "__arrow_c_stream__")


def _is_shapely_like(x) -> bool:
    # objects as_wkb_list() already handles directly (shapely / geopandas);
    # these take precedence over the Arrow path even when they also expose
    # __arrow_c_stream__, so existing behaviour is preserved.
    return (
        hasattr(x, "geom_type")
        or hasattr(x, "to_wkb")
        or (hasattr(x, "geometry") and hasattr(x, "columns"))
    )


def _arrow_chunks(x):
    """Yield ``(values, offsets, valid, n)`` per chunk from Arrow input.

    ``values`` is a uint8 numpy view of the chunk's WKB data buffer,
    ``offsets`` an int64 array of length ``n + 1``, ``valid`` a uint8 mask
    (1 = valid) or ``None`` when there are no nulls. Imports the object
    through nanoarrow's Arrow C data interface support -- accepts
    ``(large_)binary`` and ``geoarrow.wkb`` (binary storage) arrays and
    streams.
    """
    import nanoarrow as na

    for chunk in na.Array(x).iter_chunks():
        n = len(chunk)
        if n == 0:
            continue
        bufs = chunk.buffers
        if len(bufs) != 3:
            raise NotImplementedError(
                "Arrow input must be (large_)binary or geoarrow.wkb WKB; "
                f"got an array with {len(bufs)} buffers (a native GeoArrow "
                "encoding?). Convert to WKB first."
            )
        valid_b, off_b, data_b = bufs
        off = np.asarray(off_b)
        if off.dtype not in (np.int32, np.int64):
            raise NotImplementedError(
                "Arrow input must be (large_)binary or geoarrow.wkb WKB "
                f"(int32/int64 offsets); got offset dtype {off.dtype}."
            )
        offsets = np.ascontiguousarray(off, dtype=np.int64)
        if offsets.shape[0] != n + 1:
            raise ValueError(
                "unexpected offsets length; sliced arrays with a non-zero "
                "offset are not supported yet"
            )
        data = np.asarray(data_b)
        if data.dtype != np.uint8:
            data = data.view(np.uint8)
        vb = np.asarray(valid_b)
        valid = None
        if vb.size:
            valid = np.unpackbits(
                vb.astype(np.uint8, copy=False), bitorder="little"
            )[:n].astype(np.uint8)
        yield data, offsets, valid, n


def _arrow_to_wkb_list(x) -> list:
    """Materialize Arrow WKB input to a ``list[bytes | None]``.

    The compatibility fallback used by :func:`as_wkb_list`; :func:`burn`
    takes the zero-copy buffer path instead (see :func:`_burn_arrow`).
    """
    out = []
    for data, offsets, valid, n in _arrow_chunks(x):
        for i in range(n):
            if valid is not None and not valid[i]:
                out.append(None)
            else:
                out.append(bytes(data[offsets[i]:offsets[i + 1]]))
    return out


def _arrow_extent(chunks):
    """Derive an ``(xmin, xmax, ymin, ymax)`` extent from Arrow WKB chunks.

    Uses the C++ core's WKB-envelope scan (``_core.bbox_wkb_arrow``) --
    no shapely. Chunks with no finite coordinate contribute nothing; an
    all-empty input raises.
    """
    xmin = ymin = float("inf")
    xmax = ymax = float("-inf")
    found = False
    for values, offsets, valid, _n in chunks:
        bb = _core.bbox_wkb_arrow(values, offsets, valid)
        if bb is None:
            continue
        bxmin, bymin, bxmax, bymax = bb
        xmin = min(xmin, bxmin)
        ymin = min(ymin, bymin)
        xmax = max(xmax, bxmax)
        ymax = max(ymax, bymax)
        found = True
    if not found:
        raise ValueError(
            "cannot derive extent from Arrow input: no finite coordinates"
        )
    return (xmin, xmax, ymin, ymax)


def _burn_arrow(geoms, extent, shape, resolution, mode) -> BurnResult:
    """Arrow-native burn: consume WKB directly from Arrow buffers.

    The zero-copy counterpart of :func:`burn` for Arrow input (anything
    exposing the Arrow C data interface with ``(large_)binary`` or
    ``geoarrow.wkb`` geometry). No shapely, no per-geometry Python
    objects -- spans point straight into the Arrow values buffer.

    When ``extent`` is None it is derived from the geometry via a
    WKB-envelope scan in the C++ core (:func:`_arrow_extent`), so the
    Arrow path has no shapely dependency at all.
    """
    # An Arrow stream is single-pass, and deriving the extent needs a look
    # at the data before the burn -- so materialize the chunk buffers once.
    chunks = list(_arrow_chunks(geoms))

    if extent is None:
        extent = _arrow_extent(chunks)
    extent, shape = resolve_grid(None, extent, shape, resolution)
    xmin, xmax, ymin, ymax = extent
    nrow, ncol = shape

    dtypes = {"runs": _RUNS_DTYPE, "edges": _EDGES_DTYPE,
              "lines": _LINES_DTYPE, "points": _POINTS_DTYPE}
    parts = {k: [] for k in dtypes}
    id_offset = 0
    for values, offsets, valid, n in chunks:
        raw = _core.burn_wkb_arrow(
            values, offsets, valid,
            xmin, ymin, xmax, ymax, ncol, nrow, mode,
        )
        for geom_index, message in raw["notes"]:
            warnings.warn(f"geometry {geom_index + id_offset}: {message}",
                          stacklevel=3)
        for key, dt in dtypes.items():
            arr = _structured(raw[key], dt)
            if id_offset and len(arr):
                arr = arr.copy()
                arr["id"] += id_offset
            parts[key].append(arr)
        id_offset += n

    return BurnResult(
        runs=np.concatenate(parts["runs"]) if parts["runs"]
        else np.empty(0, _RUNS_DTYPE),
        edges=np.concatenate(parts["edges"]) if parts["edges"]
        else np.empty(0, _EDGES_DTYPE),
        lines=np.concatenate(parts["lines"]) if parts["lines"]
        else np.empty(0, _LINES_DTYPE),
        points=np.concatenate(parts["points"]) if parts["points"]
        else np.empty(0, _POINTS_DTYPE),
        extent=extent,
        shape=shape,
    )
