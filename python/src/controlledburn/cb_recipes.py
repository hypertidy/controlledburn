"""cb_recipes -- helper layer for controlledburn's Python bindings.

Everything here is pure Python over the existing public API. Nothing
touches the C++ core. The intent is twofold:

  1. Make the Python surface as immediately usable as the R package's.
  2. Serve as a working prototype for what should eventually be
     upstreamed into controlledburn/__init__.py.

Usage:

    import controlledburn as cb
    from cb_recipes import burn, materialize_all, summary, transform

    r = burn(some_geopandas_geoseries, resolution=0.01)
    print(summary(r))
    arr = materialize_all(r)

Requires numpy. shapely is needed only for geometry coercion and for
deriving bounds from geometry; the rest works on WKB bytes alone.
"""

from __future__ import annotations

import numpy as np

import controlledburn as _cb

__all__ = [
    "as_wkb_list",
    "resolve_grid",
    "burn",
    "summary",
    "materialize_all",
    "covered_area",
    "transform",
    "cell_centers",
    "to_dataframe",
    "expand_runs",
    "to_coo",
]


# --------------------------------------------------------------------
# input coercion -- the Python analogue of R's as_wkb_list() S3 generic
# --------------------------------------------------------------------

def as_wkb_list(x):
    """Coerce assorted geometry containers to a list of WKB blobs.

    Accepts: raw bytes, a single shapely geometry, a shapely array, a
    geopandas GeoSeries or GeoDataFrame, or any sequence mixing WKB
    bytes, shapely geometries and None.
    """
    if isinstance(x, (bytes, bytearray, memoryview)):
        return [x]

    if hasattr(x, "geom_type"):                 # single shapely geometry
        import shapely
        return [shapely.to_wkb(x)]

    if hasattr(x, "geometry") and hasattr(x, "columns"):   # GeoDataFrame
        x = x.geometry

    if hasattr(x, "to_wkb"):                    # GeoSeries / GeometryArray
        return list(x.to_wkb())

    out = []
    shapely = None
    for g in x:
        if g is None or isinstance(g, (bytes, bytearray, memoryview)):
            out.append(g)
        elif hasattr(g, "geom_type"):
            if shapely is None:
                import shapely
            out.append(shapely.to_wkb(g))
        else:
            raise TypeError(
                f"cannot interpret {type(g).__name__} as geometry; "
                "pass WKB bytes, a shapely geometry, or None"
            )
    return out


# --------------------------------------------------------------------
# grid resolution -- the analogue of R's .resolve_grid_params()
# --------------------------------------------------------------------

def resolve_grid(geoms, bounds=None, shape=None, resolution=None, size=256):
    """Fill in bounds and shape the way the R package does.

    bounds defaults to the bounding box of `geoms`. shape defaults to a
    grid of at most `size` cells on the longer axis, preserving aspect
    ratio. `shape` and `resolution` are mutually exclusive.

    Returns (bounds, shape) with bounds as (xmin, ymin, xmax, ymax) and
    shape as (nrow, ncol).
    """
    if bounds is None:
        import shapely
        wkb = as_wkb_list(geoms)
        wkb = [w for w in wkb if w is not None]
        if not wkb:
            raise ValueError("cannot derive bounds from empty geometry")
        b = np.atleast_2d(shapely.bounds(shapely.from_wkb(wkb)))
        bounds = (float(np.nanmin(b[:, 0])), float(np.nanmin(b[:, 1])),
                  float(np.nanmax(b[:, 2])), float(np.nanmax(b[:, 3])))

    xmin, ymin, xmax, ymax = (float(v) for v in bounds)
    if not (xmax > xmin and ymax > ymin):
        raise ValueError("invalid bounds: need xmax > xmin and ymax > ymin")

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
        raise ValueError("shape must be positive")

    return (xmin, ymin, xmax, ymax), (nrow, ncol)


# --------------------------------------------------------------------
# burn wrapper -- optional bounds/shape, resolution, flexible input
# --------------------------------------------------------------------

def burn(geoms, bounds=None, shape=None, resolution=None, mode="coverage",
         size=256):
    """controlledburn.burn() with R-style conveniences.

    Unlike cb.burn(), `bounds` and `shape` are optional and `geoms` may
    be shapely geometries, a GeoSeries, a GeoDataFrame, or WKB.
    """
    wkb = as_wkb_list(geoms)
    bounds, shape = resolve_grid(wkb, bounds, shape, resolution, size)
    return _cb.burn(wkb, bounds=bounds, shape=shape, mode=mode)


# --------------------------------------------------------------------
# summary -- the analogue of R's print.controlledburn
# --------------------------------------------------------------------

def summary(r):
    """One-screen summary of a BurnResult. Safe on huge results."""
    nrow, ncol = r.shape
    interior = int((r.runs["col_end"] - r.runs["col_start"] + 1).sum())
    touched = interior + len(r.edges) + len(r.lines) + len(r.points)
    ids = np.unique(np.concatenate([
        r.runs["id"], r.edges["id"], r.lines["id"], r.points["id"]]))
    total = float(nrow) * float(ncol)

    n = len(ids)
    out = [f"<BurnResult> {ncol} x {nrow} grid, {n} "
           f"{'geometry' if n == 1 else 'geometries'}",
           f"  runs:   {len(r.runs)} ({interior} interior cells)",
           f"  edges:  {len(r.edges)} polygon boundary cells"]
    if len(r.lines):
        out.append(f"  lines:  {len(r.lines)} cells")
    if len(r.points):
        out.append(f"  points: {len(r.points)} cells")
    out.append(f"  sparsity: {100.0 * (1.0 - touched / total):.1f}% empty")
    return "\n".join(out)


# --------------------------------------------------------------------
# dense materialization covering all four tables
# --------------------------------------------------------------------

def _kinds_present(r):
    kinds = []
    if len(r.runs) or len(r.edges):
        kinds.append("polygon")
    if len(r.lines):
        kinds.append("line")
    if len(r.points):
        kinds.append("point")
    return kinds


def materialize_all(r, shape=None, background=0.0, allow_mixed=False):
    """Dense sum of runs, edges, lines and points.

    The analogue of R's materialise_chunk(). Unlike cb.materialize(),
    this does not drop lines and points.

    Raises on mixed geometry kinds unless allow_mixed=True, because the
    three measures have different units: edge fractions are
    dimensionless, line lengths are in CRS units, points are counts.
    """
    kinds = _kinds_present(r)
    if len(kinds) > 1 and not allow_mixed:
        raise ValueError(
            f"mixed geometry kinds present ({', '.join(kinds)}); the "
            "measures have different units. Filter by id to one kind, "
            "burn each kind separately, or pass allow_mixed=True."
        )

    nrow, ncol = shape if shape is not None else r.shape
    out = np.full((int(nrow), int(ncol)), float(background), dtype="f8")

    for run in r.runs:
        out[run["row"] - 1, run["col_start"] - 1:run["col_end"]] += 1.0
    if len(r.edges):
        np.add.at(out, (r.edges["row"] - 1, r.edges["col"] - 1),
                  r.edges["fraction"].astype("f8"))
    if len(r.lines):
        np.add.at(out, (r.lines["row"] - 1, r.lines["col"] - 1),
                  r.lines["length"].astype("f8"))
    if len(r.points):
        np.add.at(out, (r.points["row"] - 1, r.points["col"] - 1), 1.0)
    return out


def covered_area(r):
    """Total polygon area covered, in CRS units, from runs + edges."""
    xmin, ymin, xmax, ymax = r.bounds
    nrow, ncol = r.shape
    cell = ((xmax - xmin) / ncol) * ((ymax - ymin) / nrow)
    full = int((r.runs["col_end"] - r.runs["col_start"] + 1).sum())
    return cell * (full + float(r.edges["fraction"].sum(dtype="f8")))


# --------------------------------------------------------------------
# georeferencing
# --------------------------------------------------------------------

def transform(r):
    """GDAL/rasterio affine tuple (a, b, c, d, e, f) for the burn grid."""
    xmin, ymin, xmax, ymax = r.bounds
    nrow, ncol = r.shape
    return ((xmax - xmin) / ncol, 0.0, xmin,
            0.0, -(ymax - ymin) / nrow, ymax)


def cell_centers(row, col, bounds, shape):
    """1-based (row, col) to (x, y) cell centre coordinates."""
    xmin, ymin, xmax, ymax = bounds
    nrow, ncol = shape
    x = xmin + (np.asarray(col) - 0.5) * (xmax - xmin) / ncol
    y = ymax - (np.asarray(row) - 0.5) * (ymax - ymin) / nrow
    return x, y


# --------------------------------------------------------------------
# tabular and sparse interop
# --------------------------------------------------------------------

def expand_runs(runs):
    """Expand run-length encoded interior cells to one row per cell.

    Returns a structured array with columns row, col, id.
    """
    if len(runs) == 0:
        return np.empty(0, dtype=[("row", "i4"), ("col", "i4"), ("id", "i4")])
    n = runs["col_end"] - runs["col_start"] + 1
    cols = np.concatenate([np.arange(a, b + 1) for a, b in
                           zip(runs["col_start"], runs["col_end"])])
    out = np.empty(int(n.sum()),
                   dtype=[("row", "i4"), ("col", "i4"), ("id", "i4")])
    out["row"] = np.repeat(runs["row"], n)
    out["col"] = cols
    out["id"] = np.repeat(runs["id"], n)
    return out


def to_dataframe(r, table="edges", expand=False):
    """Convert one table to a pandas DataFrame.

    table: 'runs', 'edges', 'lines' or 'points'.
    expand: for 'runs', expand to one row per cell.
    """
    import pandas as pd
    arr = getattr(r, table)
    if table == "runs" and expand:
        arr = expand_runs(arr)
    return pd.DataFrame(arr)


def to_coo(r, kind="polygon"):
    """Convert one geometry kind to a scipy.sparse.coo_array.

    Polygons give coverage (1.0 for interior, fraction for boundary),
    lines give length-in-cell, points give counts.
    """
    from scipy.sparse import coo_array

    nrow, ncol = r.shape
    rows, cols, vals = [], [], []

    if kind == "polygon":
        if len(r.runs):
            ex = expand_runs(r.runs)
            rows.append(ex["row"] - 1)
            cols.append(ex["col"] - 1)
            vals.append(np.ones(len(ex)))
        if len(r.edges):
            rows.append(r.edges["row"] - 1)
            cols.append(r.edges["col"] - 1)
            vals.append(r.edges["fraction"].astype("f8"))
    elif kind == "line":
        rows.append(r.lines["row"] - 1)
        cols.append(r.lines["col"] - 1)
        vals.append(r.lines["length"].astype("f8"))
    elif kind == "point":
        rows.append(r.points["row"] - 1)
        cols.append(r.points["col"] - 1)
        vals.append(np.ones(len(r.points)))
    else:
        raise ValueError("kind must be 'polygon', 'line' or 'point'")

    if not rows or sum(len(v) for v in vals) == 0:
        return coo_array((nrow, ncol))

    return coo_array(
        (np.concatenate(vals),
         (np.concatenate(rows), np.concatenate(cols))),
        shape=(nrow, ncol),
    )
