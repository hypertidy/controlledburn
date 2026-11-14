# 04: rasterize onto a geolocation-array (curvilinear) target

`python3 geoloc.py` builds a polar sector grid: 60 columns in angle, 40 rows in
radius. It is given only as 2D X/Y arrays of cell centres in `GEOLOCATION`
metadata, with `GEOREFERENCING_CONVENTION=PIXEL_CENTER`. Two polygons are
burned onto it: a long sliver and a quadrilateral. The truth is the cell-centre
rule: centre `(X[j,i], Y[j,i])` inside the polygon.

| method | sliver: missed / extra (of 132) | quad: missed / extra (of 676) |
|---|---|---|
| `gdal.Rasterize` onto the geolocated dataset | 132 / 0 (nothing burned) | 676 / 0 |
| same, with `-to DST_METHOD=GEOLOC_ARRAY` | 132 / 0 | 676 / 0 |
| vertices through GDAL's geoloc transformer, pixel-space burn | 129 / 100 | 114 / 313 |
| densified to 5 units first | 1 / 2 | 0 / 2 |
| densified to 1 unit (about a cell) | 0 / 0 | 0 / 0 |

Why nothing is burned:
- `apps/gdal_rasterize_lib.cpp` (`-to` path) and `alg/gdalrasterize.cpp`
  (default path) both set `DST_METHOD=NO_GEOTRANSFORM` whenever the target
  has no geotransform, GCPs or RPCs. That overrides a user's
  `DST_METHOD=GEOLOC_ARRAY` and ignores `GEOLOCATION` metadata.
- `gdalwarp` honours destination geolocation (GDAL ≥ 3.5); rasterize does not.

The vertices-only row is what GDAL would do if the override were removed: a
straight world edge is a curve in pixel space. Densifying in world units is a
stand-in for cutting edges in index space (controlledburn-rs
`design/grid-frontends.md`).
