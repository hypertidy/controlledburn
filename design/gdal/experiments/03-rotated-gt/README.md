# 03: rotated geotransform

`python3 rotated.py` (GDAL Python bindings and the `gdal` CLI) builds a 40 × 30
grid of unit cells rotated 30 degrees and one polygon, then counts its cells:

- **truth:** cell centres inside the polygon, by OGR, from world
  coordinates. 606.
- **`gdal.Rasterize`:** 606. GDAL moves vertices through the inverse
  geotransform, which is exact for an affine map.
- **`gdal raster zonal-stats --stat count`:** 0, exit status 0, no warning.
- **control:** the same polygon on the unrotated grid gives 606.

The cause, from reading `alg/zonal.cpp`:
- the per-geometry coverage window is built from the geometry's envelope
  with `xrot = yrot = 0` (in `CalculateCoverage`, both paths);
- the GEOS path takes an axis-aligned box only.

So rotation is dropped rather than refused. The fractional path (GEOS ≥ 3.14)
is untested here because the local GDAL links GEOS 3.12.
