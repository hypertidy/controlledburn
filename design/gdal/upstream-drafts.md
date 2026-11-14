# Facts for upstream reports (nothing filed)

These are lists of facts and reproductions, not text to paste. Under GDAL's
LLM tool policy (RFC 111), the reporter writes the issue in their own words
and discloses that an LLM was used in the investigation. Each item names the
experiment that reproduces it. Check that each one still reproduces on
current master before filing.

## GDAL 1: zonal-stats ignores rotation terms of the geotransform

- Command: `gdal raster zonal-stats --stat count`, default pixels, on a
  40 × 30 raster with a 30° rotated geotransform and one polygon inside it.
- Result: count 0, exit 0, no warning. The unrotated equivalent gives 606.
  `gdal_rasterize` on the rotated raster burns 606 cells, the right count.
- Cause (code reading): `alg/zonal.cpp`, `CalculateCoverage` builds the
  coverage window's geotransform with `xrot = yrot = 0` from the geometry's
  envelope. The GEOS path takes an axis-aligned box.
- Fix options:
  - error out on rotation;
  - or transform the geometry to pixel space with the inverse geotransform
    and compute coverage on a unit grid (exact for affine).
- Reproduction: `experiments/03-rotated-gt/rotated.py`.
- Seen on GDAL 3.14.0dev-9bb42f5a7b with GEOS 3.12. The fractional path is
  untested here; it is likely the same, since it uses the envelope as an
  axis-aligned box.

## GDAL 2: rasterize cannot target a dataset georeferenced by geolocation arrays

- Target: a MEM dataset with only `GEOLOCATION` metadata (X/Y arrays), no
  geotransform.
- `gdal.Rasterize` burns nothing; with `-to DST_METHOD=GEOLOC_ARRAY`, still
  nothing.
- Cause:
  - `apps/gdal_rasterize_lib.cpp`, around "If we have transformer
    options", and `alg/gdalrasterize.cpp` (`GDALRasterizeGeometries`) both
    set `DST_METHOD=NO_GEOTRANSFORM` when the target has no geotransform,
    GCPs or RPCs.
  - That overrides the user's `DST_METHOD`.
- If honoured, moving vertices alone is still wrong on a non-affine grid:
  129 of 132 cells wrong for a long sliver on a polar sector grid. Densifying
  to about one cell gives the exact cell-centre result. A densify option is
  needed along with the fix.
- `gdalwarp` honours destination geolocation (3.5+), so this is an
  inconsistency, not a missing concept.
- Reproduction: `experiments/04-geoloc-rasterize/geoloc.py`.

## GEOS 1: crash in GridIntersectionFractions on 3.14 and 3.15.0

- `GEOSGridIntersectionFractions_r` crashes (SIGSEGV in
  `Cell::traversal_in_progress`) on valid polygons:
  - 3.14.2dev: `POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))`, extent
    0 0 10 10, 100 × 100;
  - 3.15.0: the `gen140` case in
    `experiments/02-geos-crash/interior-cases.txt` (hex WKB and grid there).
- Both are fine on main (after #1589 and #1593).
- Ask: backport to 3.14 and 3.15? GDAL 3.12+ calls this from
  `gdal raster zonal-stats --pixels fractional`.
- Reproduction: `experiments/02-geos-crash/repro.c`, `interior.c`.

## GEOS 2: fully covered interior cells return 0 (main 2d12763)

- `random_39`, cell (36, 18): `GridIntersectionFractions` 2.8e-30, while
  GEOSIntersection with the cell box gives fraction 1.
- `gen140`, cell (11, 4): 0 against 1.
- Both polygons are valid. These are the same symptom as the exactextract
  0.3.0 interior-cell bug documented in burnscars' oracle notes. 4 of the 5
  burnscars cases are fixed on main; `random_39` is not.
- Reproduction: `experiments/02-geos-crash/interior.c` with
  `interior-cases.txt`.

## GEOS 3 (discussion, not a bug report): fractions on invalid polygons

- Problem: on self-intersecting polygons, fractions come back outside
  [0, 1] (−0.98 to 2.0) on 29 of the 57 invalid polygons of
  controlledburn's canon corpus.
- Why it matters: zonal-stats uses these as weights.
- Question: undefined by contract (valid input required), or should it
  return non-zero-winding or even-odd area? At least document it;
  `GDALZonalStats` could check validity or clamp.
- Reproduction: `experiments/01-geos-vs-cb` (`differences.csv`, comparison
  `canon_cb_vs_geos`).
