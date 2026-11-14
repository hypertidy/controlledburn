# Experiments behind the GDAL plan

Each directory is self-contained: a script, how to run it, and the output it
gave on 2026-11-14 (`results.txt`). Versions:

- GDAL 3.14.0dev-9bb42f5a7b (local build), linked against GEOS 3.12.1, so
  `--pixels fractional` is not available in its zonal-stats;
- GEOS 3.14.2dev (branch `3.14`), 3.15.0 (release) and 3.16.0dev `2d12763`
  (main, 2026-11-13), built for experiments 01 and 02;
- controlledburn C++ core at this repository's `cpp/` (crate 0.3.1 parity).

None of this is part of the R package (`design/` is in `.Rbuildignore`).

| # | question | result |
|---|---|---|
| 01 | Does GEOS `GridIntersectionFractions` (exactextract in GEOS, used by `gdal raster zonal-stats --pixels fractional`) agree with controlledburn? How do they compare on speed and memory? | On the 193 valid polygons of both corpora: agree on all but **2 cells**. In both, GEOS gives 0 for a fully covered interior cell. Invalid (self-intersecting) polygons: GEOS returns fractions **outside [0, 1]**, from −0.98 to 2.0 (29 of 57 invalid polygons); controlledburn stays in [0, 1] (non-zero winding). Speed: controlledburn is **2.4× to 32× faster** on one large polygon, and GEOS needs a dense float buffer of the window (1.6 GB at 20000²). |
| 02 | Minimal reproductions of the GEOS problems from 01 | GEOS **segfaults** on valid polygons: 3.14.2 on `POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))` (100 × 100 grid over [0, 10]²), and **3.15.0** (the current release) on `gen140`. Both are fixed on main by #1589/#1593 (2026-11). A fully covered interior cell comes back 0 on all three (`random_39`), and a second one on main (`gen140`, which crashes the releases). GEOS's own intersection area is 1 for both. |
| 03 | GDAL on a rotated geotransform | `gdal_rasterize` is right: 606 cells, the truth. **`gdal raster zonal-stats` returns count 0**, with exit status 0 and no warning. The same polygon unrotated gives 606. |
| 04 | GDAL onto a curvilinear (geolocation array) target | `gdal_rasterize` **burns nothing**, even with `-to DST_METHOD=GEOLOC_ARRAY`: without a geotransform it forces `DST_METHOD=NO_GEOTRANSFORM`. Moving vertices only through GDAL's geoloc transformer gets **129 of 132** sliver cells wrong. Densifying to about one cell makes the cell-centre burn exact. |
| 05 | GDAL's rasterizer against controlledburn, both into the same dense Byte buffer | Cell-centre rule: **identical cells**. controlledburn is **100× faster** on one 100k-vertex polygon (GDAL loops over every edge on every row), and 1.5× (one thread) to 7× (eight) faster on 20k small polygons. For ALL_TOUCHED against coverage > 0, GDAL's extra cells have true areas below 1e-6 of a cell: controlledburn drops edge fractions ≤ 1e-6 by design. Exact coverage costs 5 to 10× the approx burn. |
