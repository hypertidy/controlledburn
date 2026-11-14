# 05: GDAL's rasterizer vs controlledburn

`./run.sh` (needs `gdal-config`) burns two corpora on a 10000 × 10000 grid of
unit cells into a dense Byte buffer, with GDAL's `GDALRasterizeGeometries`
and with the controlledburn C++ core (burn, then materialise). Best of 3.

- **Cell-centre rule:**
  - identical cells on both corpora;
  - one 100k-vertex polygon: GDAL 1505 ms, controlledburn 15 ms;
  - 20k small polygons: GDAL 260 to 325 ms, controlledburn 167 ms on one
    thread and 44 to 47 ms on eight.
  - Why the big gap: `GDALdllImageFilledPolygon` (`alg/llrasterize.cpp`)
    loops over every edge of the geometry on every scanline, which is
    O(rows × edges). That is the gap fasterize showed in OSGeo/gdal#7200.
- **ALL_TOUCHED against controlledburn fraction > 0:**
  - 134 and 1607 cells are set only by GDAL. All of the 535 checked have a
    true intersection area above 0 and below 1e-6 of a cell.
  - controlledburn drops edge fractions ≤ 1e-6 (`cpp/src/scanline.cpp`,
    the `tol` in the edge sweep), so coverage > 0 is not an exact
    ALL_TOUCHED. A GDAL design needs an explicit threshold or an exact
    predicate (`relate_cells`).
  - 1 to 3 cells go the other way (set only by controlledburn, with
    area > 0), consistent with the ALL_TOUCHED edge cases in
    OSGeo/gdal#10606 and #14615.
- **Cost of exact coverage:** 5 to 10 times the approx burn. On many small
  polygons it is slower than GDAL's ALL_TOUCHED on one thread (1.7 s
  against 0.4 s), but not on eight (0.3 s).
