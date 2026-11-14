# 01: GEOS GridIntersectionFractions vs controlledburn

`compare.cpp` burns every polygon of two corpora with the controlledburn C++
core and with `GEOSGridIntersectionFractions_r`, the function
`gdal raster zonal-stats --pixels fractional` calls (`alg/zonal.cpp`,
`CalculateCoverage`). It compares them cell by cell, and against
exactextract 0.3.0 and GDAL's cell-centre rasterize where the corpus has them.

    GEOS_PREFIX=/path/to/geos-main ORACLE_DIR=~/Git/controlledburn-rs/fixtures/oracle ./run.sh

GEOS 3.14.x crashes on an oracle case (`snap_a`) and 3.15.0 on a canon case
(`gen140`, see 02), so this needs GEOS main (fixes of 2026-11-05/08). `results.txt` is from GEOS
3.16.0dev 2d12763. `differences.csv` lists every cell that differs by more
than 1e-6.

Reading the results:

- **Oracle corpus, coverage** (burnscars' 50 cases, all valid):
  - controlledburn and GEOS differ in 1 cell, `random_39` (36, 18). GEOS
    gives 2.8e-30; the GEOS intersection area of that cell is 1.0.
  - GEOS has fixed 4 of the 5 interior-cell zeros that exactextract 0.3.0
    shows; this one remains.
- **Oracle corpus, approx:**
  - the 64 differing cells are the documented tie convention.
    controlledburn's cell-centre rule is left-inclusive and half-open in y;
    GDAL rounds intercepts with `floor(x + 0.5)`.
  - Each such centre lies within 4e-14 of a cell of an edge (see
    controlledburn-rs `fixtures/oracle/known_diffs.toml`).
- **Canon corpus** (controlledburn's 240 generated cases, polygonal
  geometries only):
  - 142 cells differ, in 39 geometries. 38 of the 39 are invalid
    (self-intersecting, which the generator produces on purpose).
  - On those, GEOS returns fractions outside [0, 1], down to −0.98 and up
    to 2.0 (29 geometries out of range). One invalid geometry made GEOS return an error.
  - The one valid geometry that differs is `gen140#0`, cell (11, 4): GEOS
    gives 0 and the true fraction is 1.
- **Timing** (one wavy star polygon, best of 3):
  - GEOS's time includes the WKB parse and filling the dense buffer, as
    GDAL's zonal-stats does.
  - controlledburn's time is the sparse burn only.
