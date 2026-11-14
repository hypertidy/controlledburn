# 02: GEOS GridIntersectionFractions reproductions

- `repro.c`: SIGSEGV in `geos::operation::grid::Cell::traversal_in_progress`.
  - Input: `POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))`, valid, on a
    100 × 100 grid over [0, 10]². Grids of 10, 20, 50, 99, 101 and 1000
    cells do not crash.
  - Crashes on GEOS 3.14.2dev (branch `3.14`, 296b9c1); fine on 3.15.0
    and on main 2d12763.
- `interior.c` with `interior-cases.txt`: a cell fully inside a valid
  polygon comes back 0 (`random_39` cell (36, 18): 2.8e-30; `gen140` cell
  (11, 4): 0).
  - GEOSIntersection of the cell box gives a fraction of 1 in both.
  - Still present on main.
  - The `gen140` input also **crashes GEOS 3.15.0** (the current
    release) and 3.14.2. It is fine on main, after libgeos/geos#1589 and
    #1593 (merged 2026-11-05 and 11-08, after 3.15.0 of 2026-09-01).
  - GDAL 3.12+ built against GEOS 3.14 or 3.15 is exposed through
    `gdal raster zonal-stats --pixels fractional`.

Build:

    cc repro.c -I$GEOS/include -L$GEOS/lib -lgeos_c -o repro
    cc interior.c -I$GEOS/include -L$GEOS/lib -lgeos_c -o interior
    while read name args; do ./interior $args; done < interior-cases.txt

`results.txt` is the run on all three versions (in it, a run that stops after
`valid: 1` is the crash). Draft report text is in `../../upstream-drafts.md`. Nothing has been filed.
