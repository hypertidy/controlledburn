# GDAL: planning notes

- [plan.md](plan.md): the staged plan for fast (approx) and exact (coverage)
  rasterization in GDAL, including affine, curvilinear and multidim grids,
  with constraints, routes and other implications.
- [experiments/](experiments/README.md): five experiments with scripts and
  results. They compare GEOS and controlledburn, test GDAL on rotated and
  geolocated grids, and benchmark GDAL's rasterizer.
- [upstream-drafts.md](upstream-drafts.md): facts for bug reports to GDAL and
  GEOS. Nothing has been filed.

The plan and drafts were written with Claude (an LLM). GDAL's policy (RFC 111)
requires human authorship of anything submitted, so treat these as working
notes. The repository's `design/` directory is excluded from the R package
build.
