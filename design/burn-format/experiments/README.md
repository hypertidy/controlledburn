# Experiments behind the BURN notes

Each directory is self-contained: scripts, how to run them, and the output
they gave. None of this is part of the R package (`design/` is in
`.Rbuildignore`).

| # | question | result |
|---|---|---|
| 01 | Can a raster's isobands and isolines be burned cell by cell, with no polygons? Do tiles, pyramids and Tanaka contours work from the tables? | Yes. On REMA 100 m around Mt Erebus (1M cells, 75 levels), the cell-local burn equals the burn of **GDAL's whole-raster contour polygons to f32**, and isoband's with the mean saddle rule. Tiles reproduce the whole burn **exactly**. The pyramid by group-by equals a burn onto the coarse grid, for factors 2, 4 and 8. Tanaka renders need no stroking. At 67M cells the tiles run in parallel with no merge step: 10.5 s for the R prototype on 16 cores, against 23.6 s for GDAL's whole-raster polygons on one thread. The linear interpolant is **discontinuous at saddles** (one cell jumped by 0.42), which argues for bilinear. The crate's threshold is symmetric and applied to f32. |
