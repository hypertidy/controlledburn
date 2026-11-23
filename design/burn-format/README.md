# BURN: the Rasterized Vector Format (working notes)

Status: notes towards a spec, 2026-11-14. Nothing here is decided, and
the name is a placeholder. These notes were written with Claude (an LLM).
A GDAL driver in GDAL's own tree must be written by a person (RFC 111;
see "GDAL drivers" below).

## What it is

A file holding one burn: polygons, lines and points burned onto a grid,
kept as controlledburn's sparse tables. It has two readings:

- **as a raster:** the grid is implicit, and any cell value (coverage,
  feature id, count, an attribute) is computed from the tables when it is
  read;
- **as vector:** the features, each with the set of cells it covers and
  the measures of those cells (area fraction, line length, point count).

The tables are the data. Everything else, such as bounding boxes,
per-feature measures and spatial indexes, can be derived from them in one
pass. The file may cache these ("latent indexes"), but a reader must not
depend on the cache being there.

The file is meant to be ready for:

- **zonal statistics:** join the cells to a value raster on the same grid
  by (row, col), then group by id with fraction as the weight. This costs
  O(boundary), with no dense window per feature.
- **DE-9IM per cell:** the crate's `relate_cells` needs the burn *and the
  source geometries*. So the format carries the geometries, or a
  reference to them, when this reading is wanted.

## What it stores

Every record comes from the crate's four-table contract (`src/output.rs`):

- indices are 0-based, row 0 is the top row, and `col_end` is exclusive;
- each table's measure column means exactly one thing.

| table | columns | measure |
|---|---|---|
| `runs` | row, col_start, col_end, id (i32) | none: fully covered |
| `edges` | row, col, fraction (f32), id | area fraction in (0, 1) |
| `lines` | row, col, length (f32), id | length in CRS units |
| `points` | row, col, id | none: count 1 |
| `notes` | geom_index, message | the burn's notes on skipped input |

- **Decide: keep runs and edges as separate tables, or one `polygon`
  table of spans with a fraction?** A single table would be
  (row, col_start, col_end, fraction, id): a run is a span with fraction
  1, and an edge is a span of length 1. That is the shape plan.md
  Stage 4 sketches for `gdal vector rasterize` table output, and it makes
  materialising a single loop. Separate tables keep the crate's schema
  exact and keep the measure columns type-pure.
- Fractions are f32, good to about 6e-8 of a cell. This matches the
  crate's records; the crate has `CoverageWindow` for f64.

### Burn metadata (required)

Without these a reader cannot interpret the tables:

- the crate version;
- the burn options: mode (coverage or approx), merge_parts,
  coalesce_runs and id_offset;
- the fill rule (non-zero winding);
- the edge threshold: fractions of 1e-6 or less are dropped. Record it,
  because "fraction > 0" is then not an exact ALL_TOUCHED
  (`design/gdal/experiments/05-rasterize-speed`);
- the vertex snap tolerance (1e-9 of a cell).

## The grid: the implicit raster

The burn runs in index space, so the tables do not depend on the kind of
grid. Only the grid's description does (crate `design/grid-frontends.md`).

| grid kind | description | exact? |
|---|---|---|
| regular | ncol, nrow, extent | yes |
| affine | ncol, nrow, a six-term geotransform (rotation and shear) | yes, vertices only |
| rectilinear | x and y cell-edge arrays | yes, once edges are cut at cell walls |
| curvilinear | X and Y arrays of cell corners or centres, given inline or by reference (as GDAL's GEOLOCATION metadata) | to within about 1% of a cell for bilinear cells; exact with the equal-area split |

- A CRS (WKT2 or PROJJSON) goes with every grid kind.
- **Use GTI's layer-metadata names where they fit**
  (`doc/source/drivers/raster/gti.rst` in GDAL's source): `XSIZE`, `YSIZE`,
  `GEOTRANSFORM`, `SRS`, `RESX`, `RESY`. GTI is GDAL's existing "raster
  described by a vector table" driver, so GDAL users already know these
  names.
- **For a non-affine grid, record how the geometry was moved into index
  space.** This means the densify step, and whether bilinear or the
  equal-area split was used. Two burns of the same features onto the same
  curvilinear grid are only comparable if these match.

## Order and latent indexes

- **The burn's output is id-major:** geometry by geometry, then row by
  row. Writing in that order costs nothing, and a per-feature read is one
  contiguous range.
- **Raster reads want row-major order**, partitioned into tiles. The
  crate's `slice_tiles` does exactly that partitioning, and its tile
  order is deterministic.
- **Proposal: store tiles as Parquet row groups, sorted by (row,
  col_start, id) inside each tile.**
  - Row-group min/max statistics on row and col then *are* the spatial
    index, and they come for free.
  - A reader of one window or block reads only the row groups that
    overlap it.
  - **The id column keeps the input order.** So merge rules that depend
    on order (FIRST, LAST) still work after re-sorting.
- **A `features` table, one row per id,** derived in one pass:
  - the cell bounding box (row0, row1, col0, col1);
  - the measures: covered cells (runs plus summed fractions), line
    length, point count, and the record count.
  - This is the per-feature index, and it is what "bboxes and measures"
    means here.
- **Sparse awareness in GDAL:** a raster driver can answer
  `GetDataCoverageStatus` (`gcore/gdal.h`) from row-group statistics. Then
  tools that skip empty blocks (COG writing with `SPARSE_OK`, `gdalwarp`)
  can skip the empty parts of a BURN raster.

## The source features: the implicit vector

The `features` table can optionally carry:

- the source geometry as WKB;
- the source attributes;
- the source FID next to the burn id.

With these, the file is GeoParquet-shaped (GeoParquet is the
geospatial-metadata convention for Parquet), and existing readers see
it as a vector layer. They are also what makes the file DE-9IM-ready:
`relate_cells` can run later without going back to the source.

- **Optional `relate` table:** the output of `relate_cells`, as (row, col,
  id, matrix). The matrix can be a 9-character DE-9IM string, or packed
  into 18 bits.
- **Keep the relate table optional.** Most of it can be derived from the
  fraction alone (crate `design/de9im-cells.md`). Only the
  boundary-boundary entry needs the geometry.

## Burns from rasters

Everything above assumes vector input. A burn can also be made from a
raster, with no polygons in between. (Added 2026-11-23.)

### Isobands and isolines

- **A band [a, b) of a raster is a burn.** On the dual grid (below),
  each cell's share of the band depends only on its four corner values:
  - a cell whose four corners are all in the band is a run, since a
    bilinear surface takes its extremes over a cell at its corners;
  - a cell crossed by the isoline at a or b is an edge, and its fraction
    is the area in the band;
  - the length of the isoline at each level within a cell is a `lines`
    record.
  - Ids are band or level indexes, not features.
- **No ring building.** Polygon contour output joins segments into lines
  and nests rings into polygons across the whole raster (GDAL's
  `alg/marching_squares/`: `segment_merger.h` and
  `polygon_ring_appender.h`). A burn needs neither step. Tiles are
  independent, so a burn can be made in parallel or streamed, with
  nothing to stitch at tile seams.
- **Connected parts are optional.** Splitting a band into its connected
  parts, the features a polygon output would hold, is a union-find over
  the runs and edges. It is a separate step and can come later.
- **Bands partition the cell.** When the levels span the raster's values,
  the fractions of all bands in a cell sum to 1, to within the edge
  threshold.
- **Decide: a bilinear or a linear interpolant.**
  - Bilinear: the fraction is an area bounded by a hyperbola, which has a
    closed form. Saddle cells resolve themselves, and isolines of
    different levels never cross, because they are level sets of one
    surface.
  - Linear, as in marching squares (straight segments between points
    interpolated on the cell walls): cheaper, and it matches
    `gdal_contour`. But a saddle cell has two readings. GDAL always picks
    the same one (`alg/marching_squares/square.h`), because choosing by
    the centre value can make contours of different levels cross.
  - Record which was used, as the grid section asks for the densify step.
- **Decide: nodata.** A cell with a nodata corner has no bilinear
  surface. Either drop the cell, or split it into quarters around the
  mean of its valid corners, as GDAL's marching squares does
  (`square.h`).

### The grid

- **Isobands sit on the dual grid.** Marching squares treats pixel
  centres as cell corners, so the burn's grid is (nrow − 1) × (ncol − 1),
  offset by half a cell. For a regular or affine grid it is a shifted
  geotransform.
- **A curvilinear grid given by its centres already gives the dual grid's
  corners,** so the burn needs no new coordinates.
- **The isolines are defined in index space,** by interpolating between
  cell centres. So the burn has none of the index-space error that a
  vector burn onto a curvilinear grid has. The index-space burn *is* the
  definition, and the world-space view is the GEOLOCATION reading. Ocean
  model output on a tripolar grid is the case to try.
- **Decide: also offer the source grid itself?** That needs values at the
  source cells' corners, so an interpolant to the corners would have to
  be chosen first. The dual grid avoids that choice.

### Order and pyramids

- **The output comes in row order,** tile by tile. That is the order a
  raster read wants (see "Order and latent indexes"), so the tiled layout
  is the write order. The `features` table, one row per level, is the
  part that needs a pass.
- **For integer factors, coarser grids are a group-by.** Lengths in
  `lines` add up exactly. A coarse cell's fraction is the mean of its
  fine cells' fractions (a run counts as 1), to within the edge
  threshold. So the pyramid levels are
  `group by (row %/% f, col %/% f, id)` over one fine burn. Contouring
  each overview of the raster would give different contours, because an
  overview smooths the surface.
  - This gives the open question of several grids in one file a concrete
    case.

### Tanaka contours

- Tanaka's illuminated contours draw each contour lit or shaded, with a
  width proportional to |cos θ|, where θ is the angle between the
  contour's normal (the surface gradient) and the light.
- **At the native grid this needs no new measure.** The gradient comes
  from the source raster, joined by (row, col): for a dual cell, the mean
  gradient of its four corners. The weight to draw with is the length
  times cos θ.
- **A pyramid needs the weights summed at the fine level,** because θ
  cannot be recovered from a coarse cell. Keep the lit sum (cos θ > 0)
  and the shaded sum apart. Each is a sum of length × |cos θ|, so it adds
  up down the pyramid as lengths do, and a contour that bends within a
  coarse cell does not cancel itself out.
  - **Decide: store the two sums as tables of their own** beside `lines`
    (same columns, with a weighted length as the measure), so that each
    table keeps one measure.
- **The exact version** integrates cos θ along the isoline within the
  cell, using the bilinear surface's gradient rather than the cell mean.
  It is still local to the cell, so it changes the weight, not the
  layout.
- **The fill comes free.** The band burns give the hypsometric tint, and
  at the native grid the `edges` fractions anti-alias it.

### Other raster sources

- **Classified rasters and masks** (land cover, fire scars, flood
  extents, cloud masks) are `runs` only, with id = class. A row-wise
  run-length encoding of the raster is already the burn.
  - Cells with no data have no record, so `GetDataCoverageStatus`
    reports the empty blocks for free.
  - Such a burn is also a set of zones for zonal statistics.
  - Its boundaries lie on cell walls, which is the good case for the
    DE-9IM cell tables (crate `design/de9im-cells.md`).
- **Distance bands** are isobands of a distance transform, which makes
  them buffers in grid form.
- **Source cells as features.** Burn each cell of one grid (a square, or
  a quad on a curvilinear grid) onto another, with id = source cell. A
  target cell's fractions over the source cells are then its first-order
  conservative regridding weights. This is the ESMF case in crate
  `design/topology-mesh.md`.
  - The ids are source cells, so they overflow i32 on a grid of more than
    about 46,000 × 46,000 cells. This bears on the id-type question
    below.

## Container

- **A. A directory, or a zip of a directory, holding one Parquet file
  per table plus a JSON manifest** for the burn metadata and the grid.
  This is close to Zarr's layout.
  - As a zip it is a single `.burn` file.
  - GDAL already reads each table through `/vsizip/x.burn/edges.parquet`
    with the Parquet driver.
  - The manifest can be versioned on its own.
- **B. One Parquet file of the single span table,** with burn and grid
  metadata as Parquet key-value metadata and the features in a second
  file. This is the simplest to write, but it holds only one table, so
  lines and points need their own files or a kind column.
- **C. GeoPackage** (SQLite): several tables and an R-tree come
  natively. But it is row-oriented and larger, and it has no row-group
  statistics.

**Leaning A.** It keeps each table type-pure, and it is built entirely
from formats GDAL, Arrow, DuckDB, R and Python already read. GDAL's Arrow
stream API (RFC 86) then hands the tables to numpy or Arrow without a
copy.

## GDAL drivers

- **Raster:** a read-only driver.
  - It opens a `.burn` file and exposes bands computed on read:
    - coverage, the summed fraction (Float32);
    - id, by a merge rule: FIRST, LAST, or largest fraction;
    - count;
    - an attribute from `features`.
  - Open options would choose the band, the merge rule, and a threshold
    for "touched".
  - It reads blocks by row group and reports empty blocks through
    `GetDataCoverageStatus`.
  - An affine grid gives a geotransform. A curvilinear grid gives
    GEOLOCATION metadata, which `gdalwarp` honours already.
- **Vector:** the existing Parquet driver already reads the tables and
  features through `/vsizip/`, so a vector driver is optional. A combined
  BURN driver would mainly add a `cells` layer, with each record drawn as
  its cell polygon.
- **Write:** this is plan.md Stage 4: `gdal vector rasterize --output
  x.burn`, plus a materialise step. The format and the API should be
  designed together.
- **Routes, given RFC 111:**
  - **an out-of-tree C++ plugin driver** (loaded through
    `GDAL_DRIVER_PATH`; `doc/source/development/building_from_source.rst`,
    "Build drivers as plugins"). It is not a contribution to GDAL, so
    RFC 111 does not apply. It can use controlledburn's core, and it is
    the way to prototype and gather users;
  - **an RFC 76 Python driver:** vector-only and read-only, so it covers
    only the vector reading. A quick demo, not the raster driver;
  - **an in-tree GDAL driver:** written by a person, with the spec as the
    agreed reference. The spec should therefore live as a standalone,
    human-written document, outside both GDAL and controlledburn's code.

## Prior art to compare against

- **GTI:** a raster defined by a vector table. It is the nearest GDAL
  precedent.
- **GeoParquet and GeoArrow:** the vector side.
- **Run-length and compressed-sparse-row layouts** (scipy.sparse,
  AnnData): runs are row-wise run-length encoding.
- **Raster-interval approximations for spatial joins** (APRIL, cited in
  crate `design/de9im-cells.md`): the FULL and PARTIAL lists are runs and
  edges.
- **Zarr sparse proposals,** and xarray and rasterix as raster consumers.
- **exactextract's** per-feature output tables (cell, fraction).
- **GDAL's contouring** (`alg/marching_squares/`, and
  `gdal raster contour --polygonize`): the line and polygon outputs that
  a band burn would stand in for.
- **The R isoband package:** isolines and isobands from a grid.
- **Tanaka's illuminated contours** (1950), and hypsometric tinting.

## Open questions

- **The name.** `.burn` and a driver short name `BURN` need a collision
  check. `-burn` is already a `gdal_rasterize` option, which may confuse
  people.
- **Separate runs and edges, or one span table?** (See "What it stores".)
- **Id type.** i32 limits a file to about 2.1 billion features. Should it
  be i64, as OGR FIDs are? Regridding burns, whose ids are source cells,
  reach that limit on large grids.
- **Tile size for row groups.** A fixed default, or chosen from the
  record density?
- **Several burns in one file?** For example coverage and approx of the
  same features, or several grids (pyramid levels from one set of
  features).
- **Notes.** How should a burn with notes (skipped geometries) be
  flagged at open, so a reader is not misled by missing features?
- **Invalid input.** Coverage follows the winding rule, and GDAL's
  rasterize uses even-odd. Should a file record which one produced it,
  and should the raster driver refuse to mix them?
- **Updates.** Should appending features (with id_offset) be allowed
  without rewriting the file? The crate's id_offset makes appended burns
  concatenate.
