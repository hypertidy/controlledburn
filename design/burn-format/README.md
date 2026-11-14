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

## Open questions

- **The name.** `.burn` and a driver short name `BURN` need a collision
  check. `-burn` is already a `gdal_rasterize` option, which may confuse
  people.
- **Separate runs and edges, or one span table?** (See "What it stores".)
- **Id type.** i32 limits a file to about 2.1 billion features. Should it
  be i64, as OGR FIDs are?
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
