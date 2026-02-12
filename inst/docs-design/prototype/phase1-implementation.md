# Phase 1 Implementation: gridburn (working name)

## The Key Finding

`GEOSGridIntersectionFractions_r` in GEOS 3.14 has this signature:

```c
int GEOSGridIntersectionFractions_r(
    GEOSContextHandle_t handle,
    const GEOSGeometry* g,
    double xmin, double ymin,
    double xmax, double ymax,
    unsigned nx, unsigned ny,
    float* buf);
```

It fills a **dense** `float[nx * ny]` buffer with coverage fractions (0.0–1.0).
Returns per-geometry. Caller manages the buffer.

## What our package actually does

The GEOS function gives us the raw computation. Our package's value is:

1. **Bounding-box clipping**: For each geometry, compute its bbox intersection with the grid, call `GEOSGridIntersectionFractions` only on that subgrid (not the full raster). This is critical — a small polygon on a 100k × 100k grid only needs a tiny buffer.

2. **Dense → sparse conversion**: Walk the dense buffer and produce the two-table output:
   - **runs table**: consecutive cells with coverage ≈ 1.0 → `(row, col_start, col_end, id)`
   - **edges table**: cells with 0 < coverage < 1.0 → `(row, col, weight, id)`

3. **Multi-geometry loop**: Iterate over a vector of geometries, accumulate runs and edges tables.

4. **Approximate mode**: Skip the coverage computation entirely, use the fasterize/controlledburn scanline algorithm → runs table only. (Or: call the exact function and threshold at weight > 0.)

## Architecture: Two backends

### Backend A: GEOS 3.14+ (future, when libgeos catches up)

```
R: geos_geometry vector + extent + dimension
  → C++ (cpp11 + libgeos):
      for each geometry:
        1. Compute bbox intersection with grid
        2. Allocate float buffer for subgrid
        3. Call GEOSGridIntersectionFractions_r()
        4. Walk buffer → append to runs/edges vectors
      Return list(runs = data.frame, edges = data.frame)
```

No vendored code. Just libgeos linkage.

### Backend B: Vendored exactextract (immediate, for GEOS < 3.14)

Same R-level API. The C++ layer calls the exactextract `RasterCellIntersection` 
class instead of `GEOSGridIntersectionFractions`. Requires vendoring ~10 files
from exactextract src/.

### Backend C: Pure R prototype (for design validation)

We can actually prototype the entire flow in R right now using exactextractr:

```r
library(exactextractr)
library(terra)

# Create a minimal raster matching our grid spec
r <- rast(nrows = ny, ncols = nx, 
          xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
values(r) <- 1  # dummy values

# Get per-cell coverage fractions
result <- exact_extract(r, polygons, fun = NULL, 
                        include_cols = NULL,
                        coverage_area = FALSE)
# result is a list of data.frames with columns: value, coverage_fraction
# (plus row/col if we use include_xy = TRUE or include_cell = TRUE)
```

Actually, `exact_extract(r, polygons, fun = NULL, include_cell = TRUE)` returns
data.frames with `cell` and `coverage_fraction` columns. We can convert cell → 
(row, col) using vaster, then RLE-compress the weight-1 cells.

This is **perfect for prototyping** — we can validate the output format, test the
two-table design, and benchmark before writing any C++.

## Immediate Plan: R Prototype

### Step 1: Pure R function using exactextractr

```r
burn_sparse <- function(geometry, extent, dimension, exact = TRUE) {
  # geometry: sf/sfc polygon(s)
  # extent: c(xmin, xmax, ymin, ymax)  
  # dimension: c(ncol, nrow)
  
  # Create minimal raster
  r <- terra::rast(nrows = dimension[2], ncols = dimension[1],
                   xmin = extent[1], xmax = extent[2], 
                   ymin = extent[3], ymax = extent[4])
  terra::values(r) <- 1L
  
  if (exact) {
    # Get coverage fractions per cell per geometry
    results <- exactextractr::exact_extract(r, geometry, fun = NULL,
                                            include_cell = TRUE)
    # Convert to two-table format
    # ... (see Step 2)
  } else {
    # Approximate mode — use controlledburn or fasterize
    # ... 
  }
}
```

### Step 2: Dense → sparse conversion in R

```r
dense_to_sparse <- function(cells_df, ncol, id) {
  # cells_df has columns: cell, coverage_fraction
  # Convert cell index to (row, col)
  row <- ((cells_df$cell - 1) %/% ncol) + 1L
  col <- ((cells_df$cell - 1) %% ncol) + 1L
  weight <- cells_df$coverage_fraction
  
  # Split into interior (weight ≈ 1) and edge (weight < 1) cells
  is_interior <- weight >= (1.0 - 1e-6)
  
  # Edge cells
  edges <- data.frame(
    row = row[!is_interior],
    col = col[!is_interior],
    weight = weight[!is_interior],
    id = id
  )
  
  # Interior cells → RLE compress by row
  int_row <- row[is_interior]
  int_col <- col[is_interior]
  
  if (length(int_row) == 0) {
    runs <- data.frame(row = integer(), col_start = integer(), 
                       col_end = integer(), id = integer())
  } else {
    # Sort by row then col, find runs
    ord <- order(int_row, int_col)
    int_row <- int_row[ord]
    int_col <- int_col[ord]
    
    # RLE: detect breaks in consecutive columns within same row
    row_break <- c(TRUE, diff(int_row) != 0)
    col_break <- c(TRUE, diff(int_col) != 1)
    run_start <- which(row_break | col_break)
    run_end <- c(run_start[-1] - 1, length(int_row))
    
    runs <- data.frame(
      row = int_row[run_start],
      col_start = int_col[run_start],
      col_end = int_col[run_end],
      id = id
    )
  }
  
  list(runs = runs, edges = edges)
}
```

### Step 3: Multi-geometry wrapper

```r
burn_sparse_all <- function(geometry, extent, dimension, exact = TRUE) {
  r <- terra::rast(nrows = dimension[2], ncols = dimension[1],
                   xmin = extent[1], xmax = extent[2],
                   ymin = extent[3], ymax = extent[4])
  terra::values(r) <- 1L
  
  results <- exactextractr::exact_extract(r, geometry, fun = NULL,
                                          include_cell = TRUE)
  
  all_runs <- vector("list", length(results))
  all_edges <- vector("list", length(results))
  
  for (i in seq_along(results)) {
    sparse <- dense_to_sparse(results[[i]], dimension[1], i)
    all_runs[[i]] <- sparse$runs
    all_edges[[i]] <- sparse$edges
  }
  
  list(
    runs = do.call(rbind, all_runs),
    edges = do.call(rbind, all_edges)
  )
}
```

### Step 4: Validate against controlledburn

```r
# Compare approximate mode output (runs only) against controlledburn
library(controlledburn)
cb_result <- controlledburn::rasterize(polygons, raster_spec)
# Compare row/col_start/col_end/id
```

## Package Structure (Eventual)

```
gridburn/
├── DESCRIPTION
│   Imports: geos
│   LinkingTo: cpp11, libgeos
│   Suggests: exactextractr, terra, controlledburn, vaster, sf
├── R/
│   ├── burn_sparse.R          # Main user-facing function
│   ├── dense_to_sparse.R      # Dense buffer → two-table conversion
│   ├── merge_runs_edges.R     # Utility: two-table → per-cell table
│   └── materialise_chunk.R    # Utility: sparse → dense matrix block
├── src/
│   ├── init.cpp               # libgeos_init_api() + cpp11 registration
│   ├── burn_sparse.cpp        # Core: geometry loop + GEOSGridIntersectionFractions
│   ├── dense_to_sparse.cpp    # Dense buffer → runs/edges (C++ for speed)
│   └── vendor/                # exactextract headers (temporary, until libgeos 3.14)
└── vignettes/
    └── introduction.Rmd
```

## Comparison with hypertidy packages

| Feature | controlledburn | fasterize | gridburn (this) |
|---------|---------------|-----------|-----------------|
| Output format | RLE sparse | Dense raster | Two-table sparse |
| Coverage fractions | No (binary) | No (binary) | Yes (exact) |
| Interior cells | RLE runs | Full raster | RLE runs |
| Edge cells | Included in runs | Full raster | Separate table with weights |
| Lines support | Partial | No | Deferred |
| Geometry input | WKB/sf | Any polygon vector | geos_geometry |
| C++ binding | Rcpp | Rcpp | cpp11 |
| GEOS dependency | No | No | Yes (via libgeos) |
| Memory for large grids | O(sparse) | O(dense) | O(sparse) |

## Next Concrete Step

Write the R prototype (Steps 1–4 above) and test it against a simple polygon
on a small grid. Verify the two-table output is correct and the RLE compression
is working. Then test against controlledburn for consistency on the runs table.

This can be done entirely in R with exactextractr + terra — no C++ needed yet.
The C++ package is a performance optimisation of the same algorithm.
