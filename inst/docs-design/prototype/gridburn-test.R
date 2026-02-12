## gridburn prototype — test and comparison with hypertidy packages
## Run this locally where exactextractr, terra, sf, controlledburn are available

source("inst/prototype/gridburn-prototype.R")

library(sf)
library(terra)

## ---- Test 1: Simple triangle on a 10x10 grid ----

# Grid: 10 cols × 10 rows, extent [0, 10, 0, 10]
extent <- c(0, 10, 0, 10)  # xmin, xmax, ymin, ymax
dimension <- c(10L, 10L)   # ncol, nrow

# Triangle polygon
coords <- matrix(c(2, 2,  8, 2,  5, 8,  2, 2), ncol = 2, byrow = TRUE)
tri <- st_sfc(st_polygon(list(coords)))

# Exact mode
result_exact <- burn_sparse(tri, extent, dimension, exact = TRUE)
cat("=== Exact mode ===\n")
cat("Runs (interior, weight=1):\n")
print(result_exact$runs)
cat("\nEdges (boundary, weight<1):\n")
print(result_exact$edges)

# Approximate mode
result_approx <- burn_sparse(tri, extent, dimension, exact = FALSE)
cat("\n=== Approximate mode ===\n")
cat("Runs (all touched cells):\n")
print(result_approx$runs)

# Merge to per-cell table
merged <- merge_runs_edges(result_exact)
cat("\n=== Merged per-cell table ===\n")
print(merged)

# Materialise as dense matrix
mat <- materialise_chunk(result_exact, c(1L, 10L), c(1L, 10L))
cat("\n=== Dense matrix (coverage fractions) ===\n")
print(round(mat, 2))


## ---- Test 2: Compare approximate mode with controlledburn ----

if (requireNamespace("controlledburn", quietly = TRUE)) {
  library(controlledburn)

  # controlledburn wants: geometry as WKB or sf, raster spec as extent + dimension
  # Its output is: list(row, col_start, col_end, id) — same as our runs table

  # Use a realistic-ish polygon set
  nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc <- st_transform(nc, 32119)  # NC state plane

  bb <- st_bbox(nc)
  ext <- c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
  dim <- c(500L, 500L)

  # gridburn approximate
  t1 <- system.time(gb <- burn_sparse(nc, ext, dim, exact = FALSE))
  cat("\ngridburn approximate:", nrow(gb$runs), "runs in", t1["elapsed"], "s\n")

  # controlledburn
  r <- rast(nrows = dim[2], ncols = dim[1],
            xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4])

  t2 <- system.time({
    cb <- controlledburn:::burn_polygon(nc, as.vector(ext(r)), dim(r)[2:1])
    # controlledburn returns list of vectors; convert to data.frame
    cb_df <- data.frame(
      row = cb$row,
      col_start = cb$col_start,
      col_end = cb$col_end,
      id = cb$id
    )
    cb_df <- do.call(rbind, cb) |>
      as.data.frame() |>
      setNames(c("row", "col_start", "col_end", "id")) |>
      tibble::as_tibble()
  })
  cat("controlledburn:", nrow(cb_df), "runs in", t2["elapsed"], "s\n")

  # Compare: the runs should be very similar
  # (not identical because controlledburn uses scanline/centre-of-cell,
  #  gridburn approximate uses exactextract coverage > 0 which is "all touched")
  cat("\nNote: outputs won't be identical because:\n")
  cat("  - controlledburn: cell-centre intersection (like GDAL default)\n")
  cat("  - gridburn approx: any coverage > 0 (like GDAL ALL_TOUCHED)\n")
  cat("  The gridburn approx will have MORE cells at polygon edges.\n")
}


## ---- Test 3: Compare exact mode with exactextractr coverage_fraction ----

if (requireNamespace("exactextractr", quietly = TRUE)) {
  # Direct exactextractr call for comparison
  r <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  values(r) <- 1

  # exactextractr returns per-feature data.frames
  ee_result <- exactextractr::exact_extract(r, tri, fun = NULL,
                                             include_cell = TRUE,
                                             include_xy = TRUE)[[1]]

  cat("\n=== exactextractr raw output (first feature) ===\n")
  print(ee_result)

  # Our merged table should have the same cells and weights
  # (modulo floating point tolerance for the 1.0 threshold)
  cat("\nCell count: exactextractr =", nrow(ee_result),
      ", gridburn merged =", nrow(merged), "\n")
}


## ---- Test 4: Check vaster compatibility ----

if (requireNamespace("vaster", quietly = TRUE)) {
  library(vaster)

  # vaster provides cell <-> coordinate conversion for a grid spec
  # Our runs/edges tables use (row, col) which vaster can convert to
  # geographic coordinates or cell indices

  # Example: convert edge cells to coordinates
  if (nrow(result_exact$edges) > 0) {
    edge_cells <- vaster::cell_from_row_col(
      dimension, result_exact$edges$row, result_exact$edges$col
    )

    edge_xy <- vaster::xy_from_cell(dimension, extent, edge_cells)
    cat("\n=== Edge cell coordinates (via vaster) ===\n")
    print(head(edge_xy))
  }
}
