## Validation: burn_scanline() vs burn_sparse()
##
## Tests that the scanline winding-sweep algorithm produces identical output
## to the original flood-fill approach in gridburn.
##
## Run after: devtools::load_all() or R CMD INSTALL

library(controlledburn)
library(geos)

# ---- helper: compare two gridburn results ----

compare_gridburn <- function(a, b, label = "") {
  # Materialise both to dense matrices and compare
  mat_a <- materialise_chunk(a)
  mat_b <- materialise_chunk(b)

  if (!identical(dim(mat_a), dim(mat_b))) {
    cat(sprintf("FAIL [%s]: dimension mismatch %s vs %s\n",
                label, paste(dim(mat_a), collapse="x"), paste(dim(mat_b), collapse="x")))
    return(FALSE)
  }

  max_diff <- max(abs(mat_a - mat_b))
  n_diff <- sum(abs(mat_a - mat_b) > 1e-5)

  if (max_diff < 1e-5) {
    cat(sprintf("PASS [%s]: identical (max diff = %.2e)\n", label, max_diff))
    return(TRUE)
  } else {
    cat(sprintf("FAIL [%s]: %d cells differ, max diff = %.6f\n", label, n_diff, max_diff))
    # Show where differences occur
    diffs <- which(abs(mat_a - mat_b) > 1e-5, arr.ind = TRUE)
    head_diffs <- head(diffs, 10)
    for (i in seq_len(nrow(head_diffs))) {
      r <- head_diffs[i, 1]
      c <- head_diffs[i, 2]
      cat(sprintf("  [%d,%d] sparse=%.6f scanline=%.6f\n", r, c, mat_a[r,c], mat_b[r,c]))
    }
    return(FALSE)
  }
}

# Also compare the sparse representations directly
compare_sparse <- function(a, b, label = "") {
  # Sort runs for comparison
  sort_runs <- function(r) {
    if (nrow(r) == 0) return(r)
    r[order(r$id, r$row, r$col_start), ]
  }
  sort_edges <- function(e) {
    if (nrow(e) == 0) return(e)
    e[order(e$id, e$row, e$col), ]
  }

  ra <- sort_runs(a$runs)
  rb <- sort_runs(b$runs)

  cat(sprintf("  [%s] runs: %d vs %d, edges: %d vs %d\n",
              label, nrow(a$runs), nrow(b$runs), nrow(a$edges), nrow(b$edges)))
}


cat("=== Scanline burn validation ===\n\n")

results <- list()

# ---- Test 1: Simple rectangle ----
cat("--- Test 1: Simple rectangle ---\n")
poly1 <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
ext1 <- c(0, 3, 0, 3)
dim1 <- c(3L, 3L)

r1_sparse <- burn_sparse(poly1, extent = ext1, dimension = dim1)
r1_scan   <- burn_scanline(poly1, extent = ext1, dimension = dim1)
compare_sparse(r1_sparse, r1_scan, "rect 3x3")
results[["rect_3x3"]] <- compare_gridburn(r1_sparse, r1_scan, "rect 3x3")

# ---- Test 2: Rectangle at higher resolution ----
cat("\n--- Test 2: Rectangle at higher resolution ---\n")
dim2 <- c(30L, 30L)
r2_sparse <- burn_sparse(poly1, extent = ext1, dimension = dim2)
r2_scan   <- burn_scanline(poly1, extent = ext1, dimension = dim2)
compare_sparse(r2_sparse, r2_scan, "rect 30x30")
results[["rect_30x30"]] <- compare_gridburn(r2_sparse, r2_scan, "rect 30x30")

# ---- Test 3: Triangle ----
cat("\n--- Test 3: Triangle ---\n")
poly3 <- as_geos_geometry("POLYGON ((1 1, 9 1, 5 9, 1 1))")
ext3 <- c(0, 10, 0, 10)
dim3 <- c(20L, 20L)

r3_sparse <- burn_sparse(poly3, extent = ext3, dimension = dim3)
r3_scan   <- burn_scanline(poly3, extent = ext3, dimension = dim3)
compare_sparse(r3_sparse, r3_scan, "triangle 20x20")
results[["triangle"]] <- compare_gridburn(r3_sparse, r3_scan, "triangle 20x20")

# ---- Test 4: Polygon with hole (donut) ----
cat("\n--- Test 4: Polygon with hole ---\n")
poly4 <- as_geos_geometry(
  "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"
)
ext4 <- c(0, 10, 0, 10)
dim4 <- c(20L, 20L)

r4_sparse <- burn_sparse(poly4, extent = ext4, dimension = dim4)
r4_scan   <- burn_scanline(poly4, extent = ext4, dimension = dim4)
compare_sparse(r4_sparse, r4_scan, "donut 20x20")
results[["donut"]] <- compare_gridburn(r4_sparse, r4_scan, "donut 20x20")

# ---- Test 5: L-shaped polygon ----
cat("\n--- Test 5: L-shaped polygon ---\n")
poly5 <- as_geos_geometry(
  "POLYGON ((1 1, 5 1, 5 5, 3 5, 3 9, 1 9, 1 1))"
)
ext5 <- c(0, 10, 0, 10)
dim5 <- c(20L, 20L)

r5_sparse <- burn_sparse(poly5, extent = ext5, dimension = dim5)
r5_scan   <- burn_scanline(poly5, extent = ext5, dimension = dim5)
compare_sparse(r5_sparse, r5_scan, "L-shape 20x20")
results[["L_shape"]] <- compare_gridburn(r5_sparse, r5_scan, "L-shape 20x20")

# ---- Test 6: Multiple polygons ----
cat("\n--- Test 6: Multiple polygons ---\n")
polys6 <- as_geos_geometry(c(
  "POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
  "POLYGON ((6 6, 9 6, 9 9, 6 9, 6 6))"
))
ext6 <- c(0, 10, 0, 10)
dim6 <- c(20L, 20L)

r6_sparse <- burn_sparse(polys6, extent = ext6, dimension = dim6)
r6_scan   <- burn_scanline(polys6, extent = ext6, dimension = dim6)
compare_sparse(r6_sparse, r6_scan, "multi 20x20")
results[["multi_poly"]] <- compare_gridburn(r6_sparse, r6_scan, "multi 20x20")

# ---- Test 7: Diamond (diagonal edges) ----
cat("\n--- Test 7: Diamond (diagonal edges) ---\n")
poly7 <- as_geos_geometry("POLYGON ((5 1, 9 5, 5 9, 1 5, 5 1))")
ext7 <- c(0, 10, 0, 10)
dim7 <- c(40L, 40L)

r7_sparse <- burn_sparse(poly7, extent = ext7, dimension = dim7)
r7_scan   <- burn_scanline(poly7, extent = ext7, dimension = dim7)
compare_sparse(r7_sparse, r7_scan, "diamond 40x40")
results[["diamond"]] <- compare_gridburn(r7_sparse, r7_scan, "diamond 40x40")

# ---- Test 8: nc.shp counties (real data, if sf available) ----
if (requireNamespace("sf", quietly = TRUE)) {
  cat("\n--- Test 8: NC counties (first 5) ---\n")
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc_geom <- sf::st_geometry(nc[1:5, ])
  ext8 <- as.numeric(sf::st_bbox(nc_geom))[c(1, 3, 2, 4)]  # xmin, xmax, ymin, ymax
  dim8 <- c(500L, 200L)

  r8_sparse <- burn_sparse(nc_geom, extent = ext8, dimension = dim8)
  r8_scan   <- burn_scanline(nc_geom, extent = ext8, dimension = dim8)
  compare_sparse(r8_sparse, r8_scan, "nc5 500x200")
  results[["nc5"]] <- compare_gridburn(r8_sparse, r8_scan, "nc5 500x200")
} else {
  cat("\n--- Test 8: SKIPPED (sf not available) ---\n")
}

# ---- Summary ----
cat("\n=== Summary ===\n")
n_pass <- sum(unlist(results))
n_total <- length(results)
cat(sprintf("%d / %d tests passed\n", n_pass, n_total))
if (n_pass < n_total) {
  cat("Failed tests:", paste(names(results)[!unlist(results)], collapse = ", "), "\n")
}

# ---- Default parameter tests ----
cat("\n=== Default parameter tests ===\n")

# bare minimum: just geometry
cat("--- bare call (defaults) ---\n")
r_default <- burn_scanline(poly1)
cat(sprintf("  extent from bbox: [%.1f, %.1f, %.1f, %.1f]\n",
            r_default$extent[1], r_default$extent[2],
            r_default$extent[3], r_default$extent[4]))
cat(sprintf("  dimension fitted: %d x %d\n",
            r_default$dimension[1], r_default$dimension[2]))

# resolution parameter
cat("--- resolution = 0.1 ---\n")
r_res <- burn_scanline(poly1, resolution = 0.1)
cat(sprintf("  dimension: %d x %d (expect 20 x 20)\n",
            r_res$dimension[1], r_res$dimension[2]))

# resolution + explicit extent
cat("--- resolution = 0.5 with explicit extent ---\n")
r_res2 <- burn_scanline(poly1, extent = c(0, 10, 0, 5), resolution = 0.5)
cat(sprintf("  dimension: %d x %d (expect 20 x 10)\n",
            r_res2$dimension[1], r_res2$dimension[2]))

# error on both dimension and resolution
cat("--- dimension + resolution (should error) ---\n")
tryCatch({
  burn_scanline(poly1, dimension = c(10, 10), resolution = 0.1)
  cat("  FAIL: no error raised\n")
}, error = function(e) {
  cat(sprintf("  OK: %s\n", e$message))
})
