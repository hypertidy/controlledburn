## Validation: shared-boundary coverage complementarity
##
## When adjacent polygons share a boundary, coverage fractions at shared
## cells should sum to exactly 1.0 (no gaps, no overlaps).
##
## This tests the core property that makes polygon rasterization useful for
## area-weighting operations: every cell is fully accounted for.

library(controlledburn)
library(geos)

# Helper: materialise per-polygon coverage and check they sum to 1.0
check_complementarity <- function(result, label = "", tol = 1e-5) {
  ids <- unique(c(result$runs$id, result$edges$id))
  if (length(ids) < 2) {
    cat(sprintf("SKIP [%s]: only %d polygon(s)\n", label, length(ids)))
    return(TRUE)
  }

  nc <- result$dimension[1]
  nr <- result$dimension[2]

  # Sum all polygon coverages into one matrix
  total <- matrix(0, nrow = nr, ncol = nc)
  for (id in ids) {
    total <- total + materialise_chunk(result, id = id)
  }

  # Find cells that are touched by at least one polygon
  any_coverage <- total > tol

  # Check: touched cells should sum to 1.0
  bad <- any_coverage & abs(total - 1.0) > tol
  n_bad <- sum(bad)

  if (n_bad == 0) {
    n_touched <- sum(any_coverage)
    cat(sprintf("PASS [%s]: %d shared cells, all sum to 1.0 (± %.0e)\n",
                label, n_touched, tol))
    return(TRUE)
  } else {
    cat(sprintf("FAIL [%s]: %d cells have coverage sum != 1.0\n", label, n_bad))
    diffs <- which(bad, arr.ind = TRUE)
    head_diffs <- head(diffs, 10)
    for (i in seq_len(nrow(head_diffs))) {
      r <- head_diffs[i, 1]; c <- head_diffs[i, 2]
      cat(sprintf("  [%d,%d] sum=%.8f\n", r, c, total[r, c]))
    }
    return(FALSE)
  }
}

# Also check that scanline matches sparse
check_vs_sparse <- function(geoms, ext, dim, label = "") {
  r_sl <- burn_scanline(geoms, extent = ext, dimension = dim)
  r_sp <- burn_sparse(geoms, extent = ext, dimension = dim)

  mat_sl <- materialise_chunk(r_sl)
  mat_sp <- materialise_chunk(r_sp)
  max_diff <- max(abs(mat_sl - mat_sp))

  if (max_diff < 1e-5) {
    cat(sprintf("MATCH [%s]: scanline == sparse (max diff %.2e)\n", label, max_diff))
    return(TRUE)
  } else {
    n_diff <- sum(abs(mat_sl - mat_sp) > 1e-5)
    cat(sprintf("MISMATCH [%s]: %d cells differ, max diff %.6f\n", label, n_diff, max_diff))
    return(FALSE)
  }
}


cat("=== Shared boundary validation ===\n\n")

results <- list()

# ---- Test 1: Two adjacent rectangles sharing a vertical edge ----
cat("--- Test 1: Adjacent rectangles, vertical shared edge ---\n")
polys1 <- as_geos_geometry(c(
  "POLYGON ((0 0, 5 0, 5 10, 0 10, 0 0))",
  "POLYGON ((5 0, 10 0, 10 10, 5 10, 5 0))"
))
ext1 <- c(0, 10, 0, 10)
dim1 <- c(20L, 20L)

r1 <- burn_scanline(polys1, extent = ext1, dimension = dim1)
results[["adj_vert"]] <- check_complementarity(r1, "adj_vert 20x20")
results[["adj_vert_match"]] <- check_vs_sparse(polys1, ext1, dim1, "adj_vert 20x20")

# ---- Test 2: Same but boundary NOT on grid line ----
cat("\n--- Test 2: Adjacent rectangles, boundary between grid lines ---\n")
dim2 <- c(12L, 12L)  # 10/12 cells — boundary at x=5 falls mid-cell
r2 <- burn_scanline(polys1, extent = ext1, dimension = dim2)
results[["adj_vert_off"]] <- check_complementarity(r2, "adj_vert 12x12")
results[["adj_vert_off_match"]] <- check_vs_sparse(polys1, ext1, dim2, "adj_vert 12x12")

# ---- Test 3: Horizontal shared edge ----
cat("\n--- Test 3: Adjacent rectangles, horizontal shared edge ---\n")
polys3 <- as_geos_geometry(c(
  "POLYGON ((0 0, 10 0, 10 5, 0 5, 0 0))",
  "POLYGON ((0 5, 10 5, 10 10, 0 10, 0 5))"
))
dim3 <- c(15L, 15L)
r3 <- burn_scanline(polys3, extent = ext1, dimension = dim3)
results[["adj_horiz"]] <- check_complementarity(r3, "adj_horiz 15x15")
results[["adj_horiz_match"]] <- check_vs_sparse(polys3, ext1, dim3, "adj_horiz 15x15")

# ---- Test 4: Four quadrants (cross-shaped shared boundary) ----
cat("\n--- Test 4: Four quadrants ---\n")
polys4 <- as_geos_geometry(c(
  "POLYGON ((0 0, 5 0, 5 5, 0 5, 0 0))",
  "POLYGON ((5 0, 10 0, 10 5, 5 5, 5 0))",
  "POLYGON ((0 5, 5 5, 5 10, 0 10, 0 5))",
  "POLYGON ((5 5, 10 5, 10 10, 5 10, 5 5))"
))
dim4 <- c(16L, 16L)  # boundary at x=5, y=5 falls mid-cell
r4 <- burn_scanline(polys4, extent = ext1, dimension = dim4)
results[["quadrants"]] <- check_complementarity(r4, "quadrants 16x16")
results[["quadrants_match"]] <- check_vs_sparse(polys4, ext1, dim4, "quadrants 16x16")

# ---- Test 5: Diagonal shared edge (two triangles forming a square) ----
cat("\n--- Test 5: Diagonal shared edge ---\n")
polys5 <- as_geos_geometry(c(
  "POLYGON ((0 0, 10 0, 10 10, 0 0))",
  "POLYGON ((0 0, 10 10, 0 10, 0 0))"
))
dim5 <- c(20L, 20L)
r5 <- burn_scanline(polys5, extent = ext1, dimension = dim5)
results[["diagonal"]] <- check_complementarity(r5, "diagonal 20x20")
results[["diagonal_match"]] <- check_vs_sparse(polys5, ext1, dim5, "diagonal 20x20")

# ---- Test 6: Irregular tiling (L-shape + complement) ----
cat("\n--- Test 6: Irregular tiling (L-shape + complement) ---\n")
polys6 <- as_geos_geometry(c(
  "POLYGON ((0 0, 6 0, 6 4, 4 4, 4 10, 0 10, 0 0))",    # L-shape
  "POLYGON ((6 0, 10 0, 10 10, 4 10, 4 4, 6 4, 6 0))"    # complement
))
dim6 <- c(20L, 20L)
r6 <- burn_scanline(polys6, extent = ext1, dimension = dim6)
results[["irregular"]] <- check_complementarity(r6, "irregular 20x20")
results[["irregular_match"]] <- check_vs_sparse(polys6, ext1, dim6, "irregular 20x20")

# ---- Test 7: Donut — exterior polygon with interior polygon filling hole ----
cat("\n--- Test 7: Donut filled with inner polygon ---\n")
polys7 <- as_geos_geometry(c(
  "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))",  # outer with hole
  "POLYGON ((3 3, 7 3, 7 7, 3 7, 3 3))"                                 # inner fill
))
dim7 <- c(20L, 20L)
r7 <- burn_scanline(polys7, extent = ext1, dimension = dim7)
results[["donut_filled"]] <- check_complementarity(r7, "donut_filled 20x20")
results[["donut_filled_match"]] <- check_vs_sparse(polys7, ext1, dim7, "donut_filled 20x20")

# ---- Test 8: High-res adjacent triangles (many diagonal boundary cells) ----
cat("\n--- Test 8: High-res diagonal ---\n")
dim8 <- c(100L, 100L)
r8 <- burn_scanline(polys5, extent = ext1, dimension = dim8)
results[["diag_hires"]] <- check_complementarity(r8, "diagonal 100x100")
results[["diag_hires_match"]] <- check_vs_sparse(polys5, ext1, dim8, "diagonal 100x100")

# ---- Test 9: NC counties (real-world adjacent polygons) ----
if (requireNamespace("sf", quietly = TRUE)) {
  cat("\n--- Test 9: NC counties (adjacent real polygons) ---\n")
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc_geom <- sf::st_geometry(nc[1:10, ])
  ext9 <- as.numeric(sf::st_bbox(nc_geom))[c(1, 3, 2, 4)]
  dim9 <- c(500L, 200L)

  r9 <- burn_scanline(nc_geom, extent = ext9, dimension = dim9)

  # For real data, check total coverage doesn't exceed 1.0 per cell
  total9 <- matrix(0, nrow = dim9[2], ncol = dim9[1])
  ids9 <- unique(c(r9$runs$id, r9$edges$id))
  for (id in ids9) {
    total9 <- total9 + materialise_chunk(r9, id = id)
  }

  over <- sum(total9 > 1.0 + 1e-5)
  max_total <- max(total9)
  cat(sprintf("  %d cells exceed 1.0, max total = %.6f\n", over, max_total))
  results[["nc_no_overlap"]] <- (over == 0)
  if (over == 0) {
    cat(sprintf("PASS [nc10 500x200]: no cell exceeds 1.0\n"))
  } else {
    cat(sprintf("FAIL [nc10 500x200]: %d cells exceed 1.0\n", over))
  }

  results[["nc_match"]] <- check_vs_sparse(nc_geom, ext9, dim9, "nc10 500x200")
}

# ---- Summary ----
cat("\n=== Summary ===\n")
n_pass <- sum(unlist(results))
n_total <- length(results)
cat(sprintf("%d / %d tests passed\n", n_pass, n_total))
if (n_pass < n_total) {
  cat("Failed tests:", paste(names(results)[!unlist(results)], collapse = ", "), "\n")
}
