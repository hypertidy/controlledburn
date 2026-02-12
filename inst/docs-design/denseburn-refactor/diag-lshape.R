## Diagnostic: dump boundary cells and winding deltas for L-shape
##
## The L-shape is concave: POLYGON ((1 1, 5 1, 5 5, 3 5, 3 9, 1 9, 1 1))
## Grid: 20x20 over [0,10]x[0,10], cell size 0.5x0.5
##
## The vertical arm (x=1..3, y=5..9) has boundary edges at x=1 (downward)
## and x=3 (upward). Both are exactly on cell boundaries.
## The reflex vertex is at (5,5)->(3,5)->(3,9).

library(controlledburn)
library(geos)

poly <- as_geos_geometry("POLYGON ((1 1, 5 1, 5 5, 3 5, 3 9, 1 9, 1 1))")
ext <- c(0, 10, 0, 10)
dim <- c(20L, 20L)

# Reference
r_ref <- burn_sparse(poly, extent = ext, dimension = dim)
mat_ref <- materialise_chunk(r_ref)

# Scanline
r_scan <- burn_scanline(poly, extent = ext, dimension = dim)
mat_scan <- materialise_chunk(r_scan)

cat("=== L-shape grid geometry ===\n")
cat("Cell size: 0.5 x 0.5\n")
cat("Row 1 = y [9.5, 10.0] (top), Row 20 = y [0.0, 0.5] (bottom)\n")
cat("Col 1 = x [0.0, 0.5] (left), Col 20 = x [9.5, 10.0] (right)\n\n")

cat("=== Polygon boundary edges ===\n")
cat("(1,1)->(5,1): bottom, rightward, y=1 (row 18 boundary)\n")
cat("(5,1)->(5,5): right side lower arm, upward, x=5 (col 10 boundary)\n")
cat("(5,5)->(3,5): reflex, leftward, y=5 (row 10 boundary)\n")
cat("(3,5)->(3,9): right side upper arm, upward, x=3 (col 6 boundary)\n")
cat("(3,9)->(1,9): top, leftward, y=9 (row 2 boundary)\n")
cat("(1,9)->(1,1): left side, downward, x=1 (col 2 boundary)\n\n")

cat("=== Reference (burn_sparse) matrix ===\n")
# Show the interior region that should be filled
cat("Rows 3-18, cols 2-10 (covering the L-shape region):\n")
sub <- mat_ref[2:19, 2:11]
rownames(sub) <- paste0("r", 2:19)
colnames(sub) <- paste0("c", 2:11)
print(round(sub, 2))

cat("\n=== Difference (reference - scanline) ===\n")
diff <- mat_ref - mat_scan
cat("Non-zero differences:\n")
nz <- which(abs(diff) > 1e-6, arr.ind = TRUE)
if (nrow(nz) > 0) {
  for (i in seq_len(nrow(nz))) {
    r <- nz[i, 1]; c <- nz[i, 2]
    cat(sprintf("  [row=%d, col=%d] ref=%.4f scan=%.4f (y=[%.1f,%.1f] x=[%.1f,%.1f])\n",
                r, c, mat_ref[r,c], mat_scan[r,c],
                10 - r*0.5, 10 - (r-1)*0.5,
                (c-1)*0.5, c*0.5))
  }
} else {
  cat("  None!\n")
}

cat("\n=== Scanline sparse output ===\n")
cat("Runs:\n")
print(r_scan$runs)
cat("\nEdges:\n")
print(r_scan$edges)

cat("\n=== Reference sparse output ===\n")
cat("Runs:\n")
print(r_ref$runs)
