# benchmark-scaling.R — Item 4: Benchmark perimeter-proportional scaling
#
# Compare burn_scanline (O(perimeter) memory, analytical coverage) against
# burn_sparse (O(bbox area) memory, dense matrix + flood fill).
#
# Key hypothesis: for a fixed polygon, doubling grid resolution
# quadruples burn_sparse cost but only doubles burn_scanline cost
# (boundary cells scale with perimeter, interior runs are O(1) per row).

library(controlledburn)
library(geos)

# ---- Helper: robust timing ----
time_one <- function(expr_fn, times = 3) {
  # expr_fn is a zero-arg function to call
  timings <- numeric(times)
  for (i in seq_len(times)) {
    gc(FALSE)
    timings[i] <- system.time(expr_fn())[["elapsed"]]
  }
  median(timings)
}

# ---- Test geometries ----

# 1. Simple square — compact, low perimeter:area
square_wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))"

# 2. Thin horizontal sliver — high perimeter:area
sliver_wkt <- "POLYGON ((0.5 4.5, 9.5 4.5, 9.5 5.5, 0.5 5.5, 0.5 4.5))"

# 3. Star shape — complex boundary, many diagonal edges
star_coords <- function(n = 12, r_outer = 4.5, r_inner = 2.0, cx = 5, cy = 5) {
  angles <- seq(0, 2 * pi, length.out = 2 * n + 1)
  radii <- rep(c(r_outer, r_inner), length.out = 2 * n + 1)
  x <- cx + radii * cos(angles)
  y <- cy + radii * sin(angles)
  paste0("POLYGON ((", paste(sprintf("%.6f %.6f", x, y), collapse = ", "), "))")
}
star_wkt <- star_coords(12)

# 4. Donut — hole means more boundary per area
donut_wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"

# 5. Fractal-ish coastline — many vertices, high perimeter
jagged_coords <- function(n_teeth = 40, amplitude = 0.3) {
  base_x <- c(seq(1, 9, length.out = n_teeth),
               rep(9, n_teeth),
               seq(9, 1, length.out = n_teeth),
               rep(1, n_teeth))
  base_y <- c(rep(1, n_teeth),
               seq(1, 9, length.out = n_teeth),
               rep(9, n_teeth),
               seq(9, 1, length.out = n_teeth))
  n <- length(base_x)
  perturb <- rep(c(0, amplitude), length.out = n)
  px <- numeric(n); py <- numeric(n)
  side_n <- n_teeth
  px[1:side_n] <- base_x[1:side_n]
  py[1:side_n] <- base_y[1:side_n] + perturb[1:side_n]
  idx <- (side_n+1):(2*side_n)
  px[idx] <- base_x[idx] - perturb[idx]
  py[idx] <- base_y[idx]
  idx <- (2*side_n+1):(3*side_n)
  px[idx] <- base_x[idx]
  py[idx] <- base_y[idx] - perturb[idx]
  idx <- (3*side_n+1):(4*side_n)
  px[idx] <- base_x[idx] + perturb[idx]
  py[idx] <- base_y[idx]
  px <- c(px, px[1])
  py <- c(py, py[1])
  paste0("POLYGON ((", paste(sprintf("%.6f %.6f", px, py), collapse = ", "), "))")
}
jagged_wkt <- jagged_coords(40)

geoms <- list(
  square = as_geos_geometry(square_wkt),
  sliver = as_geos_geometry(sliver_wkt),
  star   = as_geos_geometry(star_wkt),
  donut  = as_geos_geometry(donut_wkt),
  jagged = as_geos_geometry(jagged_wkt)
)

# ---- Grid resolution sweep ----
#
# Scale grid from coarse to fine while holding extent fixed at (0, 10, 0, 10).
# Cell count scales as n^2 but perimeter cell count scales as n.
#
# We push to high resolutions where the dense approach starts to hurt.
# burn_sparse uses tiling (4096 default) so it won't OOM, but it does
# more work per tile.

ext <- c(0, 10, 0, 10)
resolutions <- c(100, 200, 400, 800, 1600, 3200)

cat("=== Scaling benchmark: burn_scanline vs burn_sparse ===\n")
cat(sprintf("Resolutions: %s\n", paste(resolutions, collapse = ", ")))
cat(sprintf("Geometries: %s\n", paste(names(geoms), collapse = ", ")))
cat(sprintf("Extent: [%g, %g, %g, %g]\n\n", ext[1], ext[2], ext[3], ext[4]))

results <- list()

for (gname in names(geoms)) {
  g <- geoms[[gname]]
  cat(sprintf("--- %s ---\n", gname))
  cat(sprintf("%8s  %12s  %12s  %8s  %8s  %8s\n",
              "res", "scanline_ms", "sparse_ms", "speedup",
              "edges", "runs"))

  for (res in resolutions) {
    dim <- c(res, res)

    t_scanline <- time_one(function() burn_scanline(g, extent = ext, dimension = dim)) * 1000
    t_sparse   <- time_one(function() burn_sparse(g, extent = ext, dimension = dim)) * 1000

    r_sl <- burn_scanline(g, extent = ext, dimension = dim)
    speedup <- t_sparse / t_scanline

    cat(sprintf("%8d  %12.2f  %12.2f  %8.1fx  %8d  %8d\n",
                res, t_scanline, t_sparse, speedup,
                nrow(r_sl$edges), nrow(r_sl$runs)))

    results[[length(results) + 1]] <- data.frame(
      geometry = gname, resolution = res,
      scanline_ms = t_scanline, sparse_ms = t_sparse,
      speedup = speedup,
      edges = nrow(r_sl$edges), runs = nrow(r_sl$runs),
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

results_df <- do.call(rbind, results)

# ---- Scaling analysis ----
cat("=== Scaling rates (log2 of time ratio between consecutive resolutions) ===\n")
cat("O(n) ≈ 1.0    O(n log n) ≈ 1.x    O(n²) ≈ 2.0\n\n")

for (gname in unique(results_df$geometry)) {
  sub <- results_df[results_df$geometry == gname, ]
  if (nrow(sub) < 2) next

  sl_rates <- diff(log2(sub$scanline_ms))
  sp_rates <- diff(log2(sub$sparse_ms))
  res_labels <- paste0(sub$resolution[-1], "/", sub$resolution[-nrow(sub)])

  cat(sprintf("--- %s ---\n", gname))
  cat(sprintf("%12s  %10s  %10s\n", "step", "scanline", "sparse"))
  for (i in seq_along(sl_rates)) {
    cat(sprintf("%12s  %10.2f  %10.2f\n", res_labels[i], sl_rates[i], sp_rates[i]))
  }
  cat(sprintf("%12s  %10.2f  %10.2f\n\n", "median",
              median(sl_rates), median(sp_rates)))
}

# ---- Edge count scaling ----
cat("=== Edge count scaling (log2 ratio between consecutive resolutions) ===\n")
cat("O(n) ≈ 1.0 — boundary cells grow linearly with resolution\n\n")

for (gname in unique(results_df$geometry)) {
  sub <- results_df[results_df$geometry == gname, ]
  if (nrow(sub) < 2) next

  edge_rates <- diff(log2(pmax(sub$edges, 1)))
  res_labels <- paste0(sub$resolution[-1], "/", sub$resolution[-nrow(sub)])

  cat(sprintf("%-8s  ", gname))
  for (i in seq_along(edge_rates)) {
    cat(sprintf("%s:%.2f  ", res_labels[i], edge_rates[i]))
  }
  cat(sprintf(" median: %.2f\n", median(edge_rates)))
}

# ---- Memory scaling ----
cat("\n=== Memory comparison at highest resolution ===\n")
max_res <- max(resolutions)
cat(sprintf("Grid: %d x %d = %.1fM cells\n", max_res, max_res,
            as.numeric(max_res)^2 / 1e6))
cat(sprintf("Dense matrix: %.1f MB (float)\n",
            as.numeric(max_res)^2 * 4 / 1024^2))

for (gname in names(geoms)) {
  g <- geoms[[gname]]
  r <- burn_scanline(g, extent = ext, dimension = c(max_res, max_res))
  sparse_bytes <- object.size(r$runs) + object.size(r$edges)
  cat(sprintf("%-8s sparse: %.1f KB (%d runs, %d edges)\n",
              gname, as.numeric(sparse_bytes) / 1024,
              nrow(r$runs), nrow(r$edges)))
}

cat("\n=== Summary ===\n")
cat("If scanline scaling ≈ 1.0 and sparse ≈ 2.0, the O(perimeter) hypothesis\n")
cat("is confirmed: doubling resolution doubles scanline cost but quadruples\n")
cat("dense cost. The speedup factor grows with resolution.\n")
