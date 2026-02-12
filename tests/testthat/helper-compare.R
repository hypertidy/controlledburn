# Helper: compare burn_scanline vs burn_sparse on materialised matrices
expect_scanline_matches_sparse <- function(geoms, ext, dim, tol = 1e-5,
                                           label = NULL) {
  r_sl <- burn_scanline(geoms, extent = ext, dimension = dim)
  r_sp <- burn_sparse(geoms, extent = ext, dimension = dim)
  mat_sl <- materialise_chunk(r_sl)
  mat_sp <- materialise_chunk(r_sp)
  max_diff <- max(abs(mat_sl - mat_sp))
  expect_lt(max_diff, tol, label = label)
}

# Helper: check coverage complementarity across polygon ids
expect_complementary <- function(result, tol = 1e-5, label = NULL) {
  ids <- unique(c(result$runs$id, result$edges$id))
  nc <- result$dimension[1]
  nr <- result$dimension[2]

  total <- matrix(0, nrow = nr, ncol = nc)
  for (id in ids) {
    total <- total + materialise_chunk(result, id = id)
  }

  any_coverage <- total > tol
  bad <- any_coverage & abs(total - 1.0) > tol
  expect_equal(sum(bad), 0L, label = label)
}
