#' @title Deprecated: use burn() instead
#'
#' @description
#' `burn_sparse()` has been removed. Use [burn()] instead, which
#' produces identical output for polygon input and also supports
#' lines and points.
#'
#' @param ... ignored
#' @export
burn_sparse <- function(...) {
  .Deprecated("burn", msg = paste(
    "burn_sparse() has been removed.",
    "Use burn() instead (identical output for polygons,",
    "plus line and point support)."
  ))
  stop("burn_sparse() is no longer available. Use burn().")
}

#' @export
print.controlledburn <- function(x, ...) {
  ncol <- x$dimension[1]
  nrow <- x$dimension[2]
  n_runs <- nrow(x$runs)
  n_edges <- nrow(x$edges)
  n_lines  <- if (!is.null(x$lines))  nrow(x$lines)  else 0L
  n_points <- if (!is.null(x$points)) nrow(x$points) else 0L
  n_ids <- length(unique(c(x$runs$id, x$edges$id, x$lines$id, x$points$id)))

  total_interior <- if (n_runs > 0)
    sum(as.numeric(x$runs$col_end - x$runs$col_start + 1L)) else 0
  total_cells <- total_interior + n_edges + n_lines + n_points
  grid_cells <- as.numeric(ncol) * as.numeric(nrow)
  sparsity <- 1 - total_cells / grid_cells

  cat(sprintf("<controlledburn> %d x %d grid, %d geometr%s\n",
              ncol, nrow, n_ids, if (n_ids == 1) "y" else "ies"))
  cat(sprintf("  runs:   %d (%.0f interior cells)\n", n_runs, total_interior))
  cat(sprintf("  edges:  %d polygon boundary cells\n", n_edges))
  if (n_lines > 0) cat(sprintf("  lines:  %d cells\n", n_lines))
  if (n_points > 0) cat(sprintf("  points: %d cells\n", n_points))
  cat(sprintf("  sparsity: %.1f%% empty\n", sparsity * 100))
  invisible(x)
}

#' Materialise a controlledburn sparse result into a dense matrix
#'
#' Converts the sparse run/edge/line tables from [burn()] into a dense
#' matrix. Polygon runs and edges are summed; lines contribute their
#' per-cell length; points contribute 1.
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param fun Currently only `"sum"` (default).
#' @return A numeric matrix with dimensions `nrow × ncol`.
#' @export
materialise_chunk <- function(x, fun = "sum") {
  stopifnot(inherits(x, "controlledburn"))
  ext <- x$extent
  dim <- x$dimension
  ncol_full <- dim[1]
  nrow_full <- dim[2]

  mat <- matrix(0.0, nrow = nrow_full, ncol = ncol_full)

  if (nrow(x$runs) > 0) {
    for (i in seq_len(nrow(x$runs))) {
      r <- x$runs$row[i]
      cs <- x$runs$col_start[i]
      ce <- x$runs$col_end[i]
      mat[r, cs:ce] <- 1
    }
  }

  if (nrow(x$edges) > 0) {
    for (i in seq_len(nrow(x$edges))) {
      r <- x$edges$row[i]
      c <- x$edges$col[i]
      mat[r, c] <- mat[r, c] + x$edges$fraction[i]
    }
  }

  if (!is.null(x$lines) && nrow(x$lines) > 0) {
    for (i in seq_len(nrow(x$lines))) {
      r <- x$lines$row[i]
      c <- x$lines$col[i]
      mat[r, c] <- mat[r, c] + x$lines$length[i]
    }
  }

  if (!is.null(x$points) && nrow(x$points) > 0) {
    for (i in seq_len(nrow(x$points))) {
      r <- x$points$row[i]
      c <- x$points$col[i]
      mat[r, c] <- mat[r, c] + 1
    }
  }

  mat
}

#' @rdname materialise_chunk
#' @export
materialize_chunk <- materialise_chunk
