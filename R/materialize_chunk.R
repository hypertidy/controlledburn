## prototype for working with align
materialise_chunk_align <- function(burn, target, what = c("coverage", "id")) {
  what <- match.arg(what)
  parent <- as_vast(burn)
  if (missing(target)) target <- align::as_vast(burn)


  if (!is_aligned(target, parent)) {
    stop("target is not aligned to burn grid")
  }

  off <- offset(target, parent)
  nr <- nrow(target)
  nc <- ncol(target)

  row_lo <- off[["row"]] + 1L
  row_hi <- off[["row"]] + nr
  col_lo <- off[["col"]] + 1L
  col_hi <- off[["col"]] + nc

  m <- matrix(NA_real_, nrow = nr, ncol = nc)

  ## runs: interior spans, full coverage
  runs <- burn$runs
  if (!is.null(runs) && nrow(runs) > 0L) {
    keep <- runs$row >= row_lo & runs$row <= row_hi &
      runs$col_end >= col_lo & runs$col_start <= col_hi
    if (any(keep)) {
      r <- runs[keep, , drop = FALSE]
      local_row   <- r$row - off[["row"]]
      local_start <- pmax(r$col_start - off[["col"]], 1L)
      local_end   <- pmin(r$col_end   - off[["col"]], nc)
      val <- if (what == "id") r$id else 1.0
      for (i in seq_len(nrow(r))) {
        cols <- local_start[i]:local_end[i]
        m[local_row[i], cols] <- val[if (length(val) == 1) 1L else i]
      }
    }
  }

  ## edges: boundary cells with coverage fraction
  edges <- burn$edges
  if (!is.null(edges) && nrow(edges) > 0L) {
    keep <- edges$row >= row_lo & edges$row <= row_hi &
      edges$col >= col_lo & edges$col <= col_hi
    if (any(keep)) {
      e <- edges[keep, , drop = FALSE]
      local_row <- e$row - off[["row"]]
      local_col <- e$col - off[["col"]]
      val <- if (what == "id") e$id else e$weight
      for (i in seq_len(nrow(e))) {
        m[local_row[i], local_col[i]] <- val[i]
      }
    }
  }

  m
}





## initial version
#' Materialise a controlledburn result to a dense matrix or vector
#'
#' Expands the sparse two-table representation into a dense coverage fraction
#' matrix. Primarily for visualisation and testing.
#'
#' @param x a `"controlledburn"` object from [burn_sparse()]
#' @param id integer polygon id to materialise, or `NULL` (default) for all
#' @param type character, one of `"matrix"` (default) or `"vector"`
#'
#' @return A numeric matrix (nrow × ncol) or vector (length nrow*ncol) of
#'   coverage fractions. Values range from 0 (outside) to 1 (fully inside).
#'
#' @export

materialise_chunk <- function(x, id = NULL, type = c("matrix", "vector")) {
  stopifnot(inherits(x, "controlledburn"))
  type <- match.arg(type)

  ncol <- x$dimension[1]
  nrow <- x$dimension[2]

  mat <- matrix(0, nrow = nrow, ncol = ncol)

  runs <- x$runs
  edges <- x$edges

  if (!is.null(id)) {
    runs <- runs[runs$id %in% id, , drop = FALSE]
    edges <- edges[edges$id %in% id, , drop = FALSE]
  }

  # Fill interior runs
  if (nrow(runs) > 0) {
    for (i in seq_len(nrow(runs))) {
      r <- runs$row[i]
      cs <- runs$col_start[i]
      ce <- runs$col_end[i]
      mat[r, cs:ce] <- 1
    }
  }

  # Fill edge cells
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      r <- edges$row[i]
      c <- edges$col[i]
      w <- edges$weight[i]
      # For overlapping polygons, sum weights (or take max depending on use case)
      mat[r, c] <- mat[r, c] + w
    }
  }

  if (type == "vector") {
    as.vector(t(mat))  # row-major order
  } else {
    mat
  }
}
