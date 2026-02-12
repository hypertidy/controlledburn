## gridburn prototype
## Sparse polygon rasterization with exact coverage fractions
## Uses exactextractr as the computational backend for prototyping
##
## Output: two-table sparse representation
##   runs:  (row, col_start, col_end, id) — interior cells, weight = 1
##   edges: (row, col, weight, id) — boundary cells, 0 < weight < 1

#' Convert a dense per-cell coverage table to the two-table sparse format
#'
#' @param cell integer vector of 1-based cell indices (column-major, as terra uses)
#' @param coverage_fraction numeric vector of coverage fractions (0–1)
#' @param ncol integer, number of columns in the grid
#' @param id integer, geometry index
#' @param tol numeric, tolerance for "fully covered" (default 1 - 1e-6)
#' @return named list with `runs` and `edges` data.frames
dense_to_sparse <- function(cell, coverage_fraction, ncol, id, tol = 1 - 1e-6) {
  ## terra cell indexing is 1-based, row-major:
  ##   cell = (row - 1) * ncol + col
  ## so: row = ((cell - 1) %/% ncol) + 1
  ##     col = ((cell - 1) %%  ncol) + 1
  row <- ((cell - 1L) %/% ncol) + 1L
  col <- ((cell - 1L) %%  ncol) + 1L

  is_interior <- coverage_fraction >= tol

  ## --- Edges: individual cells with fractional weight ---
  edge_idx <- which(!is_interior & coverage_fraction > 0)
  edges <- data.frame(
    row    = row[edge_idx],
    col    = col[edge_idx],
    weight = coverage_fraction[edge_idx],
    id     = rep.int(id, length(edge_idx))
  )


  ## --- Runs: RLE-compress interior cells ---
  int_idx <- which(is_interior)
  if (length(int_idx) == 0L) {
    runs <- data.frame(
      row       = integer(0),
      col_start = integer(0),
      col_end   = integer(0),
      id        = integer(0)
    )
    return(list(runs = runs, edges = edges))
  }

  int_row <- row[int_idx]
  int_col <- col[int_idx]

  ## sort by row then col (should already be, but be safe)
  ord <- order(int_row, int_col)
  int_row <- int_row[ord]
  int_col <- int_col[ord]

  ## detect run boundaries: new row OR non-consecutive column
  n <- length(int_row)
  if (n == 1L) {
    runs <- data.frame(
      row       = int_row,
      col_start = int_col,
      col_end   = int_col,
      id        = id
    )
    return(list(runs = runs, edges = edges))
  }

  row_break <- c(TRUE, diff(int_row) != 0L)
  col_break <- c(TRUE, diff(int_col) != 1L)
  is_start  <- row_break | col_break

  run_start_idx <- which(is_start)
  run_end_idx   <- c(run_start_idx[-1L] - 1L, n)

  runs <- data.frame(
    row       = int_row[run_start_idx],
    col_start = int_col[run_start_idx],
    col_end   = int_col[run_end_idx],
    id        = rep.int(id, length(run_start_idx))
  )

  list(runs = runs, edges = edges)
}


#' Sparse polygon rasterization with exact coverage fractions
#'
#' @param geometry sf/sfc POLYGON or MULTIPOLYGON geometry
#' @param extent numeric(4): c(xmin, xmax, ymin, ymax)
#' @param dimension integer(2): c(ncol, nrow)
#' @param exact logical, if TRUE compute coverage fractions (default TRUE)
#' @return named list with `runs` and `edges` data.frames
#'
#' @details
#' Uses exactextractr as the computational backend.
#' The `runs` table contains interior cells (fully covered, weight = 1)
#' in RLE format: row, col_start, col_end, id.
#' The `edges` table contains boundary cells with fractional coverage:
#' row, col, weight, id.
#'
#' In approximate mode (exact = FALSE), only the `runs` table is returned
#' (edges is empty).
burn_sparse <- function(geometry, extent, dimension, exact = TRUE) {
  stopifnot(requireNamespace("terra", quietly = TRUE))
  stopifnot(requireNamespace("exactextractr", quietly = TRUE))

  ncol_grid <- dimension[1]
  nrow_grid <- dimension[2]

  ## Create a minimal raster matching the grid spec
  r <- terra::rast(
    nrows = nrow_grid, ncols = ncol_grid,
    xmin = extent[1], xmax = extent[2],
    ymin = extent[3], ymax = extent[4],
    crs = ""
  )
  terra::values(r) <- 1L

  ## ensure geometry is sfc
  if (inherits(geometry, "sf")) {
    geometry <- sf::st_geometry(geometry)
  }

  n_geom <- length(geometry)

  ## Get per-cell coverage fractions
  results <- exactextractr::exact_extract(
    r, geometry,
    fun = NULL,
    include_cell = TRUE,
    progress = FALSE
  )

  all_runs  <- vector("list", n_geom)
  all_edges <- vector("list", n_geom)

  for (i in seq_len(n_geom)) {
    df <- results[[i]]
    if (is.null(df) || nrow(df) == 0L) {
      all_runs[[i]]  <- data.frame(row = integer(0), col_start = integer(0),
                                    col_end = integer(0), id = integer(0))
      all_edges[[i]] <- data.frame(row = integer(0), col = integer(0),
                                    weight = numeric(0), id = integer(0))
      next
    }

    sparse <- dense_to_sparse(
      cell = df$cell,
      coverage_fraction = df$coverage_fraction,
      ncol = ncol_grid,
      id = i
    )

    if (!exact) {
      ## In approximate mode, include edge cells as single-cell runs
      ## (weight threshold > 0 means "touched")
      if (nrow(sparse$edges) > 0L) {
        extra_runs <- data.frame(
          row       = sparse$edges$row,
          col_start = sparse$edges$col,
          col_end   = sparse$edges$col,
          id        = sparse$edges$id
        )
        sparse$runs <- rbind(sparse$runs, extra_runs)
        ## re-sort
        ord <- order(sparse$runs$id, sparse$runs$row, sparse$runs$col_start)
        sparse$runs <- sparse$runs[ord, , drop = FALSE]
      }
      sparse$edges <- data.frame(row = integer(0), col = integer(0),
                                  weight = numeric(0), id = integer(0))
    }

    all_runs[[i]]  <- sparse$runs
    all_edges[[i]] <- sparse$edges
  }

  list(
    runs  = do.call(rbind, all_runs),
    edges = do.call(rbind, all_edges)
  )
}


#' Merge runs and edges into a per-cell table
#'
#' @param sparse named list with `runs` and `edges` as from burn_sparse()
#' @return data.frame with columns: row, col, weight, id
merge_runs_edges <- function(sparse) {
  runs  <- sparse$runs
  edges <- sparse$edges

  ## Expand runs to per-cell
  if (nrow(runs) > 0L) {
    expanded <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
      cols <- runs$col_start[i]:runs$col_end[i]
      data.frame(
        row    = rep.int(runs$row[i], length(cols)),
        col    = cols,
        weight = rep.int(1.0, length(cols)),
        id     = rep.int(runs$id[i], length(cols))
      )
    }))
  } else {
    expanded <- data.frame(row = integer(0), col = integer(0),
                           weight = numeric(0), id = integer(0))
  }

  ## Combine with edges
  if (nrow(edges) > 0L) {
    result <- rbind(expanded, edges)
  } else {
    result <- expanded
  }

  ## Sort
  ord <- order(result$id, result$row, result$col)
  result[ord, , drop = FALSE]
}


#' Materialise a chunk of the sparse representation as a dense matrix
#'
#' @param sparse named list with `runs` and `edges` as from burn_sparse()
#' @param row_range integer(2): c(row_start, row_end) inclusive
#' @param col_range integer(2): c(col_start, col_end) inclusive
#' @param id integer, which geometry to materialise (or NULL for all)
#' @return numeric matrix of coverage fractions
materialise_chunk <- function(sparse, row_range, col_range, id = NULL) {
  nrow_chunk <- row_range[2] - row_range[1] + 1L
  ncol_chunk <- col_range[2] - col_range[1] + 1L

  mat <- matrix(0.0, nrow = nrow_chunk, ncol = ncol_chunk)

  runs  <- sparse$runs
  edges <- sparse$edges

  if (!is.null(id)) {
    runs  <- runs[runs$id == id, , drop = FALSE]
    edges <- edges[edges$id == id, , drop = FALSE]
  }

  ## Filter to chunk extent
  runs <- runs[runs$row >= row_range[1] & runs$row <= row_range[2] &
               runs$col_end >= col_range[1] & runs$col_start <= col_range[2], ,
               drop = FALSE]

  edges <- edges[edges$row >= row_range[1] & edges$row <= row_range[2] &
                 edges$col >= col_range[1] & edges$col <= col_range[2], ,
                 drop = FALSE]

  ## Fill runs (clipped to chunk)
  if (nrow(runs) > 0L) {
    for (i in seq_len(nrow(runs))) {
      r <- runs$row[i] - row_range[1] + 1L
      c_start <- max(runs$col_start[i], col_range[1]) - col_range[1] + 1L
      c_end   <- min(runs$col_end[i],   col_range[2]) - col_range[1] + 1L
      mat[r, c_start:c_end] <- 1.0
    }
  }

  ## Fill edges
  if (nrow(edges) > 0L) {
    for (i in seq_len(nrow(edges))) {
      r <- edges$row[i] - row_range[1] + 1L
      c <- edges$col[i] - col_range[1] + 1L
      ## max() handles overlapping geometries — take highest coverage
      mat[r, c] <- max(mat[r, c], edges$weight[i])
    }
  }

  mat
}
