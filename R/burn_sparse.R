#' Sparse polygon rasterization with exact coverage fractions
#'
#' Computes exact coverage fractions for polygon-grid intersections and returns
#' results in a sparse two-table format: run-length encoded interior cells and
#' individually weighted boundary cells.
#'
#' @param x geometry input, one of:
#'   - an `sfc` geometry column (from sf)
#'   - a `geos_geometry` vector (from geos)
#'   - a list of raw vectors containing WKB
#' @param extent numeric vector `c(xmin, xmax, ymin, ymax)` defining the raster
#'   extent. If `NULL` (default), derived from the bounding box of `x`.
#' @param dimension integer vector `c(ncol, nrow)` defining the grid dimensions.
#'   If `NULL` (default), fitted to the extent with at most 256 cells along the
#'   longer axis, preserving aspect ratio. Mutually exclusive with `resolution`.
#' @param resolution numeric, cell size (scalar for square cells, or
#'   `c(dx, dy)`). If supplied, `dimension` is computed from `extent / resolution`.
#'   Mutually exclusive with `dimension`.
#' @param tile_size integer, maximum tile dimension (default 4096). The grid is
#'   processed in tiles of at most `tile_size x tile_size` cells to bound memory
#'   usage. Set to `Inf` to disable tiling.
#'
#' @return A list with class `"controlledburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`, `id` —
#'       run-length encoded interior cells (coverage fraction ≈ 1.0)}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `weight`, `id` —
#'       boundary cells with partial coverage (0 < weight < 1)}
#'     \item{`extent`}{the raster extent}
#'     \item{`dimension`}{the grid dimensions}
#'   }
#'
#'   Row and column indices are 1-based. Row 1 is the top (ymax) row.
#'   The `id` column is a 1-based index into the input geometry vector.
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'   poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
#'
#'   # Explicit extent and dimension
#'   result <- burn_sparse(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'
#'   # Defaults: extent from bbox, 256-cell fitted grid
#'   result <- burn_sparse(poly)
#'
#'   # Specify resolution (cell size)
#'   result <- burn_sparse(poly, resolution = 0.1)
#' }
burn_sparse <- function(x, extent = NULL, dimension = NULL, resolution = NULL,
                        tile_size = 4096L) {
  gp <- .resolve_grid_params(x, extent, dimension, resolution)
  extent <- gp$extent
  dimension <- gp$dimension

  wkb <- as_wkb_list(x)
  tile_size <- as.integer(tile_size)

  ncol_full <- dimension[1]
  nrow_full <- dimension[2]

  # If grid fits in a single tile, no tiling needed
  if (ncol_full <= tile_size && nrow_full <= tile_size) {
    result <- cpp_burn_sparse(
      wkb,
      extent[1], extent[3], extent[2], extent[4],
      ncol_full, nrow_full
    )
    result$extent <- extent
    result$dimension <- dimension
    class(result) <- "controlledburn"
    return(result)
  }

  # Tiled processing
  xmin <- extent[1]; xmax <- extent[2]
  ymin <- extent[3]; ymax <- extent[4]
  dx <- (xmax - xmin) / ncol_full
  dy <- (ymax - ymin) / nrow_full

  # Tile breaks (column and row indices, 0-based)
  col_breaks <- seq(0L, ncol_full, by = tile_size)
  if (col_breaks[length(col_breaks)] != ncol_full)
    col_breaks <- c(col_breaks, ncol_full)
  row_breaks <- seq(0L, nrow_full, by = tile_size)
  if (row_breaks[length(row_breaks)] != nrow_full)
    row_breaks <- c(row_breaks, nrow_full)

  n_tile_cols <- length(col_breaks) - 1L
  n_tile_rows <- length(row_breaks) - 1L

  all_runs <- vector("list", n_tile_cols * n_tile_rows)
  all_edges <- vector("list", n_tile_cols * n_tile_rows)
  k <- 0L

  for (ti in seq_len(n_tile_rows)) {
    r0 <- row_breaks[ti]       # 0-based row start
    r1 <- row_breaks[ti + 1L]  # 0-based row end (exclusive)
    tile_nrow <- r1 - r0

    # Tile y extent: row 0 is at ymax
    tile_ymax <- ymax - r0 * dy
    tile_ymin <- ymax - r1 * dy

    for (tj in seq_len(n_tile_cols)) {
      c0 <- col_breaks[tj]
      c1 <- col_breaks[tj + 1L]
      tile_ncol <- c1 - c0

      tile_xmin <- xmin + c0 * dx
      tile_xmax <- xmin + c1 * dx

      tile_result <- cpp_burn_sparse(
        wkb,
        tile_xmin, tile_ymin, tile_xmax, tile_ymax,
        tile_ncol, tile_nrow
      )

      k <- k + 1L

      # Offset tile-local coords to full grid coords
      runs <- tile_result$runs
      if (nrow(runs) > 0) {
        runs$row <- runs$row + r0
        runs$col_start <- runs$col_start + c0
        runs$col_end <- runs$col_end + c0
      }
      all_runs[[k]] <- runs

      edges <- tile_result$edges
      if (nrow(edges) > 0) {
        edges$row <- edges$row + r0
        edges$col <- edges$col + c0
      }
      all_edges[[k]] <- edges
    }
  }

  result <- list(
    runs = do.call(rbind, all_runs),
    edges = do.call(rbind, all_edges),
    extent = extent,
    dimension = dimension
  )
  # Clean up row names from rbind
  if (!is.null(result$runs)) rownames(result$runs) <- NULL
  if (!is.null(result$edges)) rownames(result$edges) <- NULL
  class(result) <- "controlledburn"
  result
}

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

#' @export
print.controlledburn <- function(x, ...) {
  ncol <- x$dimension[1]
  nrow <- x$dimension[2]
  n_runs <- nrow(x$runs)
  n_edges <- nrow(x$edges)
  n_ids <- length(unique(c(x$runs$id, x$edges$id)))

  # Compute total cells represented
  total_interior <- if (n_runs > 0) sum(as.numeric(x$runs$col_end - x$runs$col_start + 1L)) else 0
  total_cells <- total_interior + n_edges
  grid_cells <- as.numeric(ncol) * as.numeric(nrow)
  sparsity <- 1 - total_cells / grid_cells

  cat(sprintf("<controlledburn> %d x %d grid, %d geometr%s\n",
              ncol, nrow, n_ids, if (n_ids == 1) "y" else "ies"))
  cat(sprintf("  runs:  %d (%.0f interior cells)\n", n_runs, total_interior))
  cat(sprintf("  edges: %d boundary cells\n", n_edges))
  cat(sprintf("  sparsity: %.1f%% empty\n", sparsity * 100))
  invisible(x)
}

# ---- internal: geometry coercion to WKB list ----

as_wkb_list <- function(x) {
  UseMethod("as_wkb_list")
}

#' @export
as_wkb_list.list <- function(x) {
  # Assume already a list of raw vectors (WKB)
  stopifnot(all(vapply(x, is.raw, logical(1))))
  wk::wkb(x)
}

#' @export
as_wkb_list.blob <- function(x) {
  wk::wkb(x)
}

#' @export
as_wkb_list.wk_wkb <- function(x) {
  # Already wk_wkb — passthrough
  x
}

#' @export
as_wkb_list.sfc <- function(x) {
  # sf geometry column → WKB via wk
  wk::as_wkb(x)
}

#' @export
as_wkb_list.geos_geometry <- function(x) {
  # geos geometry → WKB via wk
  wk::as_wkb(x)
}
