#' Sparse polygon rasterization with exact coverage fractions
#'
#' Computes exact coverage fractions for polygon-grid intersections via
#' an exactextract-style flood-fill, returning sparse two-table output:
#' run-length-encoded interior cells and individually weighted boundary
#' cells.
#'
#' Polygon-only. Line and point input are rejected with an error directing
#' the caller to [burn_scanline()], which produces the unified
#' (polygon / line / point) output schema.
#'
#' Compared to [burn_scanline()], this function is older and uses a
#' bounding-box-bounded dense intermediate (hence the `tile_size`
#' parameter). For polygon input the two functions produce numerically
#' identical results; [burn_scanline()] is preferred for new code.
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
#'   usage. Set to `Inf` to disable tiling. (Not relevant on [burn_scanline()],
#'   whose memory cost is O(perimeter).)
#'
#' @return A list with class `"controlledburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`, `id` —
#'       run-length encoded interior cells (full coverage)}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `fraction`, `id` —
#'       boundary cells with partial coverage, `fraction` in (0, 1)}
#'     \item{`extent`}{the raster extent}
#'     \item{`dimension`}{the grid dimensions}
#'   }
#'
#'   Row and column indices are 1-based. Row 1 is the top (ymax) row.
#'   The `id` column is a 1-based index into the input geometry vector.
#'
#' @details  controlledburn does not validate input geometry. It rasterizes what you pass it,
#' valid or not, and trusts you to judge the result.
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
#' Expands a sparse [burn_sparse()] / [burn_scanline()] result into a dense
#' matrix, optionally over a subwindow of the parent grid. The matrix's
#' values depend on which kind of geometry produced the result:
#'
#'   - **Polygon** (`$runs` + `$edges`): coverage fraction in `[0, 1]`.
#'     Interior cells are 1; boundary cells are the partial fraction.
#'   - **Line** (`$lines`): absolute length within each cell, in CRS units.
#'     Multiple lines crossing the same cell sum.
#'   - **Point** (`$points`): integer count of points landing in each
#'     cell. Multiple points in one cell sum.
#'
#' If a result holds tables for *more than one* geometry kind (which can
#' happen when a single burn was given a mixed-kind input vector),
#' `materialise_chunk()` errors rather than silently combining different
#' units in one matrix. Filter with `id =` to one kind first, or run
#' separate burns for each kind.
#'
#' @param x a `"controlledburn"` object from [burn_sparse()] or [burn_scanline()]
#' @param target numeric extent `c(xmin, xmax, ymin, ymax)` to materialise, or
#'   `NULL` (default) for the full grid. The target extent is snapped outward to
#'   cell boundaries of the source grid and clamped to the source extent.
#' @param id integer geometry id (or vector of ids) to materialise, or `NULL`
#'   (default) for all. Useful for filtering a mixed-kind result down to a
#'   single kind before materialising.
#' @param type character, one of `"matrix"` (default) or `"vector"`
#' @param max_cells maximum number of cells to allocate (default 1e8). Set to
#'   `Inf` to disable the safety check.
#'
#' @return A numeric matrix (nrow × ncol) or vector (length nrow*ncol).
#'   The cell unit depends on geometry kind (see Description). When `target`
#'   is specified, the matrix dimensions correspond to the snapped subwindow.
#'   The snapped extent is available as `attr(result, "extent")`.
#'
#' @export
materialise_chunk <- function(x, target = NULL, id = NULL,
                              type = c("matrix", "vector"),
                              max_cells = 1e8) {
  stopifnot(inherits(x, "controlledburn"))
  type <- match.arg(type)

  src_nc <- x$dimension[1]
  src_nr <- x$dimension[2]
  dx <- (x$extent[2] - x$extent[1]) / src_nc
  dy <- (x$extent[4] - x$extent[3]) / src_nr

  if (is.null(target)) {
    # Full grid
    col_off <- 0L
    row_off <- 0L
    out_nc <- src_nc
    out_nr <- src_nr
    out_extent <- x$extent
  } else {
    stopifnot(is.numeric(target), length(target) == 4L, target[2] > target[1], target[4] > target[3])

    # Snap outward to source cell boundaries, then clamp to source extent
    snap <- .snap_extent(target, x$extent, dx, dy)
    out_extent <- snap$extent

    col_off <- snap$col_off
    row_off <- snap$row_off
    out_nc  <- snap$ncol
    out_nr  <- snap$nrow

    if (out_nc <= 0L || out_nr <= 0L) {
      warning("Target extent does not intersect source grid")
      mat <- matrix(0, nrow = 0L, ncol = 0L)
      attr(mat, "extent") <- out_extent
      if (type == "vector") return(numeric(0))
      return(mat)
    }
  }

  n_cells <- as.double(out_nc) * as.double(out_nr)
  if (n_cells > max_cells) {
    stop(sprintf(
      "Requested window is %d x %d = %.0f cells (max_cells = %.0f). Use a smaller target extent or increase max_cells.",
      out_nc, out_nr, n_cells, max_cells), call. = FALSE)
  }

  mat <- matrix(0, nrow = out_nr, ncol = out_nc)

  runs   <- x$runs
  edges  <- x$edges
  # $lines and $points are added by burn_scanline; older burn_sparse results
  # don't have them. Guard with NULL checks throughout.
  lines  <- if (!is.null(x$lines))  x$lines  else NULL
  points <- if (!is.null(x$points)) x$points else NULL

  if (!is.null(id)) {
    runs   <- runs[runs$id %in% id, , drop = FALSE]
    edges  <- edges[edges$id %in% id, , drop = FALSE]
    if (!is.null(lines))  lines  <- lines[lines$id %in% id, , drop = FALSE]
    if (!is.null(points)) points <- points[points$id %in% id, , drop = FALSE]
  }

  # Refuse to materialise mixed-type tables. Each table type produces a
  # different kind of measure (polygon: fraction in [0,1]; line: length in
  # CRS units; point: count) and combining them in one matrix would silently
  # mix units. The walker can produce all four tables in one pass — that
  # capacity exists for downstream uses that want them separately — but
  # materialise_chunk is the consumer that would need to pick units, so
  # it errors instead.
  n_polygon <- nrow(runs) + nrow(edges)
  n_line    <- if (!is.null(lines))  nrow(lines)  else 0L
  n_point   <- if (!is.null(points)) nrow(points) else 0L
  n_kinds <- (n_polygon > 0) + (n_line > 0) + (n_point > 0)
  if (n_kinds > 1) {
    stop(
      "materialise_chunk(): result has multiple geometry types ",
      "(polygon = ", n_polygon, ", line = ", n_line, ", point = ", n_point, "). ",
      "Each type produces a different measure (fraction vs length vs count); ",
      "filter to one type with id = before calling, or run separate burns ",
      "for each type.",
      call. = FALSE
    )
  }

  # Fill interior runs (with offset + clipping for subwindow)
  if (nrow(runs) > 0) {
    for (i in seq_len(nrow(runs))) {
      r <- runs$row[i] - row_off
      if (r < 1L || r > out_nr) next
      cs <- runs$col_start[i] - col_off
      ce <- runs$col_end[i] - col_off
      cs <- max(cs, 1L)
      ce <- min(ce, out_nc)
      if (cs > ce) next
      mat[r, cs:ce] <- 1
    }
  }

  # Fill polygon edge cells with the coverage fraction (in [0, 1]).
  # The column is named `fraction` post Phase 4 to make the unit explicit;
  # tolerate the older `weight` name for forward-compat with cached results.
  if (nrow(edges) > 0) {
    fcol <- if ("fraction" %in% names(edges)) "fraction" else "weight"
    for (i in seq_len(nrow(edges))) {
      r <- edges$row[i] - row_off
      cc <- edges$col[i] - col_off
      if (r < 1L || r > out_nr || cc < 1L || cc > out_nc) next
      mat[r, cc] <- mat[r, cc] + edges[[fcol]][i]
    }
  }

  # Fill line cells with length-in-cell (absolute, CRS units).
  if (!is.null(lines) && nrow(lines) > 0) {
    for (i in seq_len(nrow(lines))) {
      r <- lines$row[i] - row_off
      cc <- lines$col[i] - col_off
      if (r < 1L || r > out_nr || cc < 1L || cc > out_nc) next
      mat[r, cc] <- mat[r, cc] + lines$length[i]
    }
  }

  # Fill point cells with implicit weight = 1.
  # Multiple points landing in the same cell accumulate (count of points
  # at this location), consistent with the line-overlap convention
  # ("lengths sum, not unioned"). Use distinct = TRUE upstream if you
  # want presence rather than count.
  if (!is.null(points) && nrow(points) > 0) {
    for (i in seq_len(nrow(points))) {
      r <- points$row[i] - row_off
      cc <- points$col[i] - col_off
      if (r < 1L || r > out_nr || cc < 1L || cc > out_nc) next
      mat[r, cc] <- mat[r, cc] + 1
    }
  }

  attr(mat, "extent") <- out_extent

  if (type == "vector") {
    v <- as.vector(t(mat))
    attr(v, "extent") <- out_extent
    v
  } else {
    mat
  }
}


# ---- Snap extent outward to source cell boundaries, clamp to source ----
#
# Returns a list with snapped extent, col/row offsets, and output dimensions.
# Snap-out: requested boundary moves outward to include any partially covered cell.
# Clamp: snapped boundary is clamped to the source grid extent (no implied expand).
# On-boundary tolerance: if a requested boundary is within tol * cell_size of a
# cell boundary, it's treated as exactly on that boundary (no snap needed).

.snap_extent <- function(target, src_extent, dx, dy) {
  tol <- 1e-8

  src_xmin <- src_extent[1]; src_xmax <- src_extent[2]
  src_ymin <- src_extent[3]; src_ymax <- src_extent[4]

  # Snap outward: floor for min edges, ceiling for max edges
  # (relative to source origin)
  x_off_min <- (target[1] - src_xmin) / dx
  x_off_max <- (target[2] - src_xmin) / dx
  y_off_min <- (target[3] - src_ymin) / dy
  y_off_max <- (target[4] - src_ymin) / dy

  # Snap: floor for min (outward = leftward/downward),
  #        ceil for max (outward = rightward/upward).
  # Tolerance: if within tol of integer, treat as on-boundary.
  snap_floor <- function(v) {
    r <- floor(v)
    if (abs(v - round(v)) < tol) r <- round(v)
    r
  }
  snap_ceil <- function(v) {
    r <- ceiling(v)
    if (abs(v - round(v)) < tol) r <- round(v)
    r
  }

  col_start_0 <- snap_floor(x_off_min)  # 0-based column start
  col_end_0   <- snap_ceil(x_off_max)    # 0-based column end (exclusive)
  row_end_0   <- snap_ceil(y_off_max)    # 0-based row end (from bottom)
  row_start_0 <- snap_floor(y_off_min)   # 0-based row start (from bottom)

  # Clamp to source grid
  src_nc <- round((src_xmax - src_xmin) / dx)
  src_nr <- round((src_ymax - src_ymin) / dy)

  col_start_0 <- max(col_start_0, 0L)
  col_end_0   <- min(col_end_0, src_nc)
  row_start_0 <- max(row_start_0, 0L)
  row_end_0   <- min(row_end_0, src_nr)

  # Output dimensions
  out_nc <- as.integer(col_end_0 - col_start_0)
  out_nr <- as.integer(row_end_0 - row_start_0)

  # Snapped extent in source coordinates
  snapped_xmin <- src_xmin + col_start_0 * dx
  snapped_xmax <- src_xmin + col_end_0 * dx
  snapped_ymin <- src_ymin + row_start_0 * dy
  snapped_ymax <- src_ymin + row_end_0 * dy

  # Row/col offsets: source row/col 1 maps to output row/col 1 - offset
  # Source rows are numbered top-to-bottom: row 1 = ymax.
  # row_start_0/row_end_0 are from bottom, so convert.
  row_off_top <- as.integer(src_nr - row_end_0)   # how many source rows above the window
  col_off <- as.integer(col_start_0)

  list(
    extent = c(snapped_xmin, snapped_xmax, snapped_ymin, snapped_ymax),
    col_off = col_off,
    row_off = row_off_top,
    ncol = out_nc,
    nrow = out_nr
  )
}

#' @export
print.controlledburn <- function(x, ...) {
  ncol <- x$dimension[1]
  nrow <- x$dimension[2]
  n_runs <- nrow(x$runs)
  n_edges <- nrow(x$edges)
  # $lines and $points are present when burn_scanline produced them; older
  # burn_sparse results don't have them.
  n_lines  <- if (!is.null(x$lines))  nrow(x$lines)  else 0L
  n_points <- if (!is.null(x$points)) nrow(x$points) else 0L
  n_ids <- length(unique(c(x$runs$id, x$edges$id, x$lines$id, x$points$id)))

  # Compute total cells represented
  total_interior <- if (n_runs > 0) sum(as.numeric(x$runs$col_end - x$runs$col_start + 1L)) else 0
  total_cells <- total_interior + n_edges + n_lines + n_points
  grid_cells <- as.numeric(ncol) * as.numeric(nrow)
  sparsity <- 1 - total_cells / grid_cells

  cat(sprintf("<controlledburn> %d x %d grid, %d geometr%s\n",
              ncol, nrow, n_ids, if (n_ids == 1) "y" else "ies"))
  cat(sprintf("  runs:   %d (%.0f interior cells)\n", n_runs, total_interior))
  cat(sprintf("  edges:  %d polygon boundary cells\n", n_edges))
  if (n_lines > 0) {
    cat(sprintf("  lines:  %d cells\n", n_lines))
  }
  if (n_points > 0) {
    cat(sprintf("  points: %d cells\n", n_points))
  }
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
