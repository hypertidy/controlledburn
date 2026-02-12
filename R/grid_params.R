# grid_params.R — shared grid parameter resolution for burn_sparse / burn_scanline

# uses wk::wk_bbox (imported in controlledburn-package.R)
.resolve_grid_params <- function(x, extent, dimension, resolution) {

  # ---- extent: derive from geometry if NULL ----
  if (is.null(extent)) {
    bb <- as.numeric(wk::wk_bbox(x))  # xmin, ymin, xmax, ymax
    extent <- bb[c(1L, 3L, 2L, 4L)]   # → xmin, xmax, ymin, ymax
  }
  extent <- as.double(extent)
  stopifnot(
    length(extent) == 4,
    extent[2] > extent[1],  # xmax > xmin
    extent[4] > extent[3]   # ymax > ymin
  )

  # ---- dimension vs resolution: mutually exclusive ----
  if (!is.null(dimension) && !is.null(resolution)) {
    stop("specify 'dimension' or 'resolution', not both", call. = FALSE)
  }

  if (!is.null(resolution)) {
    resolution <- as.double(resolution)
    if (length(resolution) == 1L) resolution <- rep(resolution, 2L)
    stopifnot(length(resolution) == 2, all(resolution > 0))
    # ncol from x-range, nrow from y-range
    dimension <- as.integer(ceiling(
      c(extent[2] - extent[1], extent[4] - extent[3]) / resolution
    ))
  }

  if (is.null(dimension)) {
    dimension <- .fit_dims(256L, wh = c(extent[2] - extent[1], extent[4] - extent[3]))
  }

  dimension <- as.integer(dimension)
  stopifnot(length(dimension) == 2, dimension[1] > 0, dimension[2] > 0)

  list(extent = extent, dimension = dimension)
}

# Fit grid dimensions to a maximum cell count along the longer axis,
# preserving aspect ratio.
.fit_dims <- function(size = 256L, wh = c(size, size)) {
  w <- wh[1L]; h <- wh[2L]
  as.integer(ceiling(size * c(w, h) / max(w, h)))
}
