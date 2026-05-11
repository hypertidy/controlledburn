#' Sparse rasterization of polygon, line, and point geometry
#'
#' Rasterizes geometry into a sparse output schema, with one type-pure
#' table per geometry kind. Polygons produce run-length-encoded interior
#' cells (`$runs`) and weighted boundary cells (`$edges`, with a
#' coverage `fraction` in `[0, 1]`). Lines produce one record per touched
#' cell with the absolute `length` of line within the cell, in CRS units
#' (`$lines`). Points produce one record per touched cell, with no
#' measure column (`$points`). See the unified geometry rasterization
#' design doc in `inst/docs-design/`.
#'
#' Mixed geometry-kind input (e.g. a vector containing both polygons
#' and lines) is permitted and produces tables populated for whichever
#' kinds are present. Each table's measure column means exactly one
#' thing — polygon `$edges$fraction` is dimensionless, `$lines$length`
#' is in CRS units, points have no measure — so combining them in a
#' single materialized matrix would silently mix units.
#' [materialize_chunk()] therefore refuses mixed-kind input; filter to
#' one kind via `id =` first, or run separate burns for each kind.
#'
#' GeometryCollection input warns and is skipped. Curved types
#' (CircularString, CompoundCurve, etc.) must be linearised by the
#' caller before calling.
#'
#' Memory usage is O(perimeter) rather than O(bounding-box area), so
#' no `tile_size` parameter is required (in contrast to [burn_sparse()],
#' which is polygon-only and tile-bounded).
#'
#' @inheritParams burn_sparse
#'
#' @return A list with class `"controlledburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`,
#'       `id` — run-length encoded polygon interior cells (always full
#'       coverage). Empty for line/point input.}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `fraction`,
#'       `id` — polygon boundary cells with partial coverage; `fraction`
#'       is in (0, 1). Empty for line/point input.}
#'     \item{`lines`}{data.frame with columns `row`, `col`, `length`,
#'       `id` — line cells; `length` is the absolute length of the line
#'       within the cell, in CRS units. Empty for polygon/point input.}
#'     \item{`points`}{data.frame with columns `row`, `col`, `id` —
#'       point cells; no measure column (a point is either in a cell or
#'       it isn't). Empty for polygon/line input.}
#'     \item{`extent`}{the raster extent (`c(xmin, xmax, ymin, ymax)`)}
#'     \item{`dimension`}{the grid dimensions (`c(ncol, nrow)`)}
#'   }
#'
#'   Row and column indices are 1-based. Row 1 is the top (ymax) row.
#'   The `id` column is a 1-based index into the input geometry vector.
#'
#'   The measure is computed in the CRS the geometry and grid are
#'   specified in — *not* the geodesic measure on the Earth. Geometry
#'   in degrees produces area in square degrees and length in degrees.
#'   For geodesic measures, project to an appropriate CRS first.
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'
#'   # Polygon: produces $runs (interior) and $edges (fraction in [0, 1])
#'   poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
#'   r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'
#'   # Line: produces $lines with absolute length in CRS units
#'   line <- as_geos_geometry("LINESTRING (0.5 0.5, 2.5 2.5)")
#'   r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'
#'   # Point: produces $points with row/col only (no measure column)
#'   pt <- as_geos_geometry("POINT (1.5 1.5)")
#'   r <- burn_scanline(pt, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'
#'   # Defaults: extent from bbox, 256-cell fitted grid
#'   r <- burn_scanline(poly)
#'
#'   # By resolution
#'   r <- burn_scanline(poly, resolution = 0.01)
#' }
burn_scanline <- function(x, extent = NULL, dimension = NULL, resolution = NULL) {
  gp <- .resolve_grid_params(x, extent, dimension, resolution)
  extent <- gp$extent
  dimension <- gp$dimension

  wkb <- as_wkb_list(x)

  ncol_full <- dimension[1]
  nrow_full <- dimension[2]

  result <- cpp_scanline_burn(
    wkb,
    extent[1], extent[3], extent[2], extent[4],
    ncol_full, nrow_full
  )

  result$extent <- extent
  result$dimension <- dimension
  class(result) <- "controlledburn"
  result
}
