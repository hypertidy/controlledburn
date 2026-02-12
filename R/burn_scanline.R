#' Scanline polygon rasterization with exact coverage fractions
#'
#' Experimental scanline sweep implementation of sparse polygon rasterization.
#' Produces identical output to [burn_sparse()] but uses a winding-number
#' scanline sweep instead of flood fill for interior classification.
#' Memory usage is O(perimeter) rather than O(bounding-box area).
#'
#' @inheritParams burn_sparse
#'
#' @return A list with class `"controlledburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`, `id`}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `weight`, `id`}
#'     \item{`extent`}{the raster extent}
#'     \item{`dimension`}{the grid dimensions}
#'   }
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'   poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
#'
#'   # Defaults: extent from bbox, 256-cell fitted grid
#'   result <- burn_scanline(poly)
#'
#'   # Explicit
#'   result <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'
#'   # By resolution
#'   result <- burn_scanline(poly, resolution = 0.01)
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
