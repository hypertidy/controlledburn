
#include "Rcpp.h"
using namespace Rcpp;
#include "edge.h"
//#include "check_inputs.h"

#include "edgelist.h"
#include "rasterize_polygon.h"
#include "laserize.h"
#include "utils.h"


// Rasterize an sf object of polygons
//
// Rasterize set of polygons
//
// This is a high-performance replacement for [raster::rasterize()].
//
// The algorithm is based on the method described in course materials provided
// by [Wayne O. Cochran](https://labs.wsu.edu/wayne-cochran/). The algorithm
// is originally attributed to
// [Wylie et al. (1967)](https://dx.doi.org/10.1145/1465611.1465619).
//
// @param sf an [sf::sf()] object with a geometry column of POLYGON and/or
// MULTIPOLYGON objects.
// @param raster A raster  object. Used as a template for the raster output.
// Can be created with [raster::raster()].
// The fasterize package provides a method to create a raster object from
// an sf object.
// @param field character. The name of a column in `sf`,
// providing a value for each of the polygons rasterized. If NULL (default),
// all polygons will be given a value of 1.
// @param fun character. The name of a function by which to combine overlapping
// polygons. Currently takes "sum", "first", "last", "min", "max", "count", or
// "any".  Future versions may include more functions or the ability to pass
// custom R/C++ functions. If you need to summarize by a different function,
// use `by=` to get a RasterBrick and then [raster::stackApply()] or
// [raster::calc()] to summarize.
// @param background numeric. Value to put in the cells that are not covered by
// any of the features of x. Default is NA.
// @param by character.  The name of a column in `sf` by which to aggregate
// layers.  If set, fasterize will return a RasterBrick with as many layers
// as unique values of the `by` column.
// @return A raster of the same size, extent, resolution and projection as the
// provided raster template.
// @references Wylie, C., Romney, G., Evans, D., & Erdahl, A. (1967).
//   Half-tone perspective drawings by computer. Proceedings of the November
//   14-16, 1967, Fall Joint Computer Conference. AFIPS '67 (Fall).
//   <https://dx.doi.org/10.1145/1465611.1465619>
// @examples
// library(sf)
// library(laserize)
// p1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60), c(-180,-20))
// hole <- rbind(c(-150,-20), c(-100,-10), c(-110,20), c(-150,-20))
// p1 <- list(p1, hole)
// p2 <- list(rbind(c(-10,0), c(140,60), c(160,0), c(140,-55), c(-10,0)))
// p3 <- list(rbind(c(-125,0), c(0,60), c(40,5), c(15,-45), c(-125,0)))
// pols <- st_sf(value = rep(1,3),
//               geometry = st_sfc(lapply(list(p1, p2, p3), st_polygon)))
// r <- raster(pols, res = 1)
// r <- laserize(pols, r, field = "value", fun="sum")
// plot(r)
// @export
// [[Rcpp::export]]
Rcpp::S4 laserize(Rcpp::DataFrame &sf,
                   Rcpp::S4 &raster,
                   Rcpp::Nullable<std::string> field = R_NilValue,
                   std::string fun = "last",
                   double background = NA_REAL,
                   Rcpp::Nullable<std::string> by = R_NilValue) {

  Rcpp::List polygons;
  Rcpp::NumericVector field_vals;
  // check_inputs(sf, raster, field, fun, background, by,
  //              polygons, field_vals);  // Also fills in polygons and field_vals

  //set up things we'll use later
  Rcpp::List::iterator p;
  Rcpp::NumericVector::iterator f;
  RasterInfo ras(raster);

    //Rasterize but always assign to the one layer
    p = polygons.begin();
    f = field_vals.begin();
     for(; p != polygons.end(); ++p, ++f) {
       rasterize_polygon( (*p), ras);
    }

    return raster;
}



