
#include "Rcpp.h"
using namespace Rcpp;
#include "edge.h"
#include "check_inputs.h"

#include "edgelist.h"
#include "rasterize.h"
#include "burn.h"
#include "utils.h"



// Rasterize an sf object of polygons without materializing grid values
//
// Rasterize set of polygons
//
// This is a  simplification of the high-performance replacement for [raster::rasterize()].
//
// The algorithm is based on the method described in course materials provided
// by [Wayne O. Cochran](https://labs.wsu.edu/wayne-cochran/). The algorithm
// is originally attributed to
// [Wylie et al. (1967)](https://dx.doi.org/10.1145/1465611.1465619).
//
// @param sf an [sf::sf()] object with a geometry column of POLYGON and/or
// MULTIPOLYGON objects.
// @param extent numeric vector c(xmin, xmax, ymin , ymax)
// @param dimension integer vector c(ncol, nrow)
// @return nothing atm
// @references Wylie, C., Romney, G., Evans, D., & Erdahl, A. (1967).
//   Half-tone perspective drawings by computer. Proceedings of the November
//   14-16, 1967, Fall Joint Computer Conference. AFIPS '67 (Fall).
//   <https://dx.doi.org/10.1145/1465611.1465619>
// [[Rcpp::export]]
Rcpp::List burn_polygon(Rcpp::DataFrame &sf,
                   Rcpp::NumericVector &extent,
                   Rcpp::IntegerVector &dimension) {

  Rcpp::List polygons;
  Rcpp::NumericVector field_vals;
  check_inputs_polygon(sf, polygons);  // Also fills in polygons

  //set up things we'll use later
  Rcpp::List::iterator p;
  Rcpp::NumericVector::iterator f;
  RasterInfo ras(extent, dimension);
  CollectorList out_vector;
    //Rasterize but always assign to the one layer
    p = polygons.begin();
    f = field_vals.begin();
     for(; p != polygons.end(); ++p, ++f) {
       rasterize_polygon( (*p), ras, out_vector);
    }

    return out_vector.vector();
}


// [[Rcpp::export]]
Rcpp::List burn_line(Rcpp::DataFrame &sf,
                        Rcpp::NumericVector &extent,
                        Rcpp::IntegerVector &dimension) {

  Rcpp::List lines;

  check_inputs_line(sf, lines);  // Also fills in polygons

  //set up things we'll use later
  Rcpp::List::iterator ln;
  Rcpp::NumericVector::iterator f;
  RasterInfo ras(extent, dimension);
  CollectorList out_vector;
  //Rasterize but always assign to the one layer
  ln = lines.begin();

  for(; ln != lines.end(); ++ln) {
    rasterize_line( (*ln), ras, out_vector);
  }

  return out_vector.vector();
}





