#ifndef RASTERIZE_POLYGON
#define RASTERIZE_POLYGON

#include "edge.h"
#include "Rcpp.h"
using namespace Rcpp;

extern void rasterize_polygon(Rcpp::RObject polygon,
                              RasterInfo &ras, IntegerVector &out_vector);
#endif
