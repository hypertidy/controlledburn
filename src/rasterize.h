#ifndef RASTERIZE_POLYGON
#define RASTERIZE_POLYGON

#include "edge.h"
#include "Rcpp.h"
#include "CollectorList.h"

using namespace Rcpp;

extern void rasterize_polygon(Rcpp::RObject polygon,
                              RasterInfo &ras, CollectorList &out_vector, unsigned int poly_id);
extern void rasterize_line(Rcpp::RObject polygon,
                              RasterInfo &ras, CollectorList &out_vector);
#endif
