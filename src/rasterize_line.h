#ifndef RASTERIZE_LINE
#define RASTERIZE_LINE

#include "edge.h"
#include "Rcpp.h"
#include "CollectorList.h"

using namespace Rcpp;

extern void rasterize_line(Rcpp::RObject polygon,
                              RasterInfo &ras, CollectorList &out_vector);
#endif
