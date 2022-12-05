#ifndef EDGELIST
#define EDGELIST

#include "edge.h"
#include "stdlib.h"
#include "Rcpp.h"
using namespace Rcpp;
extern void edgelist_polygon(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge_polygon> &edges);
//extern void edgelist_line(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge_line> &edges);

#endif
