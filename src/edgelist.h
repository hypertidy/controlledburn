#ifndef EDGELIST
#define EDGELIST

#include "edge.h"
#include "stdlib.h"
#include "Rcpp.h"
using namespace Rcpp;
extern void edgelist(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge> &edges);
extern void edgelist_line(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge> &edges);
extern NumericVector edgelist_out(std::list<Edge> edges);

#endif
