#ifndef EDGELIST
#define EDGELIST

#include "edge.h"
#include "stdlib.h"
#include "Rcpp.h"
using namespace Rcpp;
extern void edgelist(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge> &edges);
#endif
