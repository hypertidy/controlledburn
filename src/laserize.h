#ifndef CONTROLLEDBURN
#define CONTROLLEDBURN


extern List laserize(Rcpp::DataFrame &sf,
                         Rcpp::NumericVector &extent,
                         Rcpp::IntegerVector &dimension, Rcpp::LogicalVector lines);
#endif
