#ifndef LASERIZE
#define LASERIZE


extern Rcpp::List laserize(Rcpp::DataFrame &sf,
                         Rcpp::NumericVector &extent,
                         Rcpp::IntegerVector &dimension);
#endif
