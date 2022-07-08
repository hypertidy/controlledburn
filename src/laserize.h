#ifndef LASERIZE
#define LASERIZE


extern Rcpp::IntegerVector laserize(Rcpp::DataFrame &sf,
                         Rcpp::NumericVector &extent,
                         Rcpp::IntegerVector &dimension);
#endif
