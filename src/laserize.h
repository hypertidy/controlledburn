#ifndef CONTROLLEDBURN
#define CONTROLLEDBURN


extern List laserize(Rcpp::DataFrame &sf,
                         Rcpp::NumericVector &extent,
                         Rcpp::IntegerVector &dimension);

extern List laserize_line(Rcpp::DataFrame &sf,
                     Rcpp::NumericVector &extent,
                     Rcpp::IntegerVector &dimension);

extern NumericVector edge_out(Rcpp::DataFrame &sf,
                             Rcpp::NumericVector &extent,
                             Rcpp::IntegerVector &dimension);
#endif
