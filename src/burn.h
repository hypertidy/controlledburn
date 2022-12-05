#ifndef CONTROLLEDBURN
#define CONTROLLEDBURN


extern List burn_polygon(Rcpp::DataFrame &sf,
                         Rcpp::NumericVector &extent,
                         Rcpp::IntegerVector &dimension);
#endif
