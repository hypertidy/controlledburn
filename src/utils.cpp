#include "Rcpp.h"
using namespace Rcpp;

Rcpp::CharacterVector as_character(const Rcpp::RObject vec) {
  if (vec.inherits("factor")) {
    Rcpp::IntegerVector ints(vec);
    Rcpp::StringVector levels = ints.attr("levels");
    Rcpp::CharacterVector out(ints.size());
    for(int i = 0; i < ints.size(); i++) {
      out[i] =
        (ints[i] == NA_INTEGER) ? NA_STRING : levels[ints[i] - 1];
    }
    return out;
  } else {
    return Rf_coerceVector(vec, STRSXP);
  }
}

