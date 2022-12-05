#include "Rcpp.h"
using namespace Rcpp;

void check_inputs_polygon(Rcpp::DataFrame &sf,
                  Rcpp::List &polygons) {

  std::stringstream err_msg;

  if(!Rf_inherits(sf, "sf")) {
    err_msg << "sf must be of class sf." << std::endl;
  }

  polygons = sf[Rcpp::as<std::string>(sf.attr("sf_column"))];

  if(!(Rf_inherits(polygons, "sfc_MULTIPOLYGON") |
     Rf_inherits(polygons, "sfc_POLYGON"))) {
    err_msg << "sf geometry must be POLYGON or MULTIPOLYGON" << std::endl;
  }



  std::string err = err_msg.str();
  if(!err.empty()) {
    Rcpp::stop(err);
  }
}
