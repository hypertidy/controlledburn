#ifndef LASERIZE
#define LASERIZE


extern Rcpp::S4 laserize(Rcpp::DataFrame &sf,
                          Rcpp::S4 &raster,
                          std::string field,
                          std::string fun,
                          double background,
                          std::string by);
#endif
