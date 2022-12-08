#ifndef CHECK_INPUTS
#define CHECK_INPUTS


extern void check_inputs_polygon(Rcpp::DataFrame &sf,
                         Rcpp::List &polygons);

extern void check_inputs_line(Rcpp::DataFrame &sf,
                                 Rcpp::List &lines);


#endif
