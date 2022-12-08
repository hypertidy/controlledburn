#include "edge.h"
#include "edgelist.h"

//  Builds an edge list from a polygon or multipolygon
void edgelist_polygon(Rcpp::RObject polygon, RasterInfo &ras, std::list<Edge_polygon> &edges) {
  //std::list<Edge> edges;
  double y0, y1, y0c, y1c;
  //iterate recursively over the list
  switch(polygon.sexp_type()) {
  case REALSXP: {
    //if the object is numeric, it an Nx2 matrix of polygon nodes.
    Rcpp::NumericMatrix poly(polygon);
    //Add edge to list if it's not horizontal and is in the raster
    for(int i = 0; i < (poly.nrow() - 1); ++i) {
      y0 = (ras.ymax - poly(i, 1))/ras.yres - 0.5;
      y1 = (ras.ymax - poly(i+1, 1))/ras.yres - 0.5;
      if(y0 > 0 || y1 > 0) {  //only both with edges that are in the raster
        y0c = std::ceil(y0);
        y1c = std::ceil(y1);
        if(y0c != y1c) {  //only bother with non-horizontal edges
          edges.push_back(Edge_polygon(poly(i    , 0), y0,
                               poly(i + 1, 0), y1, ras, y0c, y1c));
        }
      }
    }

    break;
  };
  case VECSXP: {
    //if the object is a list, recurse over that
    Rcpp::List polylist = Rcpp::as<Rcpp::List>(polygon);
    for(Rcpp::List::iterator it = polylist.begin();
        it != polylist.end();
        ++it) {
      edgelist_polygon(Rcpp::wrap(*it), ras, edges);
    }

    break;
  }
  default: {
    Rcpp::stop("incompatible SEXP; only accepts lists and REALSXPs");
  }
  }
}



void edgelist_line(Rcpp::RObject line, RasterInfo &ras, std::list<Edge_line> &edges) {

  //iterate recursively over the list
  switch(line.sexp_type()) {
  case REALSXP: {
    //if the object is numeric, it an Nx2 matrix of line nodes.
    Rcpp::NumericMatrix lns(line);
    //Add edge to list
    for(int i = 0; i < (lns.nrow() - 1); ++i) {
// cull any that aren't in the raster at all TODO
        edges.push_back(Edge_line(lns(i, 0), lns(i, 1),
                                         lns(i + 1, 0), lns(i + 1, 1), ras));
    }

    break;
  };
  case VECSXP: {
    //if the object is a list, recurse over that
    Rcpp::List linelist = Rcpp::as<Rcpp::List>(line);
    for(Rcpp::List::iterator it = linelist.begin();
        it != linelist.end();
        ++it) {
      edgelist_line(Rcpp::wrap(*it), ras, edges);
    }

    break;
  }
  default: {
    Rcpp::stop("incompatible SEXP; only accepts lists and REALSXPs");
  }
  }
}
