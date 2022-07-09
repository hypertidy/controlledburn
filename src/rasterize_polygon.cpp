#include "edge.h"
#include "edgelist.h"
#include "rasterize_polygon.h"

// Rasterize a single polygon
// Based on https://ezekiel.encs.vancouver.wsu.edu/~cs442/lectures/rasterization/polyfill/polyfill.pdf #nolint

void record_cell(Rcpp::IntegerVector &out_vector, unsigned int xs, unsigned int xe, unsigned int y) {
  out_vector.push_back(xs);
  out_vector.push_back(xe);
  out_vector.push_back(y);
  return;
}
void rasterize_polygon(Rcpp::RObject polygon,
                       RasterInfo &ras, IntegerVector &out_vector) {

  std::list<Edge>::iterator it;
  unsigned int counter, xstart, xend; //, xpix;
  xstart = 0;

  //Create the list of all edges of the polygon, fill and sort it
  std::list<Edge> edges;
  edgelist(polygon, ras, edges);
  edges.sort(less_by_ystart());

  // Initialize an empty list of "active" edges
  std::list<Edge> active_edges;

  //Start at the top of the first edge
  unsigned int yline(edges.front().ystart);

  //Main loop
  while(
    (yline < ras.nrow) &&
      (!(active_edges.empty() && edges.empty()))
  ) {
    // Transfer any edges starting on this row from edges to active edges
    while(edges.size() && (edges.front().ystart <= yline)) {
      active_edges.splice(active_edges.end(), edges, edges.begin());
    }
    //Sort active edges by x position of their intersection with the row
    active_edges.sort(less_by_x());

    //Iterate over active edges, fill between odd and even edges.
    counter = 0;
    for(it = active_edges.begin();
        it != active_edges.end();
        it++) {
      counter++;
      if (counter % 2) {
        xstart = ((*it).x < 0.0) ?
        0.0 :
        ((*it).x > ras.ncold ?
           ras.ncold :
           std::ceil((*it).x));
      } else {
        xend = ((*it).x < 0.0) ?
        0.0 :
        ((*it).x > ras.ncold ?
           ras.ncold :
           std::ceil((*it).x));
        record_cell(out_vector, xstart, xend, yline);
        // for(xpix = xstart; xpix < xend; ++xpix) {
        //   //note x/y switched here as raster objects store values this way
        //  // pixel_function(raster, yline, xpix, poly_value);
        //
        // }
      }
    }
    //Advance the horizontal row
    yline++;

    //Drop edges above the horizontal line, update the x-positions of the
    //intercepts of edges for the next row
    it = active_edges.begin();
    while(it != active_edges.end()) {
      if((*it).yend <= yline) {
        it = active_edges.erase(it);
      } else {
        (*it).x += (*it).dxdy;
        it++;
      }
    }
  }
}
