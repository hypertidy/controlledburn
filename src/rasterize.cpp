#include "edge.h"
#include "edgelist.h"
#include "rasterize.h"

// Rasterize a single polygon
// Based on https://ezekiel.encs.vancouver.wsu.edu/~cs442/lectures/rasterization/polyfill/polyfill.pdf #nolint

void record_polygon_scanline(CollectorList &out_vector, unsigned int xs, unsigned int xe, unsigned int y) {
  if (xs == xe) return;
  out_vector.push_back(Rcpp::IntegerVector::create(xs, xe - 1, y));
  return;
}
void rasterize_polygon(Rcpp::RObject polygon,
                       RasterInfo &ras, CollectorList &out_vector) {

  std::list<Edge_polygon>::iterator it;
  unsigned int counter, xstart, xend; //, xpix;
  xstart = 0;

  //Create the list of all edges of the polygon, fill and sort it
  std::list<Edge_polygon> edges;
  edgelist_polygon(polygon, ras, edges);
  edges.sort(less_by_ystart());

  // Initialize an empty list of "active" edges
  std::list<Edge_polygon> active_edges;

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
        xstart = ((*it).x < 0.0) ? 0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
      } else {
        xend = ((*it).x < 0.0) ?  0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
        record_polygon_scanline(out_vector, xstart, xend, yline);

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


void record_column_row(CollectorList &out_vector, unsigned int x, unsigned int y) {
  out_vector.push_back(Rcpp::IntegerVector::create(x, y));
  return;
}
void rasterize_line(Rcpp::RObject line,
                    RasterInfo &ras, CollectorList &out_vector) {

  std::list<Edge_line>::iterator it;
  unsigned int counter, xs, ys; //, xpix;
  xs = 0;

  //Create the list of all edges of the line, fill and sort it
  std::list<Edge_line> edges;
  edgelist_line(line, ras, edges);
  edges.sort(less_by_ystart_line());
  edges.sort(less_by_x_line());

  //Iterate over edges
  for(it = edges.begin();
      it != edges.end();
      it++) {

    for (counter = 0; counter < (*it).nmoves; counter++) {
      xs = ((*it).x < 0.0) ?  0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
      ys = ((*it).y < 0.0) ?  0.0 : ((*it).y >= ras.nrowd ? (ras.nrowd -1) : std::ceil((*it).y));
     // Rprintf("xs,ys: %i, %i\n", xs, ys);
      record_column_row(out_vector, xs, ys);
      (*it).x += (*it).dx;
      (*it).y += (*it).dy;
    }

  }

}

