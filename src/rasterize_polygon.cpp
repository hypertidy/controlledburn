#include "edge.h"
#include "edgelist.h"
#include "rasterize_polygon.h"

// Rasterize a single polygon
// Based on https://ezekiel.encs.vancouver.wsu.edu/~cs442/lectures/rasterization/polyfill/polyfill.pdf #nolint

void record_cell(CollectorList &out_vector, unsigned int xs, unsigned int y) {
  //if (xs == xe) return;
  out_vector.push_back(Rcpp::IntegerVector::create(xs,  y));
  return;
}

void record_horiz(CollectorList &out_vector, unsigned int xs, unsigned int xe, unsigned int y) {
  if (xs == xe) return;
  int len = xe -xs;
  Rprintf("%i\n", len);
  IntegerVector out(abs(len) + 1);
  int i;
  if (len > 0) {
    for (i = 0; i < len; i++) {
      out[i] = xs + i;
    }

  } else {
    for (i = 0; i < abs(len); i++)
      out[i] = xs  -i ;
  }
  out[i] = (int)y;
  //Rprintf("%i, %i\n", out[i+1], abs(len));
 out_vector.push_back(out);
 return;
}
void rasterize_polygon(Rcpp::RObject polygon,
                       RasterInfo &ras, CollectorList &out_vector) {

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
    // note <= here, in fasterize the first row is ignored for horizontal lines
    (yline <= ras.nrow) &&
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
      if(!(*it).horiz) {
        xstart = ((*it).x < 0.0) ? 0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
        record_cell(out_vector, xstart,  yline);

    } else {
      xstart = ((*it).x < 0.0) ? 0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
      xend = xstart+(*it).dxdy;
      record_horiz(out_vector, xstart, xend, yline);
      Rprintf("horizontal edge ystart,yend: %f, %f\n", (double)(*it).ystart, (double)(*it).yend);
      Rprintf("xstart,xend: %i, %i\n", (int)xstart, (int)xend);

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
