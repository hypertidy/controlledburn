#include "edge.h"
#include "edgelist.h"
#include "rasterize_polygon.h"

// Rasterize a single polygon
// Based on https://ezekiel.encs.vancouver.wsu.edu/~cs442/lectures/rasterization/polyfill/polyfill.pdf #nolint

void record_cell_line(CollectorList &out_vector, unsigned int xs, unsigned int y) {
  out_vector.push_back(Rcpp::IntegerVector::create(xs, y));
  return;
}
void rasterize_line(Rcpp::RObject polygon,
                    RasterInfo &ras, CollectorList &out_vector) {

  std::list<Edge>::iterator it;
  unsigned int  xstart, xend, xpix; //,counter;
  xstart = 0;

  //Create the list of all edges of the polygon, fill and sort it
  std::list<Edge> edges;
  edgelist(polygon, ras, edges);
  edges.sort(less_by_ystart());

  // Initialize an empty list of "active" edges
  std::list<Edge> active_edges;

  //Start at the top of the first edge
  unsigned int yline (edges.front().ystart);


  //Main loop
  while(
    (yline < ras.nrow) &&
      (!(active_edges.empty() && edges.empty()))
  ) {
    // Transfer any edges starting on this row from edges to active edges
    // HACK: yline + 1 else horizontal lines don't get included
    while(edges.size() && (edges.front().ystart <= yline)) {
      active_edges.splice(active_edges.end(), edges, edges.begin());
    }


    //Sort active edges by x position of their intersection with the row
    active_edges.sort(less_by_x());


    //Iterate over active edges
    for(it = active_edges.begin();
        it != active_edges.end();
        it++) {
      xstart = ((*it).x < 0.0) ? 0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
      record_cell_line(out_vector, xstart, yline);



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







  xstart = 0;

  //Create the list of all edges of the polygon, fill and sort it
  std::list<Edge> edges_h;
  edgelist_line(polygon, ras, edges_h);
  edges_h.sort(less_by_ystart());

  // Initialize an empty list of "active" edges
  std::list<Edge> active_edges_h;

  //Start at the top of the first edge
  yline =  edges_h.front().ystart;


  //Main loop
  while(
    (yline < ras.nrow) &&
      (!(active_edges_h.empty() && edges_h.empty()))
  ) {
    // Transfer any edges starting on this row from edges to active edges

    while(edges_h.size() && (edges_h.front().ystart <= yline)) {
      active_edges_h.splice(active_edges_h.end(), edges_h, edges_h.begin());
    }


    //Sort active edges by x position of their intersection with the row
    active_edges_h.sort(less_by_x());


    //Iterate over active edges
    for(it = active_edges_h.begin();
        it != active_edges_h.end();
        it++) {
      xstart = ((*it).x < 0.0) ? 0.0 : ((*it).x >= ras.ncold ? (ras.ncold -1) : std::ceil((*it).x));
      for (xpix = 0; xpix < (*it).nx; xpix++) {
        record_cell_line(out_vector, xstart + xpix, yline);
        //if (xpix > 1000) break;
      }


    }
    //Advance the horizontal row
    yline++;


    //Drop edges above the horizontal line, update the x-positions of the
    //intercepts of edges for the next row
    it = active_edges_h.begin();
    while(it != active_edges_h.end()) {
      if((*it).yend <= yline) {
        it = active_edges_h.erase(it);
      } else {
        //(*it).x += (*it).dxdy;
        it++;
      }
    }
  }
}
