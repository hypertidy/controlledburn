#ifndef EDGE
#define EDGE

#include "stdlib.h"
#include "Rcpp.h"
using namespace Rcpp;
// A small object to hold basic info about raster dimensions
struct RasterInfo {
  double xmin, xmax, ymin, ymax, xres, yres;
  unsigned int nrow, ncol, ncold;

  RasterInfo(Rcpp::NumericVector extent, Rcpp::IntegerVector dimension) {
    xmin = extent[0];
    xmax = extent[1];
    ymin = extent[2];
    ymax = extent[3];
    ncol = dimension[0];
    nrow = dimension[1];

    ncold = ncol;
    xres = (xmax - xmin)/ncol;
    yres = (ymax - ymin)/nrow;
  }
};

// A data structure to hold only the neccessary information about a polygon
// edge needed to rasterize it
struct Edge {
  unsigned int ystart;  //the first matrix row intersected
  unsigned int yend;  //the matrix row below the end of the line
  long double dxdy; //change in x per y. Long helps with some rounding errors
  long double x; //the x location on the first matrix row intersected
  unsigned int nx;  // the x length (for horizontal lines)
  Edge(double x0, double y0, double x1, double y1, RasterInfo &ras,
       double y0c, double y1c) {
    //Convert from coordinate space to matrix row/column space. This is
    //already done for ys with y0c and y1c
    x0 = (x0 - ras.xmin)/ras.xres - 0.5; //convert from native to
    x1 = (x1 - ras.xmin)/ras.xres - 0.5; // units in the matrix
    //Rprintf("%f,%f\n", x0, x1);
    nx = 1;

    double dy;
    //Make sure edges run from top of matrix to bottom, calculate value
    if(y1c > y0c) {
      ystart = std::max(y0c, 0.0);
      dy = (y1 - y0);
      dxdy = (x1-x0)/dy;
      x = x0 + (ystart - y0)*dxdy;
      nx = (unsigned int)fabs(x0 - x1);

      yend = y1c;
    } else {
      ystart = std::max(y1c, 0.0);
      dy = (y0 - y1);
      dxdy = (x0-x1)/dy;
      x = x1 + (ystart - y1)*dxdy;
      nx = (unsigned int)fabs(x1 - x0);

      yend = y0c;
 }
   // Rprintf("%f %i\n", x1 - x0, nx);
  }
};

// These two structs allow us to compare and sort edges by their coordinates.
struct less_by_ystart {
  inline bool operator() (const Edge& struct1, const Edge& struct2) {
    return (struct1.ystart < struct2.ystart); // ||
    //((struct1.ystart == struct2.ystart) && (struct1.x < struct2.x)));
  }
};

struct less_by_x {
  inline bool operator() (const Edge& struct1, const Edge& struct2) {
    return (struct1.x < struct2.x);
  }
};

#endif
