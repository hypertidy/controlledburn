#ifndef EDGE
#define EDGE

#include "stdlib.h"
#include "Rcpp.h"
using namespace Rcpp;
// A small object to hold basic info about raster dimensions
struct RasterInfo {
  double xmin, xmax, ymin, ymax, xres, yres;
  unsigned int nrow, ncol, ncold, nrowd;

  RasterInfo(Rcpp::NumericVector extent, Rcpp::IntegerVector dimension) {
    xmin = extent[0];
    xmax = extent[1];
    ymin = extent[2];
    ymax = extent[3];
    ncol = dimension[0];
    nrow = dimension[1];

    ncold = ncol;
    nrowd = nrow;
    xres = (xmax - xmin)/ncold;
    yres = (ymax - ymin)/nrowd;

  }
};

// A data structure to hold only the neccessary information about a polygon
// edge needed to rasterize it
struct Edge_polygon {
  unsigned int ystart;  //the first matrix row intersected
  unsigned int yend;  //the matrix row below the end of the line
  long double dxdy; //change in x per y. Long helps with some rounding errors
  long double x; //the x location on the first matrix row intersected

  Edge_polygon(double x0, double y0, double x1, double y1, RasterInfo &ras,
       double y0c, double y1c) {
    //Convert from coordinate space to matrix row/column space. This is
    //already done for ys with y0c and y1c
    x0 = (x0 - ras.xmin)/ras.xres - 0.5; //convert from native to
    x1 = (x1 - ras.xmin)/ras.xres - 0.5; // units in the matrix

    //Make sure edges run from top of matrix to bottom, calculate value
    if(y1c > y0c) {
      ystart = std::max(y0c, 0.0);
      dxdy = (x1-x0)/(y1-y0);
      x = x0 + (ystart - y0)*dxdy;
      yend = y1c;
    } else {
      ystart = std::max(y1c, 0.0);
      dxdy = (x0-x1)/(y0-y1);
      x = x1 + (ystart - y1)*dxdy;
      yend = y0c;
    }
  }
};
struct Edge_line {
  long double nmoves; // larger number of steps required of x1-x0, y1-y0
  long double x; //the x location on the first matrix row intersected
  long double y;
  long double dx, dy, ystart;
  Edge_line(double x0, double y0, double x1, double y1, RasterInfo &ras) {

    // Rprintf("x: %f,%f\n", x0, x1);
    // Rprintf("y: %f,%f\n", y0, y1);
    // Rprintf("\n");
    // //Convert from coordinate space to matrix row/column space
    x0 = (x0 - ras.xmin)/ras.xres - 0.5; //convert from native to
    x1 = (x1 - ras.xmin)/ras.xres - 0.5; // units in the matrix
    y0 = (ras.ymax - y0)/ras.yres - 1.0;
    y1 = (ras.ymax - y1)/ras.yres - 1.0;
    double y0c = std::ceil(y0);
    double y1c = std::floor(y1);

    dx = (x1 - x0);
    dy = (y1 - y0);
    nmoves = std::max(std::max(std::fabs(dx), std::fabs(dy)), (long double)1.0) + 1.0;
    dx = dx/nmoves;
    dy = dy/nmoves;
    //Rprintf("x: %f,%f\n", x0, x1);
    //Rprintf("y: %f,%f\n", y0, y1);
    //Rprintf("m: %f\n\n", (double)nmoves);
    // Rprintf("%i\n", (int) nmoves);
    // Rprintf("%f\n", (double) dx);
    // Rprintf("%f\n", (double) dy);

    x = x0;
    y = y0;
    ystart = y0;
  }
};



// These two structs allow us to compare and sort edges by their coordinates.
struct less_by_ystart {
  inline bool operator() (const Edge_polygon& struct1, const Edge_polygon& struct2) {
    return (struct1.ystart < struct2.ystart); // ||
    //((struct1.ystart == struct2.ystart) && (struct1.x < struct2.x)));
  }
};

struct less_by_x {
  inline bool operator() (const Edge_polygon& struct1, const Edge_polygon& struct2) {
    return (struct1.x < struct2.x);
  }
};


struct less_by_ystart_line {
  inline bool operator() (const Edge_line& struct1, const Edge_line& struct2) {
    return (struct1.ystart < struct2.ystart); // ||
    //((struct1.ystart == struct2.ystart) && (struct1.x < struct2.x)));
  }
};

struct less_by_x_line {
  inline bool operator() (const Edge_line& struct1, const Edge_line& struct2) {
    return (struct1.x < struct2.x);
  }
};

#endif
