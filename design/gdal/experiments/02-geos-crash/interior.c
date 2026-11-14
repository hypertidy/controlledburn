/* GEOSGridIntersectionFractions_r gives ~0 for a cell fully inside a valid
 * polygon. Usage: interior <wkb hex> xmin ymin xmax ymax nx ny row col
 * (row 0 at ymax). Prints the cell's GEOS fraction and its exact fraction
 * from GEOSIntersection of the cell box. */
#include <geos_c.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage\n"); return 2; }
    GEOSContextHandle_t h = GEOS_init_r();
    GEOSWKBReader *rd = GEOSWKBReader_create_r(h);
    GEOSGeometry *g = GEOSWKBReader_readHEX_r(h, rd, (const unsigned char *)argv[1], strlen(argv[1]));
    double x0 = atof(argv[2]), y0 = atof(argv[3]), x1 = atof(argv[4]), y1 = atof(argv[5]);
    unsigned nx = atoi(argv[6]), ny = atoi(argv[7]);
    int row = atoi(argv[8]), col = atoi(argv[9]);
    float *buf = calloc((size_t)nx * ny, sizeof(float));
    int ok = GEOSGridIntersectionFractions_r(h, g, x0, y0, x1, y1, nx, ny, buf);
    double dx = (x1 - x0) / nx, dy = (y1 - y0) / ny;
    GEOSGeometry *cell = GEOSGeom_createRectangle_r(h, x0 + col * dx, y1 - (row + 1) * dy, x0 + (col + 1) * dx, y1 - row * dy);
    GEOSGeometry *ix = GEOSIntersection_r(h, g, cell);
    double a, ca;
    GEOSArea_r(h, ix, &a);
    GEOSArea_r(h, cell, &ca);
    printf("GEOS %s valid=%d ok=%d cell(%d,%d): GridIntersectionFractions %.9g, intersection area fraction %.9g\n",
           GEOSversion(), GEOSisValid_r(h, g), ok, row, col, buf[(size_t)row * nx + col], a / ca);
    return 0;
}
