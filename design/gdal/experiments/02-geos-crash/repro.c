/* GEOSGridIntersectionFractions_r crashes (SIGSEGV in Cell::traversal_in_progress)
 * on a valid triangle whose vertices sit near, not on, the lattice.
 * cc repro.c -I$GEOS/include -L$GEOS/lib -lgeos_c -o repro */
#include <geos_c.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    const char *wkt = argc > 1 ? argv[1] : "POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))";
    unsigned n = argc > 2 ? (unsigned)atoi(argv[2]) : 100;
    GEOSContextHandle_t h = GEOS_init_r();
    printf("GEOS %s\n%s on [0,10]x[0,10], %u x %u cells\n", GEOSversion(), wkt, n, n);
    fflush(stdout);
    GEOSWKTReader *rd = GEOSWKTReader_create_r(h);
    GEOSGeometry *g = GEOSWKTReader_read_r(h, rd, wkt);
    printf("valid: %d\n", GEOSisValid_r(h, g));
    fflush(stdout);
    float *buf = calloc((size_t)n * n, sizeof(float));
    int ok = GEOSGridIntersectionFractions_r(h, g, 0, 0, 10, 10, n, n, buf);
    double s = 0;
    for (size_t i = 0; i < (size_t)n * n; i++) s += buf[i];
    double area;
    GEOSArea_r(h, g, &area);
    printf("ok %d, sum of fractions * cell area %.12g, polygon area %.12g\n", ok, s * 100.0 / ((double)n * n), area);
    free(buf);
    GEOSGeom_destroy_r(h, g);
    GEOSWKTReader_destroy_r(h, rd);
    GEOS_finish_r(h);
    return 0;
}
