// predicates.h -- exact orientation test, shared by the approx sweep
// (scanline.cpp) and relate_cells (relate.cpp). Internal to the core.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_PREDICATES_H
#define CONTROLLEDBURN_PREDICATES_H

#include <cmath>

namespace controlledburn {
namespace detail {

// a + b as the rounded sum s and its exact error e.
inline void two_sum(double a, double b, double& s, double& e) {
  s = a + b;
  double bb = s - a;
  e = (a - (s - bb)) + (b - bb);
}

// Sign of (bx - ax)(cy - ay) - (cx - ax)(by - ay), exactly: positive when c
// is left of the line a -> b. The determinant is expanded into six products
// of input coordinates, each split exactly with an FMA, and the twelve terms
// are summed as a non-overlapping expansion (Shewchuk's grow-expansion),
// whose sign is that of its largest component. A float filter answers the
// common case. Same as orient2d_sign in controlledburn-rs.
inline int orient2d_sign(double ax, double ay, double bx, double by,
                         double cx, double cy) {
  double l = (bx - ax) * (cy - ay);
  double r = (cx - ax) * (by - ay);
  double det = l - r;
  double bound = 4.0e-16 * (std::fabs(l) + std::fabs(r));  // ccwerrboundA, padded
  if (det > bound) return 1;
  if (det < -bound) return -1;
  const double prods[6][3] = {
    {bx, cy, 1.0}, {bx, ay, -1.0}, {ax, cy, -1.0},
    {cx, by, -1.0}, {cx, ay, 1.0}, {ax, by, 1.0}
  };
  double e[12];
  int n = 0;
  for (const auto& pq : prods) {
    double hi = pq[0] * pq[1];
    double lo = std::fma(pq[0], pq[1], -hi);
    const double terms[2] = {pq[2] * lo, pq[2] * hi};
    for (double t : terms) {
      // grow-expansion: add t to e[0..n), keeping increasing magnitude
      double q = t;
      for (int k = 0; k < n; k++) {
        double sum, err;
        two_sum(q, e[k], sum, err);
        e[k] = err;
        q = sum;
      }
      e[n++] = q;
    }
  }
  for (int k = n - 1; k >= 0; k--) {
    if (e[k] > 0.0) return 1;
    if (e[k] < 0.0) return -1;
  }
  return 0;
}

} // namespace detail
} // namespace controlledburn

#endif // CONTROLLEDBURN_PREDICATES_H
