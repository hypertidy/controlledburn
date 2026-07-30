#!/bin/sh
# sync-core.sh -- sync the canonical C++ core (cpp/) into the R package
# build locations. Run from the repository root after changing anything
# under cpp/. The copies are committed; this script is the single source
# of truth for how they are derived.
#
#   cpp/include/controlledburn/*.hpp -> inst/include/controlledburn/
#       (public C++ API; also available to downstream R packages via
#        LinkingTo: controlledburn)
#   cpp/src/scanline.cpp -> src/core_scanline.cpp
#   cpp/src/wkb.cpp      -> src/core_wkb.cpp
#
# Include-path rewrites: the standalone core vendors its GEOS-free
# exactextract subset under cpp/src/ee/, while the R package already
# compiles the full set under src/exactextract/ (needed by the older
# burn_sparse path). The rewrite points the synced core at the existing
# objects so each exactextract translation unit is compiled exactly once:
#   "ee/analytical_coverage.h" -> "analytical_coverage.h"  (src/ root)
#   "ee/<file>"                -> "exactextract/<file>"
#
# cpp/src/ee/analytical_coverage.h and src/analytical_coverage.h must be
# kept in sync by hand (they differ only in include paths); the diff
# check at the bottom enforces this.

set -e
cd "$(dirname "$0")/.."

mkdir -p inst/include/controlledburn
cp cpp/include/controlledburn/*.hpp inst/include/controlledburn/

sed -e 's|#include "ee/analytical_coverage.h"|#include "analytical_coverage.h"|' \
    -e 's|#include "ee/|#include "exactextract/|' \
    cpp/src/scanline.cpp > src/core_scanline.cpp

cp cpp/src/wkb.cpp src/core_wkb.cpp

# Enforce analytical_coverage.h parity (modulo the include-path lines)
a=$(sed 's|#include "exactextract/|#include "|' src/analytical_coverage.h)
b=$(cat cpp/src/ee/analytical_coverage.h)
if [ "$a" != "$b" ]; then
  echo "WARNING: src/analytical_coverage.h and cpp/src/ee/analytical_coverage.h have diverged" >&2
  exit 1
fi

echo "core synced: inst/include/controlledburn/, src/core_scanline.cpp, src/core_wkb.cpp"
