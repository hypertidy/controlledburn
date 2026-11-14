#!/bin/sh
# Experiment 01: GEOS GridIntersectionFractions vs the controlledburn C++ core.
# Needs GEOS >= 3.14 (GEOS_PREFIX) and the controlledburn-rs checkout (ORACLE_DIR).
#   GEOS_PREFIX=/path/to/geos-3.14 ORACLE_DIR=~/Git/controlledburn-rs/fixtures/oracle ./run.sh
set -e
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
: "${GEOS_PREFIX:?set GEOS_PREFIX to a GEOS install with the Nov 2026 GridIntersection fixes (main; 3.14.x and 3.15.0 crash on one oracle case)}"
: "${ORACLE_DIR:=$HOME/Git/controlledburn-rs/fixtures/oracle}"
build=${BUILD_DIR:-$here/build}
mkdir -p "$build"
cmake -S "$repo/cpp" -B "$build/cb" -DCMAKE_BUILD_TYPE=Release -DCONTROLLEDBURN_BUILD_TESTS=OFF >/dev/null
cmake --build "$build/cb" -j 8 >/dev/null
g++ -O2 -std=c++17 "$here/compare.cpp" -I"$repo/cpp/include" -I"$GEOS_PREFIX/include" \
  "$build/cb/libcontrolledburn.a" -L"$GEOS_PREFIX/lib" -lgeos_c -Wl,--disable-new-dtags,-rpath,"$GEOS_PREFIX/lib" -lpthread \
  -o "$build/compare"
"$build/compare" "$ORACLE_DIR" "$repo/fixtures/canon-tables.txt" "$here/differences.csv" | tee "$here/results.txt"
