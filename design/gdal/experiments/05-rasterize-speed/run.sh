#!/bin/sh
# Experiment 05: GDAL rasterize vs controlledburn, same grid, same output buffer.
set -e
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../../../.." && pwd)
build=${BUILD_DIR:-$here/build}; mkdir -p "$build"
cmake -S "$repo/cpp" -B "$build/cb" -DCMAKE_BUILD_TYPE=Release -DCONTROLLEDBURN_BUILD_TESTS=OFF >/dev/null
cmake --build "$build/cb" -j 8 >/dev/null 2>&1
g++ -O2 -std=c++17 "$here/bench.cpp" -I"$repo/cpp/include" $(gdal-config --cflags) "$build/cb/libcontrolledburn.a" \
  $(gdal-config --libs) -lpthread -o "$build/bench"
"$build/bench" | tee "$here/results.txt"
