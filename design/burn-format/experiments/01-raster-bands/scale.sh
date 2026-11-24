#!/bin/sh
# Scaling: REMA 100 m, 8193 x 8193 window (x 0 to 819 km, y -980 to -1800 km
# in EPSG:3031), 91 levels 0 to 4500 m by 50.
#
# ./scale.sh <outdir> [cores]   (needs GNU time at /usr/bin/time)
set -e
out=$1
cores=${2:-16}
here=$(dirname "$0")
mkdir -p "$out"
Rscript "$here/scale.R" prep "$out"
levels=$(Rscript -e 'cat(seq(0, 4500, 50), sep = ",")')
rm -f "$out/big-bands.gpkg" "$out/big-lines.gpkg"
echo "GDAL polygons, whole raster (one thread):"
/usr/bin/time -f "  %e s, max RSS %M KB" gdal raster contour --levels "MIN,$levels,MAX" \
  --polygonize --min-name lo --max-name hi "$out/big.tif" "$out/big-bands.gpkg" > /dev/null
echo "GDAL lines, whole raster (one thread):"
/usr/bin/time -f "  %e s, max RSS %M KB" gdal raster contour --levels "$levels" \
  --elevation-name elev "$out/big.tif" "$out/big-lines.gpkg" > /dev/null
echo "controlledburn burn of GDAL's polygons and lines:"
/usr/bin/time -f "  %e s, max RSS %M KB" Rscript "$here/scale.R" burn "$out" "$cores"
echo "band burn, tiles in parallel (R prototype):"
/usr/bin/time -f "  %e s, max RSS %M KB" Rscript "$here/scale.R" bands "$out" "$cores"
echo "band burn, one core (R prototype):"
/usr/bin/time -f "  %e s, max RSS %M KB" Rscript "$here/scale.R" bands "$out" 1
Rscript "$here/scale.R" compare "$out"
