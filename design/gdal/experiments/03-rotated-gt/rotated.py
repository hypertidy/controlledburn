"""Experiment 03: GDAL on a rotated geotransform.

A 40 x 30 grid of unit cells rotated 30 degrees. A polygon is rasterized with
gdal.Rasterize (cell-centre rule) and its cells counted by
`gdal raster zonal-stats --stat count` (default pixels). The truth is the
number of cell centres inside the polygon, computed with OGR from the centres'
world coordinates. GDAL transforms vertices to pixel space through the inverse
geotransform, which is exact for an affine map, so rasterize should agree.
zonal-stats builds its coverage window with the rotation terms set to zero
(alg/zonal.cpp, CalculateCoverage), so it is expected not to.
"""
import json, math, os, subprocess, sys, tempfile
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
tmp = tempfile.mkdtemp()
ncol, nrow, t = 40, 30, math.radians(30)
x0, y0 = 1000.0, 2000.0
gt = (x0, math.cos(t), math.sin(t), y0, math.sin(t), -math.cos(t))

def world(col, row):
    return (gt[0] + col * gt[1] + row * gt[2], gt[3] + col * gt[4] + row * gt[5])

srs = osr.SpatialReference(); srs.ImportFromEPSG(3857)
src = gdal.GetDriverByName("GTiff").Create(os.path.join(tmp, "src.tif"), ncol, nrow, 1, gdal.GDT_Float32)
src.SetGeoTransform(gt); src.SetSpatialRef(srs)
src.GetRasterBand(1).Fill(1.0); src = None

ring = [world(c, r) for c, r in [(3.3, 2.1), (31.7, 4.4), (36.2, 25.9), (12.5, 27.3), (3.3, 2.1)]]
wkt = "POLYGON ((" + ", ".join(f"{x!r} {y!r}" for x, y in ring) + "))"
poly = ogr.CreateGeometryFromWkt(wkt)
vec = os.path.join(tmp, "zones.gpkg")
ds = ogr.GetDriverByName("GPKG").CreateDataSource(vec)
lyr = ds.CreateLayer("zones", srs, ogr.wkbPolygon)
f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(poly); lyr.CreateFeature(f); ds = None

truth = sum(poly.Contains(ogr.CreateGeometryFromWkt("POINT (%r %r)" % world(c + 0.5, r + 0.5)))
            for r in range(nrow) for c in range(ncol))

dst = gdal.GetDriverByName("MEM").Create("", ncol, nrow, 1, gdal.GDT_Byte)
dst.SetGeoTransform(gt); dst.SetSpatialRef(srs)
gdal.Rasterize(dst, vec, burnValues=[1])
rast = int(dst.GetRasterBand(1).ReadAsArray().sum())

out = os.path.join(tmp, "zs.csv")
p = subprocess.run(["gdal", "raster", "zonal-stats", "--input", os.path.join(tmp, "src.tif"), "--zones", vec,
                    "--stat", "count", "--output", out, "--of", "CSV"], capture_output=True, text=True)
zs = open(out).read().strip().splitlines()[-1] if p.returncode == 0 else ("ERROR: " + p.stderr.strip())

print("GDAL", gdal.__version__)
print("rotated 30 deg, %d x %d unit cells, polygon area %.3f cells" % (ncol, nrow, poly.GetArea()))
print("cell centres inside the polygon (truth):", truth)
print("gdal.Rasterize cells burned:            ", rast)
print("gdal raster zonal-stats count:          ", zs)

# Control: the same polygon in index space on the unrotated grid (a rotation is
# rigid, so the truth is the same count).
gt0 = (x0, 1.0, 0.0, y0, 0.0, -1.0)
s0 = gdal.GetDriverByName("GTiff").Create(os.path.join(tmp, "src0.tif"), ncol, nrow, 1, gdal.GDT_Float32)
s0.SetGeoTransform(gt0); s0.SetSpatialRef(srs); s0.GetRasterBand(1).Fill(1.0); s0 = None
ring0 = [(x0 + c, y0 - r) for c, r in [(3.3, 2.1), (31.7, 4.4), (36.2, 25.9), (12.5, 27.3), (3.3, 2.1)]]
poly0 = ogr.CreateGeometryFromWkt("POLYGON ((" + ", ".join(f"{x!r} {y!r}" for x, y in ring0) + "))")
vec0 = os.path.join(tmp, "zones0.gpkg")
ds = ogr.GetDriverByName("GPKG").CreateDataSource(vec0)
lyr = ds.CreateLayer("zones", srs, ogr.wkbPolygon)
f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(poly0); lyr.CreateFeature(f); ds = None
out0 = os.path.join(tmp, "zs0.csv")
subprocess.run(["gdal", "raster", "zonal-stats", "--input", os.path.join(tmp, "src0.tif"), "--zones", vec0,
                "--stat", "count", "--output", out0, "--of", "CSV"], check=True, capture_output=True)
print("control, unrotated grid, zonal-stats:   ", open(out0).read().strip().splitlines()[-1])
print("zonal-stats exit status on rotated grid:", p.returncode, "stderr:", repr(p.stderr.strip()))
