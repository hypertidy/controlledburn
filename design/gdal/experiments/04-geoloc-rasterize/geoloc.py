"""Experiment 04: gdal.Rasterize onto a curvilinear (geolocation array) target.

A polar sector grid: 60 columns in angle (0 to 60 degrees), 40 rows in radius
(100 to 140), given as 2D X/Y arrays of cell centres in a projected CRS, as
GEOLOCATION metadata. Two polygons: a long sliver and a fat quadrilateral.

Truth (the cell-centre rule on this grid): cells whose centre (X[j,i], Y[j,i])
is inside the polygon.

GDAL moves only the vertices into pixel/line space through the inverse
geolocation transformer, then rasterizes straight edges in pixel space. On a
non-affine grid a straight world edge is a curve in pixel space, so long edges
cut the wrong cells. Densifying the polygon first (OGR Segmentize) is the
front end that fixes the shape.
"""
import math, os, tempfile
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
tmp = tempfile.mkdtemp()
ncol, nrow = 60, 40
th = np.radians(np.linspace(0, 60, ncol, endpoint=False) + 0.5)
rr = 140 - (np.arange(nrow) + 0.5)  # row 0 at the outer edge
X = rr[:, None] * np.cos(th[None, :])
Y = rr[:, None] * np.sin(th[None, :])

def write(name, a):
    p = os.path.join(tmp, name)
    d = gdal.GetDriverByName("GTiff").Create(p, ncol, nrow, 1, gdal.GDT_Float64)
    d.GetRasterBand(1).WriteArray(a); d = None
    return p

xp, yp = write("x.tif", X), write("y.tif", Y)
srs = osr.SpatialReference(); srs.ImportFromEPSG(3857)

def target():
    ds = gdal.GetDriverByName("MEM").Create("", ncol, nrow, 1, gdal.GDT_Byte)
    ds.SetMetadata({"SRS": srs.ExportToWkt(), "X_DATASET": xp, "X_BAND": "1", "Y_DATASET": yp, "Y_BAND": "1",
                    "PIXEL_OFFSET": "0", "LINE_OFFSET": "0", "PIXEL_STEP": "1", "LINE_STEP": "1",
                    "GEOREFERENCING_CONVENTION": "PIXEL_CENTER"}, "GEOLOCATION")
    return ds

def polar(r, deg):
    return (r * math.cos(math.radians(deg)), r * math.sin(math.radians(deg)))

polys = {
    "sliver": [polar(138, 2), polar(103, 57), polar(104.5, 58), polar(138.5, 4), polar(138, 2)],
    "quad": [polar(135, 5), polar(130, 50), polar(108, 45), polar(112, 12), polar(135, 5)],
}

def truth(g):
    inside = np.zeros((nrow, ncol), bool)
    for j in range(nrow):
        for i in range(ncol):
            inside[j, i] = g.Contains(ogr.CreateGeometryFromWkt("POINT (%r %r)" % (float(X[j, i]), float(Y[j, i]))))
    return inside

def burn(g, to=None):
    """gdal.Rasterize into the geolocated target; `to` adds transformer options."""
    ds = target()
    vec = os.path.join(tmp, "p.gpkg")
    d = ogr.GetDriverByName("GPKG").CreateDataSource(vec)
    lyr = d.CreateLayer("p", srs, ogr.wkbPolygon)
    f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(g); lyr.CreateFeature(f); d = None
    opts = "-burn 1" + "".join(" -to " + o for o in (to or []))
    gdal.Rasterize(ds, vec, options=opts)
    gdal.Unlink(vec)
    return ds.GetRasterBand(1).ReadAsArray().astype(bool)

print("GDAL", gdal.__version__)
print("polar sector, %d x %d cells, geolocation arrays of cell centres" % (ncol, nrow))

def pixel_burn(g):
    """What GDAL would do if it honoured the geolocation arrays: move each vertex
    to pixel/line through GDAL's own inverse geolocation transformer, then
    rasterize straight edges in pixel space (cell-centre rule)."""
    tr = gdal.Transformer(None, target(), ["DST_METHOD=GEOLOC_ARRAY"])
    ring = g.GetGeometryRef(0)
    pts = [ring.GetPoint_2D(k) for k in range(ring.GetPointCount())]
    px, ok = tr.TransformPoints(0, pts)
    assert all(ok), "inverse geolocation failed for a vertex"
    pg = ogr.CreateGeometryFromWkt("POLYGON ((" + ", ".join("%r %r" % (p[0], p[1]) for p in px) + "))")
    ds = gdal.GetDriverByName("MEM").Create("", ncol, nrow, 1, gdal.GDT_Byte)
    ds.SetGeoTransform((0, 1, 0, 0, 0, 1))  # pixel/line space, y down
    mem = ogr.GetDriverByName("MEM").CreateDataSource("")
    lyr = mem.CreateLayer("p", None, ogr.wkbPolygon)
    f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(pg); lyr.CreateFeature(f)
    gdal.RasterizeLayer(ds, [1], lyr, burn_values=[1])
    return ds.GetRasterBand(1).ReadAsArray().astype(bool)

tr = gdal.Transformer(None, target(), ["DST_METHOD=GEOLOC_ARRAY"])
chk = tr.TransformPoint(0, float(X[20, 30]), float(Y[20, 30]))[1]
print("check: centre of cell (row 20, col 30) maps to pixel/line (%.4f, %.4f); expected (30.5, 20.5)" % chk[:2])
print()
print("%-8s %-36s %6s %7s %7s %6s" % ("polygon", "method", "truth", "burned", "missed", "extra"))
for name, ring in polys.items():
    g = ogr.CreateGeometryFromWkt("POLYGON ((" + ", ".join(f"{x!r} {y!r}" for x, y in ring) + "))")
    t = truth(g)
    runs = [("gdal.Rasterize, geolocated target", lambda: burn(g)),
            ("gdal.Rasterize, -to DST_METHOD=GEOLOC", lambda: burn(g, ["DST_METHOD=GEOLOC_ARRAY"])),
            ("vertices moved, pixel-space burn", lambda: pixel_burn(g))]
    for step in (5.0, 1.0, 0.25):
        def f(step=step):
            d = g.Clone(); d.Segmentize(step); return pixel_burn(d)
        runs.append((f"densified to {step}, pixel-space burn", f))
    for label, fn in runs:
        b = fn()
        print("%-8s %-36s %6d %7d %7d %6d" % (name, label, t.sum(), b.sum(), (t & ~b).sum(), (b & ~t).sum()))
