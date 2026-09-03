
## -- rasterize with sedonadb
import pyarrow
import sedonadb
import numpy as np

polygon = 'POLYGON((150 150, 220 260, 190 300, 300 220, 150 150))'
ref = '<VRTDataset rasterXSize="400" rasterYSize="400"><GeoTransform>0,1,0,400,0,-1</GeoTransform><VRTRasterBand dataType="Byte" band="1"/></VRTDataset>'  
q = f"SELECT rs_asraster(st_geomfromwkt('{polygon}'), rs_frompath('{ref}'), 'B', false, 230, 0, false) AS rasterized_geom"

sd = sedonadb.connect()
d = sd.sql(q)   
d.show() 
tbl = d.to_arrow_table()
col = tbl.column("rasterized_geom").combine_chunks()

st = col.storage
dims  = st.field("spatial_dims")[0].as_py()    # e.g. ['y', 'x'] - check the order
shape = st.field("spatial_shape")[0].as_py()   # e.g. [400, 400]
tr    = st.field("transform")[0].as_py()       # affine, 6 numbers
bands = st.field("bands").values               # flattened band structs
b0    = bands.field("data")[0]
print(bands.field("data_type")[0].as_py(), bands.field("view")[0].as_py())
a = np.frombuffer(b0.as_buffer(), dtype=np.uint8).reshape(shape)
print(a.shape, (a == 230).sum())


## -- rasterize with controlledburn
import shapely
import controlledburn as cb

geom = shapely.from_wkt(polygon)
geom.wkb
r = cb.burn([geom.wkb],
            extent=(0, 400, 0, 400),   # (xmin, xmax, ymin, ymax), matches R's extent
            shape=(400, 400))          # (nrow, ncol), numpy-style

r.runs    # interior RLE:       (row, col_start, col_end, id)
r.edges   # boundary fractions: (row, col, fraction, id)
r.lines   # line lengths:       (row, col, length, id)
r.points  # point cells:        (row, col, id)

import pandas as pd
pd.DataFrame(r.edges)  # structured arrays convert directly

# Optional dense consumer (fasterize-style pixel functions):
m = cb.materialize(r, fn="sum", edge_policy="fraction")

# Tile workflow: burn once, then crop the sparse result to a sub-window and
# materialize just that tile. crop() returns a BurnResult, so it chains with
# materialize() -- and nothing dense is allocated until the final step. This is
# the Python counterpart of the R package's crop_burn() + materialize_chunk().
tile = r.crop((150, 250, 150, 250)).materialize(fn="sum", edge_policy="fraction")
#          (xmin, xmax, ymin, ymax) window, same order as bounds
# equivalently: r.crop(win) then .materialize(...) in two steps



