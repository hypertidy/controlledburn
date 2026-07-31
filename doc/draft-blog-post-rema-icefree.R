## 10m
#dsn <- "/vsicurl/https://raw.githubusercontent.com/mdsumner/rema-ovr/main/rema-vrt/10m_dem_tiles.vrt"
## 32m
dsn <- "/vsicurl/https://raw.githubusercontent.com/mdsumner/rema-ovr/main/rema-vrt/32m_dem_tiles.vrt"
shp <- "/vsizip/{/vsicurl/https://github.com/AustralianAntarcticDivision/rema.proc/raw/refs/heads/master/01_rock_classification/Medium_resolution_vector_polygons_of_Antarctic_rock_outcrop_-_VERSION_7.3.zip}/Medium resolution vector polygons of Antarctic rock outcrop - VERSION 7.3/add_rock_outcrop_medium_res_polygon_v7.3.gpkg"
rema_tiles <- "/vsizip/{/vsicurl/https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/REMA_Mosaic_Index_latest_shp.zip}/REMA_Mosaic_Index_v2_shp"
library(terra)
tiles <- terra::vect(rema_tiles)

raster_info <- vapour::vapour_raster_info(dsn)

## the crs is the same as the raster
rock_info <- vapour::vapour_layer_info(shp)
rock <- wk::wkb(vapour::vapour_read_geometry(shp), crs = rock_info$projection$Wkt)


library(controlledburn)
system.time({
  ba <- burn(rock, extent = raster_info$extent, dimension = raster_info$dimension, mode = "approx")
})
# user  system elapsed
# 0.469   0.000   0.468

system.time({
  be <- burn(rock, extent = raster_info$extent, dimension = raster_info$dimension, mode = "coverage")
})
# user  system elapsed
# 6.659   0.247   6.893


#now with a tile at
ll <- cbind(lon = 169.6539, lat = -71.50035)
(xy <- reproj::reproj_xy(ll, crs(tiles), source = "EPSG:4326"))

idx <- which(terra::is.related(tiles, vect(xy), "intersects"))
plot(tiles[idx, ])
rock_tile_index <- unlist(geos::geos_strtree_query( geos::geos_strtree(geos::as_geos_geometry(rock)), tiles[idx, ] ))
plot(rock[rock_tile_index], add = T)


tile_ex <- as.vector(ext(tiles[idx, ]))
ma <- materialize_chunk(crop_burn(ba, tile_ex))
dim(ma)
#[1] 1570 1570
ximage::ximage(ma, tile_ex, add = TRUE, col = c("black", "white"))
plot(rock[rock_tile_index], add = T, border = "hotpink", lwd = 2)


mc <- materialize_chunk(crop_burn(be, tile_ex))

ximage::ximage(mc, tile_ex, add = TRUE, col = c("black", "white"))
plot(rock[rock_tile_index], add = T, border = "hotpink", lwd = 2)

## tiny differences at this resolution
ximage::ximage(mc - ma)
