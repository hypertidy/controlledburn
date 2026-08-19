## Write controlledburn sparse output to a SPARSE_OK GeoTIFF
##
## Follows from the CGAZ example in the README. The sparse burn result
## is written tile by tile to a SPARSE_OK GeoTIFF using gdalraster.
## Tiles with no data are never allocated on disk.
##
## This is the downstream pattern: controlledburn produces the sparse
## intersection, then any tile-oriented consumer (GeoTIFF, COG, WMTS,
## on-the-fly zonal stats) reads from it without re-burning.

library(controlledburn)
library(gdalraster)

## --- Burn (from README) --- pretty big raster of the world
g <- geos::as_geos_geometry(wk::wkb(vapour::vapour_read_geometry(sds::CGAZ())))
ba <- burn(g, dimension = c(25600L, 12800L), mode = "approx")

## Burn from REMA/icefree example - much larger, the write takes several minutes at least for ~90000 tiles
# REMA 32m DEM mosaic — we read only the grid spec, no pixel values
# dsn <- paste0(
#   "/vsicurl/https://raw.githubusercontent.com/mdsumner/rema-ovr/",
#   "main/rema-vrt/32m_dem_tiles.vrt")
# raster_info <- vapour::vapour_raster_info(dsn)
# ext <- raster_info$extent
# dm <- raster_info$dimension
#
# shp <- paste0(
#   "/vsizip/{/vsicurl/https://github.com/AustralianAntarcticDivision/",
#   "rema.proc/raw/refs/heads/master/01_rock_classification/",
#   "Medium_resolution_vector_polygons_of_Antarctic_rock_outcrop_-_",
#   "VERSION_7.3.zip}/Medium resolution vector polygons of Antarctic ",
#   "rock outcrop - VERSION 7.3/",
#   "add_rock_outcrop_medium_res_polygon_v7.3.gpkg")
#
# rock_info <- vapour::vapour_layer_info(shp)
# rock <- wk::wkb(vapour::vapour_read_geometry(shp),
#                 crs = rock_info$projection$Wkt)
#
# ba <- burn(rock, extent = ext,
#            dimension = dm, mode = "approx")

## --- Create SPARSE_OK GeoTIFF ---

outfile <- file.path(tempdir(), "rema_approx.tif")
block_size <- 512L

ext <- ba$extent
dm <- ba$dimension

ds <- create(
  format = "COG",
  dst_filename = outfile,
  xsize = dm[1], ysize = dm[2],
  nbands = 1,
  dataType = "Int32",
  options = c("TILED=YES",
              paste0("BLOCKSIZE=", block_size),
              #paste0("BLOCKYSIZE=", block_size),
              "SPARSE_OK=YES",
              "COMPRESS=DEFLATE"),
  return_obj = TRUE
)

ds$setGeoTransform(vaster::extent_dim_to_gt(ext, dm))
ds$setNoDataValue(band = 1, 0)

## --- Tile index from grout ---

gg <- grout::grout(dm, ext, blocksize = block_size)
tiles <- grout::tile_index(gg)

# Pre-filter to tile-rows that have any data — skip entire rows of
# tiles where no geometry touched any pixel row in that band.
data_rows <- unique(c(ba$runs$row, ba$edges$row, ba$lines$row, ba$points$row))
tiles <- tiles[tiles$tile_row %in% unique(ceiling(data_rows / block_size)), ]

## --- Write tiles ---
## Each row of tiles has the extent, pixel offset, and dimensions
## for one tile. crop_burn extracts the sparse data; materialize_chunk
## makes it dense; gdalraster writes the block. Empty tiles are skipped
## — SPARSE_OK means they cost nothing on disk.

tiles_written <- 0L

for (i in seq_len(nrow(tiles))) {
  sub <- crop_burn(ba, c(tiles$xmin[i], tiles$xmax[i],
                         tiles$ymin[i], tiles$ymax[i]))

  if (nrow(sub$runs) == 0 && nrow(sub$edges) == 0) next

  mat <- materialize_chunk(sub, fun = "id")
  mat[is.na(mat)] <- 0

  ds$write(band = 1,
           xoff = tiles$offset_x[i], yoff = tiles$offset_y[i],
           xsize = tiles$ncol[i], ysize = tiles$nrow[i],
           rasterData = as.integer(t(mat)))

  tiles_written <- tiles_written + 1L
}

ds$flushCache()
ds$close()

cat(sprintf("tiles written: %d of %d\n", tiles_written, nrow(tiles)))
cat(sprintf("file size: %.1f KB\n", file.size(outfile) / 1024))

## --- Verify ---

ds2 <- new(GDALRaster, outfile)
cat(sprintf("output: %d × %d %s, nodata=%s\n",
            ds2$getRasterXSize(), ds2$getRasterYSize(),
            ds2$getDataTypeName(1), ds2$getNoDataValue(1)))
ds2$close()


#library(terra)
#plot(rast(outfile))

