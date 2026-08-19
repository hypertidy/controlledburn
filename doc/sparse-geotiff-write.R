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

## --- Burn (from README) ---

g <- geos::as_geos_geometry(wk::wkb(vapour::vapour_read_geometry(sds::CGAZ())))
ba <- burn(g, dimension = c(2560L, 1280L), mode = "approx")
ba

## --- Create SPARSE_OK GeoTIFF ---

outfile <- file.path(tempdir(), "cgaz_approx.tif")
block_size <- 256L

ext <- ba$extent
dm <- ba$dimension

ds <- create(
  format = "GTiff",
  dst_filename = outfile,
  xsize = dm[1], ysize = dm[2],
  nbands = 1,
  dataType = "Int32",
  options = c("TILED=YES",
              paste0("BLOCKXSIZE=", block_size),
              paste0("BLOCKYSIZE=", block_size),
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
