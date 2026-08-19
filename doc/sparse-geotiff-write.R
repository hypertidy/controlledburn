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
dim <- ba$dimension
dx <- (ext[2] - ext[1]) / dim[1]
dy <- (ext[4] - ext[3]) / dim[2]

ds <- create(
  format = "GTiff",
  dst_filename = outfile,
  xsize = dim[1], ysize = dim[2],
  nbands = 1,
  dataType = "Int32",
  options = c("TILED=YES",
              paste0("BLOCKXSIZE=", block_size),
              paste0("BLOCKYSIZE=", block_size),
              "SPARSE_OK=YES",
              "COMPRESS=DEFLATE"),
  return_obj = TRUE
)

ds$setGeoTransform(c(ext[1], dx, 0, ext[4], 0, -dy))
ds$setNoDataValue(band = 1, 0)

## --- Write tiles ---

n_tiles_x <- ceiling(dim[1] / block_size)
n_tiles_y <- ceiling(dim[2] / block_size)

# Pre-compute which rows have data for fast tile skipping
rows_with_data <- sort(unique(ba$runs$row))

tiles_written <- 0L
tiles_skipped <- 0L

for (ty in seq_len(n_tiles_y)) {
  row_lo <- (ty - 1L) * block_size + 1L
  row_hi <- min(ty * block_size, dim[2])

  # Skip entire row of tiles if no runs touch these rows
  if (!any(rows_with_data >= row_lo & rows_with_data <= row_hi)) {
    tiles_skipped <- tiles_skipped + n_tiles_x
    next
  }

  for (tx in seq_len(n_tiles_x)) {
    col_lo <- (tx - 1L) * block_size + 1L
    col_hi <- min(tx * block_size, dim[1])

    tile_xmin <- ext[1] + (col_lo - 1L) * dx
    tile_xmax <- ext[1] + col_hi * dx
    tile_ymin <- ext[4] - row_hi * dy
    tile_ymax <- ext[4] - (row_lo - 1L) * dy

    sub <- crop_burn(ba, c(tile_xmin, tile_xmax, tile_ymin, tile_ymax))

    if (nrow(sub$runs) == 0 && nrow(sub$edges) == 0) {
      tiles_skipped <- tiles_skipped + 1L
      next
    }

    # Materialise with polygon IDs (or use fun = "sum" for a mask)
    mat <- materialize_chunk(sub, fun = "id")
    mat[is.na(mat)] <- 0

    ds$write(band = 1,
             xoff = col_lo - 1L, yoff = row_lo - 1L,
             xsize = sub$dimension[1], ysize = sub$dimension[2],
             rasterData = as.integer(t(mat)))

    tiles_written <- tiles_written + 1L
  }
}

ds$flushCache()
ds$close()

cat(sprintf("tiles written: %d, skipped: %d\n", tiles_written, tiles_skipped))
cat(sprintf("file size: %.1f KB\n", file.size(outfile) / 1024))

## --- Verify ---

ds2 <- new(GDALRaster, outfile)
cat(sprintf("output: %d × %d %s, nodata=%s\n",
            ds2$getRasterXSize(), ds2$getRasterYSize(),
            ds2$getDataTypeName(1), ds2$getNoDataValue(1)))
ds2$close()
