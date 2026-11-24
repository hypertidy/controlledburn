# Scaling: one large raster, band burn by tiles in parallel against GDAL's
# whole-raster polygons burned by controlledburn. Run by scale.sh.
#
# Rscript scale.R <step> <outdir> [cores]
#   prep:  cut the REMA window, fill nodata, write <outdir>/big.tif
#   bands: tile-parallel band burn; totals to <outdir>/bands-totals.rds
#   burn:  burn <outdir>/big-bands.gpkg and big-lines.gpkg; totals to
#          <outdir>/gdal-totals.rds
#   compare: per-band area and per-level length, both routes

args <- commandArgs(TRUE)
step <- args[1]
out <- args[2]
cores <- if (length(args) > 2) as.integer(args[3]) else 1L
here <- dirname(normalizePath(sub("^--file=", "",
                                  grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "bands.R"))
suppressPackageStartupMessages(library(controlledburn))
levels <- seq(0, 4500, by = 50)
tif <- file.path(out, "big.tif")

# totals per id: cell-areas per band (runs count 1), metres per level
totals <- function(runs, edges, lines) {
  K <- length(levels) + 1
  area <- numeric(K)
  a <- rowsum(as.numeric(runs$col_end - runs$col_start + 1L), runs$id)
  area[as.integer(rownames(a))] <- a[, 1]
  e <- rowsum(edges$fraction, edges$id)
  area[as.integer(rownames(e))] <- area[as.integer(rownames(e))] + e[, 1]
  len <- numeric(K - 1)
  l <- rowsum(lines$length, lines$id)
  len[as.integer(rownames(l))] <- l[, 1]
  list(area = area, length = len,
       records = c(runs = nrow(runs), edges = nrow(edges), lines = nrow(lines)))
}

if (step == "prep") {
  rema <- if (length(args) > 3) args[4] else raadtools::topofile("rema_100m")
  raw <- file.path(out, "big-raw.tif")
  system2("gdal_translate", c("-q", "-projwin", 0, -980700, 819300, -1800000, rema, raw))
  r <- terra::rast(raw)
  z <- terra::as.matrix(r, wide = TRUE)
  cat(sprintf("%d x %d, valid %.1f%%, range %.1f to %.1f m\n", ncol(z), nrow(z),
              100 * mean(!is.na(z)), min(z, na.rm = TRUE), max(z, na.rm = TRUE)))
  z[is.na(z)] <- -100
  terra::writeRaster(terra::rast(z, extent = terra::ext(r), crs = terra::crs(r)),
                     tif, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("TILED=YES", "COMPRESS=NONE"))
  unlink(raw)
}

if (step == "bands") {
  r <- terra::rast(tif)
  ext <- as.vector(terra::ext(r))
  dx <- terra::res(r)[1]
  dy <- terra::res(r)[2]
  nc <- ncol(r) - 1L
  nr <- nrow(r) - 1L
  ts <- 1024L
  grid <- expand.grid(i = seq(0L, nr - 1L, by = ts), j = seq(0L, nc - 1L, by = ts))
  t0 <- Sys.time()
  # each tile reads its own window (one-pixel halo) and burns it
  res <- parallel::mclapply(seq_len(nrow(grid)), function(k) {
    r0 <- grid$i[k] + 1L; r1 <- min(grid$i[k] + ts, nr)
    c0 <- grid$j[k] + 1L; c1 <- min(grid$j[k] + ts, nc)
    zt <- terra::as.matrix(terra::rast(tif)[r0:(r1 + 1L), c0:(c1 + 1L), drop = FALSE], wide = TRUE)
    et <- c(ext[1] + (c0 - 1) * dx, ext[1] + (c1 + 1) * dx,
            ext[4] - (r1 + 1) * dy, ext[4] - (r0 - 1) * dy)
    b <- band_burn(zt, et, levels, "gdal")
    b$runs$row <- b$runs$row + r0 - 1L
    b$runs$col_start <- b$runs$col_start + c0 - 1L
    b$runs$col_end <- b$runs$col_end + c0 - 1L
    b$edges$row <- b$edges$row + r0 - 1L
    b$edges$col <- b$edges$col + c0 - 1L
    b$lines$row <- b$lines$row + r0 - 1L
    b$lines$col <- b$lines$col + c0 - 1L
    b[c("runs", "edges", "lines")]
  }, mc.cores = cores, mc.preschedule = FALSE)
  runs <- do.call(rbind, lapply(res, `[[`, "runs"))
  edges <- do.call(rbind, lapply(res, `[[`, "edges"))
  lines <- do.call(rbind, lapply(res, `[[`, "lines"))
  el <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("band burn: %d tiles of %d x %d on %d cores, %.1f s; %d runs, %d edges, %d lines\n",
              nrow(grid), ts, ts, cores, el, nrow(runs), nrow(edges), nrow(lines)))
  saveRDS(totals(runs, edges, lines), file.path(out, "bands-totals.rds"))
}

if (step == "burn") {
  r <- terra::rast(tif)
  ext <- as.vector(terra::ext(r))
  dx <- terra::res(r)[1]
  dext <- c(ext[1] + dx / 2, ext[2] - dx / 2, ext[3] + dx / 2, ext[4] - dx / 2)
  dim <- c(ncol(r) - 1L, nrow(r) - 1L)
  pg <- file.path(out, "big-bands.gpkg")
  lg <- file.path(out, "big-lines.gpkg")
  t0 <- Sys.time()
  g <- vapour::vapour_read_geometry(pg)
  band <- findInterval(vapour::vapour_read_fields(pg)$lo, levels) + 1L
  l <- vapour::vapour_read_geometry(lg)
  level <- match(vapour::vapour_read_fields(lg)$elev, levels)
  t1 <- Sys.time()
  bp <- burn(g, extent = dext, dimension = dim, threads = cores)
  bl <- burn(l, extent = dext, dimension = dim, threads = cores)
  t2 <- Sys.time()
  bp$runs$id <- band[bp$runs$id]
  bp$edges$id <- band[bp$edges$id]
  bl$lines$id <- level[bl$lines$id]
  cat(sprintf("read GDAL output %.1f s; burn %.1f s on %d threads; %d runs, %d edges, %d lines\n",
              as.numeric(t1 - t0, units = "secs"), as.numeric(t2 - t1, units = "secs"),
              cores, nrow(bp$runs), nrow(bp$edges), nrow(bl$lines)))
  saveRDS(totals(bp$runs, bp$edges, bl$lines), file.path(out, "gdal-totals.rds"))
}

if (step == "compare") {
  a <- readRDS(file.path(out, "bands-totals.rds"))
  g <- readRDS(file.path(out, "gdal-totals.rds"))
  cat(sprintf("band areas (cells): total %.6g / %.6g, max |diff| per band %.3g cells\n",
              sum(a$area), sum(g$area), max(abs(a$area - g$area))))
  cat(sprintf("level lengths (m): total %.6g / %.6g, max relative diff per level %.2g\n",
              sum(a$length), sum(g$length),
              max(abs(a$length - g$length) / pmax(g$length, 1))))
  cat("records, band burn:", a$records, " GDAL route:", g$records, "\n")
}
