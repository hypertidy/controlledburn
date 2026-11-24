# Burns from a raster: checks against GDAL and isoband, tiles, pyramids, and a
# Tanaka render. See README.md.
#
# Rscript run.R <outdir> [rema_100m_dem.tif]
#
# Needs: controlledburn, terra, vapour, isoband, sf, png; GDAL >= 3.11 on the
# PATH (`gdal raster contour`). The REMA path defaults to
# raadtools::topofile("rema_100m").

args <- commandArgs(TRUE)
out <- args[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rema <- if (length(args) > 1) args[2] else raadtools::topofile("rema_100m")
here <- dirname(normalizePath(sub("^--file=", "",
                                  grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "bands.R"))
suppressPackageStartupMessages(library(controlledburn))

key <- function(d) paste(d$row, d$col, d$id)
cellsv <- function(b) {
  x <- burn_cells(b)
  data.frame(row = x$row, col = x$col, id = x$id, v = x$fraction)
}
linesv <- function(l) data.frame(row = l$row, col = l$col, id = l$id, v = l$length)
compare <- function(mine, ref, what, scale = 1) {
  a <- setNames(mine$v, key(mine))
  b <- setNames(ref$v, key(ref))
  k <- union(names(a), names(b))
  d <- abs(ifelse(is.na(a[k]), 0, a[k]) - ifelse(is.na(b[k]), 0, b[k])) / scale
  cat(sprintf("  %-34s records %8d / %8d   max |diff| %.2g   > 1e-6: %d\n",
              what, length(a), length(b), max(d), sum(d > 1e-6)))
  invisible(k[d > 1e-6])
}
gdal_bands <- function(tif, levels, gpkg) {
  unlink(gpkg)
  t <- system.time(system2("gdal", c("raster", "contour", "--levels",
    paste(c("MIN", levels, "MAX"), collapse = ","), "--polygonize",
    "--min-name", "lo", "--max-name", "hi", tif, gpkg), stdout = FALSE))
  list(geom = vapour::vapour_read_geometry(gpkg),
       band = findInterval(vapour::vapour_read_fields(gpkg)$lo, levels) + 1L,
       time = t[["elapsed"]])
}
gdal_lines <- function(tif, levels, gpkg) {
  unlink(gpkg)
  t <- system.time(system2("gdal", c("raster", "contour", "--levels",
    paste(levels, collapse = ","), "--elevation-name", "elev", tif, gpkg),
    stdout = FALSE))
  list(geom = vapour::vapour_read_geometry(gpkg),
       level = match(vapour::vapour_read_fields(gpkg)$elev, levels),
       time = t[["elapsed"]])
}
burn_ref <- function(g, ids, extent, dimension) {
  r <- burn(g, extent = extent, dimension = dimension)
  r$runs$id <- ids[r$runs$id]
  r$edges$id <- ids[r$edges$id]
  r$lines$id <- ids[r$lines$id]
  # GDAL writes one feature per band but many per level: sum the lines to one
  # record per (cell, level), as a band burn has
  if (nrow(r$lines)) r$lines <- sum_by(r$lines$row, r$lines$col, r$lines$id, r$lines$length, "length")
  r
}
isoband_ref <- function(z, extent, levels) {
  dx <- (extent[2] - extent[1]) / ncol(z)
  dy <- (extent[4] - extent[3]) / nrow(z)
  x <- extent[1] + dx * (seq_len(ncol(z)) - 0.5)
  y <- extent[4] - dy * (seq_len(nrow(z)) - 0.5)
  ib <- isoband::isobands(x, y, z, c(-1e12, levels), c(levels, 1e12))
  unclass(sf::st_as_binary(sf::st_sfc(isoband::iso_to_sfg(ib))))
}
write_tif <- function(z, extent, crs, path) {
  r <- terra::rast(z, extent = terra::ext(extent), crs = crs)
  terra::writeRaster(r, path, overwrite = TRUE, datatype = "FLT4S")
  path
}

cat("controlledburn", format(packageVersion("controlledburn")), "/",
    system2("gdal", "--version", stdout = TRUE), "/ isoband",
    format(packageVersion("isoband")), "\n\n")

## 1. volcano --------------------------------------------------------------
cat("1. volcano (87 x 61), 11 levels\n")
a <- above_level(matrix(c(0, 1, 1, 0), 1), 0.5)
b1 <- above_level(matrix(c(1, 0, 0, 0), 1), 0.5)
cat(sprintf("  hand cells: half %.6g (len %.6g); corner %.6g (len %.6g), expect 0.125 (%.6g)\n",
            a$fraction, a$length, b1$fraction, b1$length, sqrt(0.5)))
z <- volcano
ext <- c(0, ncol(z), 0, nrow(z))
tif <- write_tif(z, ext, "", file.path(out, "volcano.tif"))
for (lv in list(seq(95.3, 195.3, 10), seq(95.5, 195.5, 10))) {
  cat(sprintf("  levels %s, %s, ..., %s\n", lv[1], lv[2], lv[length(lv)]))
  b <- band_burn(z, ext, lv, "gdal")
  tot <- rowsum(burn_cells(b)$fraction, key(transform(burn_cells(b), id = 0)))
  cat(sprintf("  %-34s %d cells, max |sum - 1| %.2g\n", "partition", length(tot), max(abs(tot - 1))))
  g <- gdal_bands(tif, lv, file.path(out, "volcano-bands.gpkg"))
  compare(cellsv(b), cellsv(burn_ref(g$geom, g$band, b$extent, b$dimension)), "bands (gdal rule) vs GDAL")
  bm <- band_burn(z, ext, lv, "mean")
  compare(cellsv(bm), cellsv(burn(isoband_ref(z, ext, lv), extent = b$extent, dimension = b$dimension)),
          "bands (mean rule) vs isoband")
}

## 2. REMA 100 m, Ross Island ---------------------------------------------
cat("\n2. REMA 100 m around Mt Erebus (1025 x 1025; nodata filled with -100), 75 levels 0 to 3700 by 50\n")
src <- file.path(out, "erebus.tif")
system2("gdal_translate", c("-q", "-projwin", 251000, -1274000, 353500, -1376500, rema, src))
r <- terra::rast(src)
z <- terra::as.matrix(r, wide = TRUE)
cat(sprintf("  valid %.1f%%, range %.1f to %.1f m\n", 100 * mean(!is.na(z)),
            min(z, na.rm = TRUE), max(z, na.rm = TRUE)))
z[is.na(z)] <- -100
ext <- as.vector(terra::ext(r))
crs <- terra::crs(r)
tif <- write_tif(z, ext, crs, file.path(out, "erebus-filled.tif"))
levels <- seq(0, 3700, by = 50)
tb <- system.time(b <- band_burn(z, ext, levels, "gdal"))[["elapsed"]]
cat(sprintf("  band burn: %d runs, %d edges, %d lines; %d saddle (cell, level) pairs; %d negative fractions; %.2f s (R)\n",
            nrow(b$runs), nrow(b$edges), nrow(b$lines), nrow(b$notes$saddles),
            b$notes$negative_fractions, tb))
tot <- rowsum(burn_cells(b)$fraction, key(transform(burn_cells(b), id = 0)))
cat(sprintf("  %-34s %d cells, max |sum - 1| %.2g\n", "partition", length(tot), max(abs(tot - 1))))
gb <- gdal_bands(tif, levels, file.path(out, "erebus-bands.gpkg"))
tr <- system.time(ref <- burn_ref(gb$geom, gb$band, b$extent, b$dimension))[["elapsed"]]
cat(sprintf("  GDAL polygons %.2f s, their burn %.2f s\n", gb$time, tr))
bad <- compare(cellsv(b), cellsv(ref), "bands (gdal rule) vs GDAL")
gl <- gdal_lines(tif, levels, file.path(out, "erebus-lines.gpkg"))
lref <- burn_ref(gl$geom, gl$level, b$extent, b$dimension)
compare(linesv(b$lines), linesv(lref$lines), "lines vs GDAL lines (per cell size)", scale = 100)
ib <- isoband_ref(z, ext, levels)
iref <- burn(ib, extent = b$extent, dimension = b$dimension)
bm <- band_burn(z, ext, levels, "mean")
cat(sprintf("  mean rule: %d negative fractions\n", bm$notes$negative_fractions))
compare(cellsv(bm), cellsv(iref), "bands (mean rule) vs isoband")
diffs <- compare(cellsv(b), cellsv(iref), "bands (gdal rule) vs isoband (saddle rules differ)")
sc <- unique(paste(b$notes$saddles$row, b$notes$saddles$col))
dc <- unique(vapply(strsplit(diffs, " "), function(s) paste(s[1:2], collapse = " "), ""))
cat(sprintf("  ... those differences are in %d cells, all saddle cells: %s\n",
            length(dc), all(dc %in% sc)))

## 3. Tiles -----------------------------------------------------------------
cat("\n3. Tiles: 4 x 4 tiles of the dual grid, each burned from its own window (one-pixel halo)\n")
dx <- (ext[2] - ext[1]) / ncol(z)
dy <- (ext[4] - ext[3]) / nrow(z)
edges <- seq(0, b$dimension[1], length.out = 5)
tiles <- list()
for (i in 1:4) for (j in 1:4) {
  r0 <- edges[i] + 1; r1 <- edges[i + 1]
  c0 <- edges[j] + 1; c1 <- edges[j + 1]
  zt <- z[r0:(r1 + 1), c0:(c1 + 1)]
  et <- c(ext[1] + (c0 - 1) * dx, ext[1] + (c1 + 1) * dx,
          ext[4] - (r1 + 1) * dy, ext[4] - (r0 - 1) * dy)
  bt <- band_burn(zt, et, levels, "gdal")
  ct <- cellsv(bt); ct$row <- ct$row + r0 - 1; ct$col <- ct$col + c0 - 1
  lt <- linesv(bt$lines); lt$row <- lt$row + r0 - 1; lt$col <- lt$col + c0 - 1
  tiles[[length(tiles) + 1]] <- list(c = ct, l = lt)
}
compare(do.call(rbind, lapply(tiles, `[[`, "c")), cellsv(b), "tiles vs whole: bands")
compare(do.call(rbind, lapply(tiles, `[[`, "l")), linesv(b$lines), "tiles vs whole: lines")

## 4. Pyramid ---------------------------------------------------------------
cat("\n4. Pyramid: group-by of the 100 m burn against a burn of GDAL's polygons and lines onto the coarse grid\n")
for (f in c(2, 4, 8)) {
  ta <- system.time(bf <- aggregate_burn(b, f))[["elapsed"]]
  cat(sprintf("  factor %d (%d x %d, %.2f s):\n", f, bf$dimension[1], bf$dimension[2], ta))
  compare(cellsv(bf), cellsv(burn_ref(gb$geom, gb$band, b$extent, bf$dimension)), "bands")
  compare(linesv(bf$lines), linesv(burn_ref(gl$geom, gl$level, b$extent, bf$dimension)$lines),
          "lines (per cell size)", scale = 100 * f)
}

## 5. Continuity at a saddle ---------------------------------------------------
cat("\n5. One cell from scale.sh (row 906, col 1264 of the large window), level 0, as its upper-left corner crosses the level\n")
bilinear_above <- function(z, lev, n = 4000) {
  u <- (seq_len(n) - 0.5) / n
  g <- expand.grid(u = u, v = u)
  f <- with(g, z[1] * (1 - u) * (1 - v) + z[2] * u * (1 - v) + z[3] * u * v + z[4] * (1 - u) * v)
  mean(f >= lev)
}
for (ul in c(-1e-3, -2.38418573772e-09, 2.38418573772e-09, 1e-3)) {
  zc <- c(-1.07340037823e-01, 0.0839747786522, -0.0348985232413, ul) # LL, LR, UR, UL
  cat(sprintf("  UL %+.3g: above 0, linear (gdal rule) %.4f, linear (mean rule) %.4f, bilinear (4000^2 samples) %.4f\n",
              ul, above_level(matrix(zc, 1), 0, "gdal")$fraction,
              above_level(matrix(zc, 1), 0, "mean")$fraction, bilinear_above(zc, 0)))
}

## 6. Tanaka ----------------------------------------------------------------
cat("\n6. Tanaka contours\n")
w <- tanaka_weights(b, z, ext, azimuth = 315)
cat(sprintf("  weights: %d lit, %d shaded records\n", nrow(w$lit), nrow(w$shade)))
palette_for <- function(lv) {
  # band 1 is below sea level (the filled ocean); the rest a narrow mid-grey
  # ramp, as Tanaka's grey ground
  c("#2c3e50", grDevices::gray(seq(0.42, 0.62, length.out = length(lv))))
}
# native: 100 m interval from the 50 m burn, by group-by on id; strokes
# 90 m wide at |cos theta| = 1
n100 <- merge_levels(b, seq(0, 3700, 100), weights = w)
img <- render_tanaka(n100, n100$weights, palette_for(n100$source$levels), width = 90)
png::writePNG(img, file.path(out, "tanaka-100m.png"))
png::writePNG(img[300:811, 300:811, ], file.path(out, "tanaka-100m-erebus.png"))
# coarse: factor 4 (400 m cells) by group-by, then a 200 m interval by
# group-by on id; strokes scaled with the cell
c4 <- aggregate_burn(b, 4, weights = w)
c4 <- merge_levels(c4, seq(0, 3600, 200), weights = c4$weights)
img4 <- render_tanaka(c4, c4$weights, palette_for(c4$source$levels), width = 360)
png::writePNG(img4, file.path(out, "tanaka-400m.png"))
cat(sprintf("  wrote tanaka-100m.png (%d x %d), tanaka-100m-erebus.png (512 x 512 crop) and tanaka-400m.png (%d x %d) to %s\n",
            n100$dimension[1], n100$dimension[2], c4$dimension[1], c4$dimension[2], out))
