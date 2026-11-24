# Burns from a raster, cell by cell: isobands as runs + edges, isolines as
# lines, Tanaka weights, and pyramid levels by group-by.
#
# Everything here is local to one cell of the dual grid (the cell whose four
# corners are four neighbouring pixel centres). Nothing joins segments into
# lines or nests rings, so any tiling of the raster (with a one-pixel halo)
# gives the same tables.
#
# Conventions follow controlledburn's R tables: 1-based row and col, row 1 at
# the top, run col_end inclusive.
#
# Interpolant: linear on the cell walls with straight chords between the wall
# crossings, as marching squares (GDAL's contouring, isoband). Saddle cells
# have two readings, chosen by `saddle`:
#   "gdal": always cut off the lower-left and upper-right corners, as GDAL's
#           alg/marching_squares/square.h does;
#   "mean": join the high corners when the mean of the four corners is at or
#           above the level, otherwise the low ones.

# Area (as a fraction of the cell) of the part of each cell where the linear
# interpolant is above `lev` (keep = "high") or below it (keep = "low"), and
# the length of the isoline chords, in CRS units.
#
# Z is an n x 4 matrix of corner values in counter-clockwise order:
# lower-left, lower-right, upper-right, upper-left. lev is a vector of n.
#
# The clip walks the cell's perimeter (Sutherland-Hodgman against a scalar
# threshold): each corner kept, then each wall crossing, up to 8 vertices.
# Slots with no vertex take a copy of the previous vertex, which adds nothing
# to the shoelace sum or the perimeter. The isoline length is the clipped
# polygon's perimeter less the length of cell wall it keeps.
clip_cell <- function(Z, lev, keep = c("high", "low"), dx = 1, dy = 1) {
  keep <- match.arg(keep)
  n <- nrow(Z)
  H <- if (keep == "high") Z >= lev else Z < lev
  px <- c(0, dx, dx, 0)
  py <- c(0, 0, dy, dy)
  elen <- c(dx, dy, dx, dy)
  X <- matrix(NA_real_, n, 8)
  Y <- matrix(NA_real_, n, 8)
  wall <- numeric(n)
  for (k in 1:4) {
    k2 <- k %% 4L + 1L
    hk <- H[, k]
    hk2 <- H[, k2]
    X[hk, 2 * k - 1] <- px[k]
    Y[hk, 2 * k - 1] <- py[k]
    cr <- hk != hk2
    t <- (lev[cr] - Z[cr, k]) / (Z[cr, k2] - Z[cr, k])
    X[cr, 2 * k] <- px[k] + t * (px[k2] - px[k])
    Y[cr, 2 * k] <- py[k] + t * (py[k2] - py[k])
    w <- numeric(n)
    w[hk & hk2] <- elen[k]
    out <- hk & !hk2
    w[out] <- t[out[cr]] * elen[k]
    inn <- !hk & hk2
    w[inn] <- (1 - t[inn[cr]]) * elen[k]
    wall <- wall + w
  }
  for (s in 2:8) {
    na <- is.na(X[, s])
    X[na, s] <- X[na, s - 1]
    Y[na, s] <- Y[na, s - 1]
  }
  for (s in 1:7) {
    na <- is.na(X[, s])
    X[na, s] <- X[na, 8]
    Y[na, s] <- Y[na, 8]
  }
  X[is.na(X)] <- 0
  Y[is.na(Y)] <- 0
  s2 <- c(2:8, 1)
  area <- 0.5 * rowSums(X * Y[, s2] - X[, s2] * Y)
  perim <- rowSums(sqrt((X[, s2] - X)^2 + (Y[, s2] - Y)^2))
  list(fraction = area / (dx * dy), length = pmax(perim - wall, 0))
}

# Fraction of each cell above `lev`, and the isoline length, with saddle
# cells resolved by `saddle`.
above_level <- function(Z, lev, saddle = c("gdal", "mean"), dx = 1, dy = 1) {
  saddle <- match.arg(saddle)
  hi <- clip_cell(Z, lev, "high", dx, dy)
  H <- Z >= lev
  sad <- (H[, 1] == H[, 3]) & (H[, 2] == H[, 4]) & (H[, 1] != H[, 2])
  # The "high" clip joins the high corners. Where the rule says join the low
  # corners instead, use the complement of the "low" clip.
  join_low <- if (saddle == "gdal") {
    sad & H[, 1] # lower-left and upper-right high: GDAL cuts them off
  } else {
    sad & (rowMeans(Z) < lev)
  }
  if (any(join_low)) {
    lo <- clip_cell(Z[join_low, , drop = FALSE], lev[join_low], "low", dx, dy)
    hi$fraction[join_low] <- 1 - lo$fraction
    hi$length[join_low] <- lo$length
  }
  hi$saddle <- sad
  hi
}

# The dual grid of a raster: corners at pixel centres.
# z: matrix, row 1 at the top. extent: c(xmin, xmax, ymin, ymax) of the pixels.
dual_grid <- function(z, extent) {
  nr <- nrow(z)
  nc <- ncol(z)
  dx <- (extent[2] - extent[1]) / nc
  dy <- (extent[4] - extent[3]) / nr
  list(
    extent = c(extent[1] + dx / 2, extent[2] - dx / 2,
               extent[3] + dy / 2, extent[4] - dy / 2),
    dimension = c(nc - 1L, nr - 1L),
    dx = dx, dy = dy,
    # corners as (nr - 1) x (nc - 1) matrices, counter-clockwise from lower-left
    LL = z[-1, -nc], LR = z[-1, -1], UR = z[-nr, -1], UL = z[-nr, -nc]
  )
}

# Round to the nearest f32, as the crate stores its measures.
f32 <- function(x) readBin(writeBin(x, raw(), size = 4), "numeric", size = 4, n = length(x))

# Run-length encode cells (row, col, id) along rows.
encode_runs <- function(row, col, id) {
  if (!length(row)) {
    return(data.frame(row = integer(), col_start = integer(),
                      col_end = integer(), id = integer()))
  }
  o <- order(row, col)
  row <- row[o]; col <- col[o]; id <- id[o]
  n <- length(row)
  start <- c(TRUE, row[-1] != row[-n] | col[-1] != col[-n] + 1L | id[-1] != id[-n])
  g <- cumsum(start)
  data.frame(
    row = row[start],
    col_start = col[start],
    col_end = as.integer(tapply(col, g, max)),
    id = id[start]
  )
}

# The band burn of raster z on its dual grid.
#
# levels: sorted finite levels L1 < ... < LK.
# Band ids: 1 is (-Inf, L1), j + 1 is [Lj, Lj+1), K + 1 is [LK, Inf).
# Line ids: j is the isoline at Lj.
# Cells with a missing corner are skipped (counted in $notes).
band_burn <- function(z, extent, levels, saddle = c("gdal", "mean"),
                      threshold = 1e-6) {
  saddle <- match.arg(saddle)
  stopifnot(!is.unsorted(levels, strictly = TRUE))
  g <- dual_grid(z, extent)
  nr <- g$dimension[2]
  K <- length(levels)
  Zall <- cbind(as.vector(g$LL), as.vector(g$LR), as.vector(g$UR), as.vector(g$UL))
  ok <- stats::complete.cases(Zall)
  cell <- which(ok)
  Zall <- Zall[ok, , drop = FALSE]
  row <- (cell - 1L) %% nr + 1L
  col <- (cell - 1L) %/% nr + 1L
  cmin <- do.call(pmin, as.data.frame(Zall))
  cmax <- do.call(pmax, as.data.frame(Zall))
  bmin <- findInterval(cmin, levels)
  bmax <- findInterval(cmax, levels)

  # Whole cells: every corner in one band.
  full <- bmin == bmax

  # Mixed cells: one (cell, level) pair per level the cell straddles.
  m <- which(!full)
  nlev <- bmax[m] - bmin[m]
  pc <- rep(m, nlev)
  pk <- bmin[pc] + sequence(nlev)
  a <- above_level(Zall[pc, , drop = FALSE], levels[pk], saddle, g$dx, g$dy)

  # Band j of a mixed cell has fraction A(j - 1) - A(j), with A(bmin) = 1 and
  # A(bmax + 1) = 0.
  key <- c((pc - 1) * (K + 1) + pk + 1, (pc - 1) * (K + 1) + pk,
           (m - 1) * (K + 1) + bmin[m] + 1)
  val <- c(a$fraction, -a$fraction, rep(1, length(m)))
  s <- rowsum(val, key, reorder = TRUE)
  skey <- as.numeric(rownames(s))
  frac <- s[, 1]
  ecell <- (skey - 1) %/% (K + 1) + 1
  eband <- as.integer((skey - 1) %% (K + 1) + 1)
  negative <- sum(frac < -1e-12)
  # Records are f32, as the crate's are, and the crate's threshold is
  # symmetric and applied to the f32 value: a fraction of `threshold` or less
  # is dropped, and one of 1 - threshold or more is a whole cell (a run).
  frac <- f32(frac)
  whole <- frac >= f32(1 - f32(threshold))
  keep <- frac > f32(threshold) & !whole
  runs <- encode_runs(c(row[full], row[ecell[whole]]),
                      c(col[full], col[ecell[whole]]),
                      c(bmin[full] + 1L, eband[whole]))
  edges <- data.frame(row = row[ecell[keep]], col = col[ecell[keep]],
                      fraction = frac[keep], id = eband[keep])
  edges <- edges[order(edges$row, edges$col, edges$id), ]
  rownames(edges) <- NULL

  lk <- a$length > 0
  lines <- data.frame(row = row[pc[lk]], col = col[pc[lk]],
                      length = f32(a$length[lk]), id = pk[lk])
  lines <- lines[order(lines$row, lines$col, lines$id), ]
  rownames(lines) <- NULL

  structure(list(
    runs = runs, edges = edges, lines = lines,
    points = data.frame(row = integer(), col = integer(), id = integer()),
    extent = g$extent, dimension = g$dimension,
    options = list(mode = "coverage", merge_parts = TRUE,
                   coalesce_runs = TRUE, id_offset = 0L),
    source = list(levels = levels, saddle = saddle, interpolant = "linear",
                  grid = "dual"),
    notes = list(skipped_cells = sum(!ok), negative_fractions = negative,
                 saddles = data.frame(row = row[pc[a$saddle]],
                                      col = col[pc[a$saddle]],
                                      level = pk[a$saddle]))
  ), class = "controlledburn")
}

# Tanaka weights for the lines of a band burn: length times cos(theta), split
# into lit (cos > 0) and shaded parts. theta is between the downhill direction
# of the cell's mean gradient (from its four corners) and the direction
# towards the light (azimuth in degrees clockwise from north).
tanaka_weights <- function(b, z, extent, azimuth = 315) {
  g <- dual_grid(z, extent)
  nr <- g$dimension[2]
  i <- (b$lines$col - 1L) * nr + b$lines$row
  gx <- ((g$LR[i] + g$UR[i]) - (g$LL[i] + g$UL[i])) / (2 * g$dx)
  gy <- ((g$UR[i] + g$UL[i]) - (g$LL[i] + g$LR[i])) / (2 * g$dy)
  gn <- sqrt(gx^2 + gy^2)
  az <- azimuth * pi / 180
  cs <- ifelse(gn > 0, (-gx * sin(az) - gy * cos(az)) / gn, 0)
  w <- b$lines
  lit <- data.frame(row = w$row, col = w$col, weight = w$length * pmax(cs, 0), id = w$id)
  shade <- data.frame(row = w$row, col = w$col, weight = w$length * pmax(-cs, 0), id = w$id)
  list(lit = lit[lit$weight > 0, ], shade = shade[shade$weight > 0, ])
}

# Every polygon cell as (row, col, id, fraction): runs expanded, plus edges.
burn_cells <- function(b) {
  r <- b$runs
  n <- r$col_end - r$col_start + 1L
  rbind(
    data.frame(row = rep(r$row, n),
               col = sequence(n, from = r$col_start),
               id = rep(r$id, n), fraction = 1),
    b$edges[, c("row", "col", "id", "fraction")]
  )
}

# Sum v by (row, col, id), sorted by row, col, id.
sum_by <- function(row, col, id, v, name = "v") {
  if (!length(row)) {
    d <- data.frame(row = integer(), col = integer(), v = numeric(), id = integer())
    names(d)[3] <- name
    return(d)
  }
  nc <- max(col) + 1
  ni <- max(id) + 1
  key <- (as.numeric(row) * nc + col) * ni + id
  s <- rowsum(v, key, reorder = TRUE)
  k <- as.numeric(rownames(s))
  d <- data.frame(row = as.integer(k %/% (nc * ni)),
                  col = as.integer((k %/% ni) %% nc),
                  v = s[, 1],
                  id = as.integer(k %% ni))
  names(d)[3] <- name
  d
}

# Polygon cells (row, col, id, fraction) back to runs and edges, with the
# crate's symmetric threshold.
cells_to_burn <- function(b, cells, threshold = 1e-6) {
  cells$fraction <- f32(cells$fraction)
  whole <- cells$fraction >= f32(1 - f32(threshold))
  b$runs <- encode_runs(cells$row[whole], cells$col[whole], cells$id[whole])
  e <- cells[!whole & cells$fraction > f32(threshold), ]
  b$edges <- e[, c("row", "col", "fraction", "id")]
  rownames(b$edges) <- NULL
  b
}

# A coarser level by group-by: integer factor f, grid dimensions divisible by f.
# Polygon cells: fraction is the sum of fine fractions over f^2. Lines and
# any weight tables (a list of row, col, weight, id): sums.
aggregate_burn <- function(b, f, weights = NULL, threshold = 1e-6) {
  stopifnot(all(b$dimension %% f == 0))
  cr <- function(x) (x - 1L) %/% f + 1L
  cells <- burn_cells(b)
  s <- sum_by(cr(cells$row), cr(cells$col), cells$id, cells$fraction / f^2, "fraction")
  out <- cells_to_burn(b, s, threshold)
  out$lines <- sum_by(cr(b$lines$row), cr(b$lines$col), b$lines$id, b$lines$length, "length")
  out$dimension <- b$dimension %/% f
  if (!is.null(weights)) {
    out$weights <- lapply(weights, function(t) sum_by(cr(t$row), cr(t$col), t$id, t$weight, "weight"))
  }
  out
}

# A coarser contour interval by group-by on id: keep the levels `keep` (a
# subset of the burn's levels). Bands merge, lines and weights are filtered.
merge_levels <- function(b, keep, weights = NULL, threshold = 1e-6) {
  L <- b$source$levels
  ki <- match(keep, L)
  stopifnot(!anyNA(ki))
  # fine band j is [L[j - 1], L[j]); its new band is the kept interval
  # holding L[j - 1]
  lower <- c(-Inf, L)
  newband <- findInterval(lower, keep) + 1L
  cells <- burn_cells(b)
  s <- sum_by(cells$row, cells$col, newband[cells$id], cells$fraction, "fraction")
  out <- cells_to_burn(b, s, threshold)
  relabel <- function(t) {
    t <- t[t$id %in% ki, ]
    t$id <- match(t$id, ki)
    t
  }
  out$lines <- relabel(b$lines)
  out$source$levels <- keep
  if (!is.null(weights)) out$weights <- lapply(weights, relabel)
  out
}

# Render a band burn (tint) and Tanaka weights (strokes) to an RGB array.
# palette: one colour per band. width: stroke width in CRS units at
# |cos theta| = 1, so the ink in a cell is weight * width / cell area.
render_tanaka <- function(b, weights, palette, width,
                          lit_col = c(1, 1, 1), shade_col = c(0.1, 0.1, 0.15)) {
  nc <- b$dimension[1]
  nr <- b$dimension[2]
  area <- ((b$extent[2] - b$extent[1]) / nc) * ((b$extent[4] - b$extent[3]) / nr)
  pal <- grDevices::col2rgb(palette) / 255
  cells <- burn_cells(b)
  idx <- (cells$col - 1L) * nr + cells$row
  img <- array(0, c(nr, nc, 3))
  for (ch in 1:3) {
    s <- rowsum(cells$fraction * pal[ch, cells$id], idx)
    m <- matrix(0, nr, nc)
    m[as.integer(rownames(s))] <- s[, 1]
    img[, , ch] <- m
  }
  ink <- function(t) {
    m <- matrix(0, nr, nc)
    if (nrow(t)) {
      s <- rowsum(t$weight, (t$col - 1L) * nr + t$row)
      m[as.integer(rownames(s))] <- s[, 1] * width / area
    }
    pmin(m, 1)
  }
  a_shade <- ink(weights$shade)
  a_lit <- ink(weights$lit)
  for (ch in 1:3) {
    img[, , ch] <- img[, , ch] * (1 - a_shade) + shade_col[ch] * a_shade
    img[, , ch] <- img[, , ch] * (1 - a_lit) + lit_col[ch] * a_lit
  }
  pmin(pmax(img, 0), 1)
}
