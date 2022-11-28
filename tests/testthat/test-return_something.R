library(sfheaders)
library(minorburn)
# p1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60), c(-180,-20))
# hole <- rbind(c(-150,-20), c(-100,-10), c(-110,20), c(-150,-20))
# p1 <- list(p1, hole)
# p2 <- list(rbind(c(-10,0), c(140,60), c(160,0), c(140,-55), c(-10,0)))
# p3 <- list(rbind(c(-125,0), c(0,60), c(40,5), c(15,-45), c(-125,0)))
# pols <- st_sf(value = rep(1,3),
#               geometry = st_sfc(lapply(list(p1, p2, p3), st_polygon)))
pdata <-
  structure(list(sfg_id = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
                 polygon_id = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
                 linestring_id = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                 x = c(-180, -140, 10, -140, -180, -150, -100, -110, -150, -10, 140, 160, 140, -10, -125, 0, 40, 15, -125),
                 y = c(-20, 55, 0, -60, -20, -20, -10, 20, -20, 0, 60, 0, -55, 0, 0, 60, 5, -45, 0)),
            class = "data.frame", row.names = c(NA, 19L))

pols <- sfheaders::sf_polygon(pdata, x= "x", y = "y", polygon_id = "polygon_id", linestring_id = "linestring_id")
ex <- c(range(pdata$x), range(pdata$y))
dm <- c(50, 40)
r <- laserize(pols, extent = ex,
                               dimension =dm)

## we are triplets xstart, xend, yline (0-based atm)
index <- do.call(rbind, r)

## is this the old problem in fasterize, it's assuming centres?
m <- matrix(FALSE, 50, 40)
m[index[,c(1, 3)] + 1]<- TRUE
m[index[,c(2, 3)] + 1]<- TRUE
ximage::ximage(t(m), extent = ex)
# library(basf)
# plot(pols, add = TRUE)
# rast(ext(ex), nrows = dm[2], ncols = dm[1], vals = m)
#
# library(basf)
# plot(pols, add = T)
test_that("multiplication works", {
  expect_true(is.integer(r[[1]]))
})
