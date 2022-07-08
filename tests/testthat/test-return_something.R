library(sfheaders)
library(laserize)
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
system.time({
r <- laserize:::laserize(pols, extent = c(range(pdata$x), range(pdata$y)),
                               dimension = c(50000, 40000))
})


test_that("multiplication works", {
  expect_equal(r, 1L)
})
