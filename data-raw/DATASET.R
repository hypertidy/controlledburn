## code to prepare `DATASET` dataset goes here

info <- vapour::vapour_raster_info("https://i.imgflip.com/23jjd0.jpg")
d <- vapour::vapour_warp_raster_hex("https://i.imgflip.com/23jjd0.jpg",
                                    extent = info$extent, dimension = info$dimension,
                                    transformation_options = "SRC_METHOD=NO_GEOTRANSFORM")

hex <- matrix(d, info$dimension[2], byrow = T)[1700:info$dimension[2], 665:info$dimension[1]]
ximage::ximage(hex, asp = 1)

# library(hexSticker)
# s <- sticker(~ximage::ximage(hex, asp = 1),
#              package="laserize", #p_size=20, s_x=.8, s_y=.6, s_width=1.4, s_height=1.2,
#              filename="man/figures/logo.png")
png::writePNG(scales::rescale(array(t(col2rgb(hex)), c(dim(hex), 3))), "logo.png")
abline(v = 1120, h = 8)
usethis::use_data(DATASET, overwrite = TRUE)
