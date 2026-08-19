## code to prepare `RUSANT` dataset goes here
## sds::CGAZ()
cgaz <- "/vsicurl/https://github.com/mdsumner/geoboundaries/releases/download/latest/geoBoundariesCGAZ_ADM0.parquet"
#sds::CGAZ_sql(c("Antarctica", "Russia"))
sql <- "SELECT shapeGroup FROM geoBoundariesCGAZ_ADM0 WHERE shapeGroup IN ('ATA','RUS')"
wkb <- wk::wkb(vapour::vapour_read_geometry(cgaz, sql = sql))
rusant <- wk::as_wkb(geos::geos_simplify_preserve_topology(wkb, .01))

usethis::use_data(rusant, internal = TRUE)
