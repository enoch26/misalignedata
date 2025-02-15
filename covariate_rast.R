# Covariate Raster layer------------------------------------------------------------------
# for terra
# Create a SpatRaster from scratch

resol <- resol # controls the horizontal resolution

# create raster
bnd_out_bbox <- st_bbox(bnd_out)
xrange <- as.numeric(bnd_out_bbox$xmax - bnd_out_bbox$xmin)
yrange <- as.numeric(bnd_out_bbox$ymax - bnd_out_bbox$ymin)
xbuff <- 0
ybuff <- 0

df <- expand.grid(
  x = seq(bnd_out_bbox$xmin - xbuff,
    bnd_out_bbox$xmax + xbuff,
    length.out = resol
  ),
  y = seq(bnd_out_bbox$ymin - ybuff,
    bnd_out_bbox$ymax + ybuff,
    length.out = round(resol * yrange / xrange)
  )
)

nepal_rast <- rast(df,
  type = "xyz",
  crs = crs_nepal$input
)


cov_val <- expand.grid(
  x = seq(-4, 4, length.out = dim(nepal_rast)[2]),
  y = seq(-2, 2, length.out = dim(nepal_rast)[1])
)


cov_val$z <- with(cov_val, cov_fun(x, y))

# nonlinear transform
cov_val$znl <- exp_nl(cov_val$z)

nepal_rast$z <- cov_val$z
nepal_rast$znl <- cov_val$znl

bnd_buff <- fm_nonconvex_hull(bnd, 5) # 20240613

if (to_plot) {
  ggplot() +
    geom_spatraster(data = nepal_rast, aes(fill = z))

  ggplot() +
    geom_spatraster(data = nepal_rast, aes(fill = z)) +
    gg(data = st_as_sf(bnd_out), col = "blue") +
    gg(data = st_as_sf(bnd_buff), col = "red") +
    gg(data = bnd) +
    gg(data = nepal_nc) +
    coord_sf()

  ggplot() +
    geom_spatraster(data = nepal_rast, aes(fill = z)) +
    gg(data = nepal_nc)
}
# we need a bit bigger for integration point ips
nepal_rast <- crop(nepal_rast, vect(bnd_out)[1, ], mask = TRUE) # make sure ips does not have NA

if (to_plot) {
  # terra has smaller pdf size
  nepal_rast1 <- crop(nepal_rast, vect(bnd)[1, ], mask = TRUE)
  p_covrast <- ggplot() +
    gg(data = nepal_rast1$z, na.rm = TRUE) +
    # gg(bnd) +
    # scale_fill_continuous(na.value = "transparent") +
    scale_fill_viridis_c(na.value = "transparent") +
    coord_sf()
  p_covrast
}

nepal_rast_agg <- terra::aggregate(nepal_rast,
  fact = factor, fun = mean,
  na.rm = TRUE
)

if (to_plot) {
  nepal_rast1 <- crop(nepal_rast, vect(bnd)[1, ], mask = TRUE)
  nepal_rast_agg1 <- crop(nepal_rast_agg, vect(bnd)[1, ], mask = TRUE)
  prast <- ggplot() +
    geom_spatraster(data = nepal_rast1, aes(fill = z)) +
    coord_sf() +
    scale_fill_viridis_c(na.value = "transparent") +
    theme(legend.position = "none") +
    ggtitle("True")
  prast_agg <- ggplot() +
    geom_spatraster(data = nepal_rast_agg1, aes(fill = z)) +
    coord_sf() +
    scale_fill_viridis_c(na.value = "transparent") +
    theme(legend.position = "none") +
    ggtitle("Mean Agg")
  (prast + prast_agg) +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "bottom") # it does not go well with spatraster, can turn off other legends
}

# Regional Aggregate  ------------------------------------------------------
nepal_poly <- cbind(nepal_nc, extract(nepal_rast, nepal_nc,
  fun = "mean",
  na.rm = TRUE
)) %>% rename(cov = z)
if (to_plot) {
  ppoly <- ggplot() +
    geom_sf(data = nepal_poly, aes(fill = cov)) +
    scale_fill_viridis_c()
  ppoly
}
