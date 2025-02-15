
# description -------------------------------------------------------------

# This script loads the data and creates the covariate raster. It also creates 
# the mesh and simulates from a Matern GRF. The script also evaluates the 
# covariate at the points and integration points. The script also evaluates the 
# linear predictor and nonlinear transformation at the points and mesh vertices. 
# The script also samples from the LGCP model and prepares the data for scoring.

# pts: pointwise simulated data, evenly distributed sf points across space
# ips_nc: integration points 
# bnd: boundary
# bnd_out: outer boundary
# bnd_buff: smoothed boundary

# read data ---------------------------------------------------------------
source(here("read_data.R"))

# create covariate raster -------------------------------------------------
source(here("covariate_rast.R"))

# Mesh --------------------------------------------------------------------
source(here("mesh.R"))

# Simulate from Matern GRF --------------------------------------------------------------------
# regular pts
pts <- st_sf(
  geometry = st_sample(bnd_buff,
    type = "regular",
    size = 250 * 250, 
    crs = crs_nepal
  )
)

if (to_plot) {
  ggplot() +
    geom_spatraster(data = nepal_rast, aes(fill = z)) +
    gg(data = pts)
}


# matern covariance parameters
rho <- 50
sigma <- 1 / 2
samp <- fm_matern_sample(mesh_fm, rho = rho, sigma = sigma)
# [, 1]

# evaluate pts with matern field
pts$u <- fm_evaluate(mesh_fm, loc = pts, field = samp)

# extract covariate field with pts in a finer resolution, ie repeated sampling
pts_vect <- vect(st_as_sfc(pts))

# true covariate pointwise
pts$cov_pts <- cov_fun(
  inv_scale(st_coordinates(pts)[, 1],
    xmin = bnd_out_bbox$xmin,
    xmax = bnd_out_bbox$xmax,
    xmin_new = -4,
    xmax_new = 4
  ),
  inv_scale(st_coordinates(pts)[, 2],
    xmin = bnd_out_bbox$ymin,
    xmax = bnd_out_bbox$ymax,
    xmin_new = -2,
    xmax_new = 2
  )
)
# raster
pts$cov <- unlist(terra::extract(nepal_rast$z, pts_vect,
  ID = FALSE,
  na.rm = TRUE
))
# aggregated raster
pts$cov_agg <- unlist(terra::extract(nepal_rast_agg$z, pts_vect,
  ID = FALSE,
  na.rm = TRUE
))


# evaluate covariate at pts
system.time({
  pts$cov_poly <- eval_spatial(
    data = nepal_poly,
    where = pts, layer = "cov"
  )
})

# integration points ------------------------------------------------------


# this create the .block for each pixel
system.time({
  ips_nc <- fm_int(mesh_fm, samplers = nepal_nc)
}) # this takes a while
# user  system elapsed
# 129.140   0.553 128.728
block_nc <- fm_block(block = ips_nc$.block, weights = ips_nc$weight)

ips_nc$cov_pts <- cov_fun(
  inv_scale(st_coordinates(ips_nc)[, 1],
    xmin = bnd_out_bbox$xmin,
    xmax = bnd_out_bbox$xmax,
    xmin_new = -4,
    xmax_new = 4
  ),
  inv_scale(st_coordinates(ips_nc)[, 2],
    xmin = bnd_out_bbox$ymin,
    xmax = bnd_out_bbox$ymax,
    xmin_new = -2,
    xmax_new = 2
  )
)

ips_nc$cov <- extract(nepal_rast, ips_nc,
  fun = "mean",
  ID = FALSE, na.rm = TRUE
)

ips_nc$cov_agg <- extract(nepal_rast_agg, ips_nc,
  fun = "mean", ID = FALSE, na.rm = TRUE
)

ips_nc$covly <- eval_spatial(nepal_poly, ips_nc, layer = "cov")

# Formula -----------------------------------------------------------------
# linear predictor parameters
intercept <- -7
b_x <- -6 # coeff of covariate
pts$loglambda <- intercept + b_x * pts$cov_pts + pts$u
pts$lambda <- exp(intercept + b_x * pts$cov_pts + pts$u)

# nonlinear transformation
pts$loglambda_nl <- intercept + b_x * exp_nl(pts$cov_pts) + pts$u
pts$lambda_nl <- exp(intercept + b_x * exp_nl(pts$cov_pts) + pts$u)

## nepal_nc state ----------------------------------------------------------
pts$ID <- st_join(pts, nepal_nc)$ID

# Evaluate mesh vertex
# field: Basis function weights, one per mesh basis function, describing the function to be evaluated at the projection locations
mesh_fm_vtx <- fm_vertices(mesh_fm, format = "sf")

if (to_plot) {
  ggplot() +
    geom_spatraster(data = nepal_rast, aes(fill = z)) +
    gg(mesh_fm) +
    gg(mesh_fm_vtx, size = .5)
}

mesh_fm_vtx$cov_pts <- cov_fun(
  inv_scale(st_coordinates(mesh_fm_vtx)[, 1],
    xmin = bnd_out_bbox$xmin,
    xmax = bnd_out_bbox$xmax,
    xmin_new = -4,
    xmax_new = 4
  ),
  inv_scale(st_coordinates(mesh_fm_vtx)[, 2],
    xmin = bnd_out_bbox$ymin,
    xmax = bnd_out_bbox$ymax,
    xmin_new = -2,
    xmax_new = 2
  )
)


if (to_plot) {
  ggplot() +
    gg(mesh_fm) +
    gg(mesh_fm_vtx, aes(col = cov), size = 1)
}

mesh_fm_vtx$u <- fm_evaluate(mesh = mesh_fm, loc = mesh_fm_vtx, field = samp)

if (to_plot) {
  ggplot() +
    gg(mesh_fm) +
    gg(mesh_fm_vtx, aes(col = spde), size = 1)
}

mesh_fm_vtx$loglambda <- as.vector(with(mesh_fm_vtx, intercept + b_x * cov_pts + u))

mesh_fm_vtx$loglambda_nl <-
  as.vector(with(mesh_fm_vtx, intercept + b_x * exp_nl(cov_pts) / b + u))

# sample pts -------------------------------------------------------------

pts_samp <- sample.lgcp(
  mesh = mesh_fm,
  loglambda = get("loglambda", mesh_fm_vtx),
  samplers = bnd
)

pts_samp_sf <- st_as_sf(pts_samp)

nepal_nc$pts_count <- lengths(st_intersects(nepal_nc, pts_samp_sf))

# nonlinear transformation
pts_samp_nl <- sample.lgcp(
  mesh = mesh_fm,
  loglambda = get("loglambda_nl", mesh_fm_vtx),
  samplers = bnd
)
pts_samp_nl_sf <- st_as_sf(pts_samp_nl)
nepal_nc$pts_count_nl <- lengths(st_intersects(nepal_nc, pts_samp_nl_sf))

# pts_int for scoring -----------------------------------------------------
pts_inside <- lengths(st_intersects(pts, bnd)) != 0
pts_int <- pts[pts_inside, , drop = FALSE]
st_geometry(pts_int) <- "geometry"



# prepare for incomplete data --------------------------------------------
z_locs <- st_sf(
  geometry = st_sample(st_cast(nepal_nc, "MULTIPOLYGON"),
    size = rep(1, nrow(st_cast(nepal_nc, "MULTIPOLYGON")))
  )
)

# sample SPDE field
z_locs$spde <- fm_evaluate(mesh_fm, loc = z_locs, field = samp)
z_locs_vect <- vect(st_as_sfc(z_locs))
# sample cov field
z_locs$cov_x <- unlist(terra::extract(nepal_rast$z, z_locs_vect,
  ID = FALSE,
  xy = FALSE,
  na.rm = TRUE
))

z_locs$cov_xnl <- unlist(terra::extract(nepal_rast$znl, z_locs_vect,
  ID = FALSE,
  xy = FALSE,
  na.rm = TRUE
))
z_locs$z <- z_locs$cov_x + rnorm(100, sd = .05)
z_locs$zsd01 <- z_locs$cov_x + rnorm(100, sd = .1)

rm(
  block_nc,
  cov_val,
  boundary,
  bnd1,
  bnd2,
  xbuff,
  xrange,
  ybuff,
  yrange,
  pts_inside,
  pts_samp,
  pts_samp_nl,
  df,
  idx
)

# plot
to_plot <- FALSE
tw <- 15.55528 # page width in cm
# tw <- 15.55528 / 2.54 # page width in inches

