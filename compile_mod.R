
# description -------------------------------------------------------------

# This script compiles everything in the simulation study.

# Prep -----------------------------------------------------------------------

## global_option --------------------------
# rerun everything
rerun <- TRUE
# rerun <- FALSE
# plot
to_plot <- FALSE
tw <- 15.55528 # page width in cm
# tw <- 15.55528 / 2.54 # page width in inches


## library --------------------------
if (rerun) {
  # library
  # core
  pkgs <- c(
    "INLA", "inlabru", "fmesher", "sf", # core
    "ggplot2", "here", "dplyr", # general
    "terra", "stars", # raster
    "patchwork", "tidyterra", # plot
    "future"
  )

  for (i in 1:length(pkgs)) {
    suppressMessages(library(pkgs[i], character.only = TRUE))
  }
  
  # set seed
  seed <- 2626
  # seed <- 3738 # another realisation
  set.seed(seed)

  # covariate field
  factor <- 10 # factor of aggregation
  resol <- 1000 

  # mesh
  x_bin <- 200
  edge_len_n <- 2

  # nonlinear transformation function
  a <- 3
  b <- 9
  c_ <- 3

  # source function for hexagon mesh
  source(here("function.R"))
}


## data --------------------------
if (rerun) {
  if (file.exists(here("misc", "RData", "load_data.RData"))) {
    load(here("misc", "RData", "load_data.RData"))
    
  } else {
    system.time({
      source(here("load_data.R"))
    })

    rm(nepal_rast, nepal_rast_agg)

    save.image("misc/RData/load_data.RData")
  }

  if (!exists("nepal_rast")) {
    source(here("covariate_rast.R"))
  }
  if (!exists("mesh_fm")) {
    source(here("mesh.R"))
  }
}


## matern --------------------------
if (rerun) {
  matern <- inla.spde2.pcmatern(mesh_fm,
    prior.range = c(4, 0.1),
    prior.sigma = c(2, 0.1)
  )
  maternc <- inla.spde2.pcmatern(mesh_fmc,
    prior.range = c(8, 0.1),
    prior.sigma = c(2, 0.1)
  )
}


### aggregate_mapper --------------------------
if (rerun) {
  agg_nc <- bru_mapper_logsumexp(rescale = FALSE, n_block = nrow(nepal_nc))
  agg_nc_rescale <- bru_mapper_aggregate(rescale = TRUE, n_block = nrow(nepal_nc))
}


# Modelling --------------------------------------------------------------------
ncore <- 10
inla.setOption(num.threads = ncore)
## Joint -------------------------------------------------------------

source("joint_mod.R")

## Incomplete covariate field for Linear Y --------------------------------------------------

source("juvpup.R")

## Nonlinear misspecification for NonLinear Y --------------------------------------------------

source("juvpup_nl.R")
