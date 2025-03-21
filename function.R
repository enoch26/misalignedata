
# description -------------------------------------------------------------

# This script contains the functions used in the simulation study.
# The functions are self-explanatory and are used in the main script.

# Covariate deterministic function ----------------------------------------


cov_fun <- function(x, y) {
  (x^2 - y^2) * exp(-.5 * x^2 - .5 * y^2)
}


# Nonlinear transformation function ---------------------------------------


exp_nl <- function(x) (exp(a * x) - c_) / b



# invert scale ------------------------------------------------------------

inv_scale <- function(x, xmin, xmax, xmin_new, xmax_new) {
  x_new <- (x - xmin) / (xmax - xmin) * (xmax_new - xmin_new) + xmin_new
  return(x_new)
}


# hexagon lattice  --------------------------------------------------------


# bnd: boundary
# x_bin: # of bin in x axis
# edge_len_n: # of edge length of mesh from the boundary to create hexagon mesh using x_bin
# S2 not supported
hexagon_lattice <- function(bnd = bnd,
                            x_bin = 250,
                            edge_len_n = 1) {
  stopifnot(x_bin / 2 > edge_len_n)
  crs <- fm_crs(bnd)
  fm_crs(bnd) <- NA
  # two separate grid and combine
  edge_len <- as.numeric((fm_bbox(bnd)[[1]][2] - fm_bbox(bnd)[[1]][1])) / (x_bin + 2 * edge_len_n) # two ends
  # sf_buffer to work on negative buffer to stay a distance from the boundary
  # Turn off S2 to avoid zig zag
  # suppressMessages(sf::sf_use_s2(FALSE))
  bnd_inner <- st_buffer(bnd, dist = -edge_len_n * edge_len) # st_buffer for edge_len x1
  y_diff <- fm_bbox(bnd_inner)[[2]][2] - fm_bbox(bnd_inner)[[2]][1]
  x_diff <- fm_bbox(bnd_inner)[[1]][2] - fm_bbox(bnd_inner)[[1]][1]
  y_bin <- as.integer(y_diff / (sqrt(3) / 2 * edge_len))
  # TODO rep n, n-1, length
  h <- (sqrt(3) / 2 * edge_len) # height
  x_adj <- .5 * (x_diff - x_bin * edge_len)
  y_adj <- .5 * (y_diff - y_bin * h)
  # x
  x_1_ <- seq(
    fm_bbox(bnd_inner)[[1]][1] + x_adj,
    fm_bbox(bnd_inner)[[1]][2] - x_adj, edge_len
  )
  x_2_ <- seq(
    (fm_bbox(bnd_inner)[[1]][1] + x_adj + .5 * edge_len),
    (fm_bbox(bnd_inner)[[1]][2] - x_adj - .5 * edge_len),
    edge_len
  )
  y_1_ <- seq(
    fm_bbox(bnd_inner)[[2]][1] + y_adj,
    fm_bbox(bnd_inner)[[2]][2] - y_adj,
    by = 2 * h
  )
  y_2_ <- seq(
    fm_bbox(bnd_inner)[[2]][1] + y_adj + h,
    fm_bbox(bnd_inner)[[2]][2] - y_adj + h, 2 * h
  )

  x_1 <- rep(x_1_, times = length(y_1_))
  x_2 <- rep(x_2_, times = length(y_2_))
  y_1 <- rep(y_1_, each = length(x_1_))
  y_2 <- rep(y_2_, each = length(x_2_))

  mesh_df <- data.frame(x = c(x_1, x_2), y = c(y_1, y_2))
  # turn the mesh nodes into lattice sf
  lattice_sf <- st_as_sf(mesh_df, coords = c("x", "y"), crs = st_crs(bnd))
  lattice_sfc <- lattice_sf %>% st_as_sfc()
  pts_inside <- lengths(st_intersects(lattice_sfc, bnd_inner)) != 0
  pts_lattice_sfc <- lattice_sfc[pts_inside]
  fm_crs(pts_lattice_sfc) <- fm_crs(bnd_inner) <- crs
  return(list(
    lattice = pts_lattice_sfc,
    edge_len = edge_len,
    bnd_inner = bnd_inner
  ))
}

# score -------------------------------------------------------------------
# A function to compute the SE and DS score for bru object(s)
# obs: the response name in the newdata to compare with the prediction
# formula(e): the formula to predict the response;
#             if NULL, the formula in the model is used;
# seed: seed for reproducibility in the the global seed and prediction;
#       seed restored after the computation
# plot: logical, whether to plot the SE and DS score
# inla.link: the INLA link function suffix (can be found in INLA::link), eg "log"
# to inverse the prediction
# set seed issue: https://github.com/inlabru-org/inlabru/issues/88
# now only work for Gaussian (and Poisson, because Poisson we can turn into Gaussian)
dse_score <- function(object,
                      newdata = NULL,
                      formula = NULL,
                      n.samples = 100,
                      obs = "lambda",
                      inla.link = NULL,
                      seed = NULL,
                      plot = FALSE,
                      ...) {
  if (length(seed) > 2) {
    stop("At most two seeds should be provided")
  }

  if (length(seed) >= 1) {
    set.seed(seed[1])
  }

  if (length(seed) >= 2) {
    seed2 <- seed[2]
  } else {
    seed2 <- as.integer(runif(1) * .Machine$integer.max) # this is how INLA set seed
  }

  if (is.null(inla.link)) {
    inla.link <- object$misc$linkfunctions[["names"]] # TODO what if multiple likelihood
    if (length(inla.link) > 1) {
      warning(paste0(
        "Multiple link functions (",
        paste0(inla.link, collapse = ", "),
        ") detected. Please provide an inla.link argument"
      ))
    }
  }

  # Use the above seed to reproduce the same result and into the predict call
  if (is.null(formula)) {
    formula <- object$bru_info$lhoods[[1]]$formula
    formula <- as.formula(paste0(
      "~ inla.link.", inla.link,
      "(", labels(terms(formula)), ", inverse = TRUE)"
    ))
  }
  pred <- predict(object,
    newdata = newdata, formula = formula,
    n.samples = n.samples,
    seed = seed2,
    ...
  )
  post_E <- get("mean", pred)
  post_Var <- (get("sd", pred))^2
  # TODO this only applies to Poisson lambda/ Gaussian mean, count needa add post_E to the variance, other cases need to do case by case.
  pred$SE_score <- SE_score <- (get(obs, newdata) - post_E)^2
  pred$DS_score <- DS_score <- SE_score / post_Var + log(post_Var)
  SE_score_mean <- mean(SE_score, na.rm = TRUE)
  DS_score_mean <- mean(DS_score, na.rm = TRUE)
  cat(
    "Mean SE score: ", SE_score_mean, "\n",
    "Mean DS score: ", DS_score_mean, "\n"
  )
  # An example of how to plot the SE and DS score
  if (plot) {
    p1 <- ggplot() +
      gg(data = pred["SE_score"], geom = "tile") +
      ggtitle("Squared Error Score")
    p2 <- ggplot() +
      gg(data = pred["DS_score"], geom = "tile") +
      ggtitle("Dawid-Sebastiani Score")
    print(p1)
    print(p2)
  }

  return(list(
    pred = pred,
    mean_score = tibble::tibble(
      SE_score_mean = SE_score_mean,
      DS_score_mean = DS_score_mean
    ),
    score = tibble::tibble(SE_score = SE_score, DS_score = DS_score)
  ))
}


# make Q symmetric --------------------------------------------------------
# thanks to Stephen Villejo for sharing

# A function to make a symmetric matrix for the Q matrix from the INLA model
make_sym <- function(Q) {
  d <- diag(Q)
  Q <- (Q + t(Q))
  diag(Q) <- d
  return(Q)
}


# Subheading for ggplot ---------------------------------------------------

gg_subhd <- function(string, x = 0, y = 1,...){
  subhd <- paste0("<b>", sprintf("%s", string), "</b>")
  ggplot(data = tibble::tibble(x = x, y = y, label = subhd)) +
    aes(x = x, y = y, label = label) +
    ggtext::geom_textbox(
      box.color = NA,
      fill = NA,
      width = unit(10, "cm"),
      ...
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}
