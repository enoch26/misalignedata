
# description -------------------------------------------------------------

# This script fits the JU: Joint Uncertainty, VP: Value Plugin, UP: Uncertainty Plugin models. 
# See scenario 4 in Figure 2 in the main text. 

## Incomplete pts obs ----------------------------------------------------------------------
inla.setOption(num.threads = ncore)
### onestage_pts --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_ic_pts_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_ic_pts_sd01.RData"))
  } else {
    z_fml_sd01 <- zsd01 ~ cov_x

    z_lik_sd01 <- bru_obs("Gaussian",
      formula = z_fml_sd01,
      data = z_locs
    )

    # exponential distribution
    cmp_ic <- ~ Intercept(1) + beta_x(
      1,
      mean.linear = 0,
      prec.linear = 1,
      marginal = bru_mapper_marginal(qexp, rate = 1)
    ) +
      cov_x(main = geometry, model = matern) +
      spde(main = geometry, model = maternc)

    fml_ic_pts <- geometry ~ Intercept - beta_x * cov_x + spde

    lik_ic_pts <- bru_obs(
      formula = fml_ic_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_ic_pts_sd01 <- bru(
        components = cmp_ic, z_lik_sd01, lik_ic_pts,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    # iinla: Iteration 6 [max:10]
    #     user   system  elapsed
    # 1294.156   45.851  193.623
    save(fit_ic_pts_sd01, file = "misc/RData/fit_ic_pts_sd01.RData")
  }
}


### onestage_count --------------------------
# sd01
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_ic_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_ic_sd01.RData"))
  } else {
    z_fml_sd01 <- zsd01 ~ cov_x

    z_lik_sd01 <- bru_obs("Gaussian",
      formula = z_fml_sd01,
      data = z_locs
    )

    # exponential distribution
    cmp_ic <- ~ Intercept(1) + beta_x(
      1,
      mean.linear = 0,
      prec.linear = 1,
      marginal = bru_mapper_marginal(qexp, rate = 1)
    ) +
      cov_x(main = geometry, model = matern) +
      spde(main = geometry, model = maternc)

    fml_ic <- pts_count ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - beta_x * cov_x + spde
    )

    lik_ic <- bru_obs(
      formula = fml_ic,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_ic_sd01 <- bru(
        components = cmp_ic, z_lik_sd01, lik_ic,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    # SD0.1
    # iinla: Iteration 4 [max:10]
    #     user   system  elapsed
    # 1967.371   54.760  811.024
    save(fit_ic_sd01, file = "misc/RData/fit_ic_sd01.RData")
  }

}


## Incomplete PolyAgg obs ----------------------------------------------------------------------


### onestage_pts_ap --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_ic_ap_pts.RData"))) {
    load(here::here("misc", "RData", "fit_ic_ap_pts.RData"))
  } else {
    z_fml <- cov ~ ibm_eval(agg_nc_rescale,
      input = list(block = .block, weights = weight),
      state = cov_x
    )

    z_lik <- bru_obs("Gaussian",
      formula = z_fml,
      data = ips_nc,
      response_data = nepal_poly["cov"]
    )

    # exponential distribution
    cmp_ic <- ~ Intercept(1) +
      beta_x(1,
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      cov_x(main = geometry, model = matern) +
      spde(main = geometry, model = maternc)

    fml_ic_pts <- geometry ~ Intercept - beta_x * cov_x + spde

    lik_ic_pts <- bru_obs(
      formula = fml_ic_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_ic_ap_pts <- bru(
        components = cmp_ic, z_lik, lik_ic_pts,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })

    save(fit_ic_ap_pts, file = "misc/RData/fit_ic_ap_pts.RData")
  }
}

### onestage_count_ap --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_ic_ap.RData"))) {
    load(here::here("misc", "RData", "fit_ic_ap.RData"))
  } else {
    z_fml <- cov ~ ibm_eval(agg_nc_rescale,
      input = list(block = .block, weights = weight),
      state = cov_x
    )

    z_lik <- bru_obs("Gaussian",
      formula = z_fml,
      data = ips_nc,
      response_data = nepal_poly["cov"]
    )

    # exponential distribution
    cmp_ic <- ~ Intercept(1) + beta_x(
      1,
      mean.linear = 0,
      prec.linear = 1,
      marginal = bru_mapper_marginal(qexp, rate = 1)
    ) + cov_x(main = geometry, model = matern) +
      spde(main = geometry, model = maternc)

    fml_ic <- pts_count ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - beta_x * cov_x + spde
    )

    lik_ic <- bru_obs(
      formula = fml_ic,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_ic_ap <- bru(
        components = cmp_ic, z_lik, lik_ic,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })

    save(fit_ic_ap, file = "misc/RData/fit_ic_ap.RData")
  }
}

# Two stage ---------------------------------------------------------------

### prep for incomplete PolyAgg ----------------------------------------------------

if (rerun) {
  cmp_ap <- ~ cov_x(main = geometry, model = matern)

  fml_z <- cov ~ ibm_eval(agg_nc_rescale, # make it linear rescale TRUE
    input = list(block = .block, weights = weight),
    state = cov_x
  )

  lik_z <- bru_obs("gaussian",
    formula = fml_z,
    data = ips_nc,
    response_data = nepal_poly["cov"]
  )

  if (file.exists(here::here("misc", "RData", "fit_ap.RData"))) {
    load(here::here("misc", "RData", "fit_ap.RData"))
  } else {
    system.time({
      fit_ap <- bru(
        components = cmp_ap, lik_z,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 5
        )
      )
    })
    save(fit_ap, file = here::here("misc", "RData", "fit_ap.RData"))
  }
  pred_ap <- predict(fit_ap,
    newdata = pts,
    formula = ~ list(
      cov = cov_x
    ),
    n.samples = 100
  )

  pts$icap <- pred_ap$cov$mean
  nepal_rast_icap <- rast(st_rasterize(pts["icap"]))

  fit_apq <- fit_ap$misc$configs$config[[1]]$Q
  fit_apqq <- make_sym(fit_apq)
}

### prep for incomplete PointVal ---------------------------------------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov_sd01.RData"))
  } else {
    cmp_icov <- ~ cov_x(main = geometry, model = matern)

    z_fml_sd01 <- zsd01 ~ cov_x

    z_lik_sd01 <- bru_obs("Gaussian",
      formula = z_fml_sd01,
      data = z_locs
    )

    system.time({
      fit_icov_sd01 <- bru(
        components = cmp_icov, z_lik_sd01,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })
    save(fit_icov_sd01, file = "misc/RData/fit_icov_sd01.RData")
  }


  # prediction step
  # sd01

  pred_icovsd01 <- predict(fit_icov_sd01,
    newdata = pts,
    formula = ~cov_x,
    n.samples = 100,
    seed = seed
  )

  pts$cov_postsd01 <- pred_icovsd01$mean
  nepal_rast_covsd01 <- rast(st_rasterize(pts["cov_postsd01"] %>%
    dplyr::select(cov_postsd01, geometry))) # vect then rast does not work
}



## VP for incomplete pts obs --------------------------------------------------------

### valueplugin_pts --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_icov2_pts_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_pts_sd01.RData"))
  } else {
    cmp_icov2_sd01 <- ~ Intercept(1) +
      cov(nepal_rast_covsd01$cov_postsd01,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2_pts <- geometry ~ Intercept - cov + spde

    lik_icov2_pts <- bru_obs(
      formula = fml_icov2_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov2_pts_sd01 <- bru(
        components = cmp_icov2_sd01, lik_icov2_pts,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })
    save(fit_icov2_pts_sd01, file = here::here("misc", "RData", "fit_icov2_pts_sd01.RData"))
  }
}


### valueplugin_count --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_icov2_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_sd01.RData"))
  } else {
    cmp_icov2_sd01 <- ~ Intercept(1) +
      cov(nepal_rast_covsd01$cov_postsd01,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2 <- pts_count ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - cov + spde
    )

    lik_icov2 <- bru_obs(
      formula = fml_icov2,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icov2_sd01 <- bru(
        components = cmp_icov2_sd01, lik_icov2,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })

    save(fit_icov2_sd01, file = here::here("misc", "RData", "fit_icov2_sd01.RData"))
  }
}


## VP for incomplete PolyAgg obs --------------------------------------------------------

### valueplugin_pts_ap --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov2_ap_pts.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_ap_pts.RData"))
  } else {
    cmp_icov2_ap <- ~ Intercept(1) +
      cov(nepal_rast_icap$icap,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2_pts <- geometry ~ Intercept - cov + spde

    lik_icov2_pts <- bru_obs(
      formula = fml_icov2_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov2_ap_pts <- bru(
        components = cmp_icov2_ap, lik_icov2_pts,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })
    save(fit_icov2_ap_pts, file = here::here("misc", "RData", "fit_icov2_ap_pts.RData"))
  }
}


### valueplugin_count_ap --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov2_ap.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_ap.RData"))
  } else {
    cmp_icov2_ap <- ~ Intercept(1) +
      cov(nepal_rast_icap$icap,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2 <- pts_count ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - cov + spde
    )

    lik_icov2 <- bru_obs(
      formula = fml_icov2,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icov2_ap <- bru(
        components = cmp_icov2_ap, lik_icov2,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })

    save(fit_icov2_ap, file = here::here("misc", "RData", "fit_icov2_ap.RData"))
  }
}


## UP for incomplete pts obs ------------------------------------------------------------


### unplugin_pts --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov3_pts_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov3_pts_sd01.RData"))
  } else {
    # sd01
    fit_qsd01 <- fit_icov_sd01$misc$configs$config[[1]]$Q
    fit_qqsd01 <- make_sym(fit_qsd01)
    cmp_icov3sd01 <- ~ Intercept(1) +
      beta_x(1,
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) + cov_uncertainty(geometry,
        mapper = bru_get_mapper(matern),
        model = "generic0", Cmatrix = fit_qqsd01
      ) +
      spde(main = geometry, model = maternc) +
      cov_pt_est(eval_spatial(nepal_rast_covsd01, geometry), model = "const")

    fml_icov3_pts <- geometry ~ Intercept + spde - beta_x * (cov_pt_est + cov_uncertainty)

    lik_icov3_pts <- bru_obs(
      formula = fml_icov3_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov3_pts_sd01 <- bru(
        components = cmp_icov3sd01, lik_icov3_pts,
        options = list(
          bru_initial = list(beta_x = 0), # exp(0) = 1
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })
    # iinla: Iteration 3 [max:10]
    #    user  system elapsed
    # 363.497  14.858 109.124

    save(fit_icov3_pts_sd01, file = here::here("misc", "RData", "fit_icov3_pts_sd01.RData"))
  }
}


### unplugin_count --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov3_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov3_sd01.RData"))
  } else {
    # sd01
    fit_qsd01 <- fit_icov_sd01$misc$configs$config[[1]]$Q
    fit_qqsd01 <- make_sym(fit_qsd01)
    cmp_icov3sd01 <- ~ Intercept(1) +
      beta_x(1,
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) + cov_uncertainty(geometry,
        mapper = bru_get_mapper(matern),
        model = "generic0", Cmatrix = fit_qqsd01
      ) +
      spde(main = geometry, model = maternc) +
      cov_pt_est(eval_spatial(nepal_rast_covsd01, geometry), model = "const")

    fml_icov3 <- pts_count ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept +
        spde -
        beta_x * (cov_pt_est + cov_uncertainty)
    )

    lik_icov3 <- bru_obs(
      formula = fml_icov3,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc # response data forces allow combine to TRUE
    )

    system.time({
      fit_icov3_sd01 <- bru(
        components = cmp_icov3sd01, lik_icov3,
        options = list(
          bru_initial = list(beta_x = 0),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 10
        )
      )
    })

    save(fit_icov3_sd01, file = here::here("misc", "RData", "fit_icov3_sd01.RData"))
  }
}


## UP for incomplete PolyAgg obs ----------------------------------------------------------------
### unplugin_pts_aggpoly --------------------------

if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icovap_pts.RData"))) {
    load(here::here("misc", "RData", "fit_icovap_pts.RData"))
  } else {
    cmp_icovap <- ~ Intercept(1) +
      beta_x(
        1,
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      cov_uncertainty(geometry,
        mapper = bru_get_mapper(matern),
        model = "generic0", Cmatrix = fit_apqq
      ) +
      spde(main = geometry, model = maternc) +
      cov_pt_est(eval_spatial(nepal_rast_icap, geometry), model = "const")

    fml_icovap_pts <- geometry ~ Intercept + spde - beta_x * (cov_pt_est + cov_uncertainty)

    lik_icovap_pts <- bru_obs(
      formula = fml_icovap_pts,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icovap_pts <- bru(
        components = cmp_icovap, lik_icovap_pts,
        options = list(
          bru_initial = list(beta_x = -1),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    save(fit_icovap_pts, file = here::here("misc", "RData", "fit_icovap_pts.RData"))
  }
}


### unplugin_count_aggpoly --------------------------

if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icovap.RData"))) {
    load(here::here("misc", "RData", "fit_icovap.RData"))
  } else {
    cmp_icovap <- ~ Intercept(1) +
      beta_x(
        1,
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      cov_uncertainty(geometry,
        mapper = bru_get_mapper(matern),
        model = "generic0", Cmatrix = fit_apqq
      ) +
      spde(main = geometry, model = maternc) +
      cov_pt_est(eval_spatial(nepal_rast_icap, geometry), model = "const")

    fml_icovap <- pts_count ~ ibm_eval(agg_nc,
      input = list(
        block = .block,
        weights = weight
      ),
      state = Intercept +
        spde -
        beta_x * (cov_pt_est + cov_uncertainty)
    )

    lik_icovap <- bru_obs(
      formula = fml_icovap,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icovap <- bru(
        components = cmp_icovap, lik_icovap,
        options = list(
          bru_initial = list(beta_x = -1),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    save(fit_icovap, file = here::here("misc", "RData", "fit_icovap.RData"))
  }
}
