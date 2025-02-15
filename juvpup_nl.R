
# description -------------------------------------------------------------

# This script fits the JU: Joint Uncertainty, VP: Value Plugin, UP: Uncertainty Plugin 
# models under nonlinear misspecification. 
# See scenario 5 in Figure 2 in the main text. 

## One stage NL ----------------------------------------------------------------------
inla.setOption(num.threads = ncore)
### onestage_pts_nl --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_ic_pts_sd01_nl.RData"))) {
    load(here::here("misc", "RData", "fit_ic_pts_sd01_nl.RData"))
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

    lik_ic_pts_nl <- bru_obs(
      formula = fml_ic_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_ic_pts_sd01_nl <- bru(
        components = cmp_ic, z_lik_sd01, lik_ic_pts_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    save(fit_ic_pts_sd01_nl, file = "misc/RData/fit_ic_pts_sd01_nl.RData")
  }
}


### onestage_count_nl --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_ic_sd01_nl.RData"))) {
    load(here::here("misc", "RData", "fit_ic_sd01_nl.RData"))
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

    fml_ic_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - beta_x * cov_x + spde
    )

    lik_ic_nl <- bru_obs(
      formula = fml_ic_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_ic_sd01_nl <- bru(
        components = cmp_ic, z_lik_sd01, lik_ic_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    save(fit_ic_sd01_nl, file = "misc/RData/fit_ic_sd01_nl.RData")
  }
}






## Incomplete PolyAgg obs ----------------------------------------------------------------------
### onestage_pts_ap_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_ic_ap_pts_nl.RData"))) {
    load(here::here("misc", "RData", "fit_ic_ap_pts_nl.RData"))
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

    lik_ic_pts_nl <- bru_obs(
      formula = fml_ic_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_ic_ap_pts_nl <- bru(
        components = cmp_ic, z_lik, lik_ic_pts_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })

    save(fit_ic_ap_pts_nl, file = "misc/RData/fit_ic_ap_pts_nl.RData")
  }
}
#
#
### onestage_count_ap_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_ic_ap_nl.RData"))) {
    load(here::here("misc", "RData", "fit_ic_ap_nl.RData"))
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

    fml_ic_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept - beta_x * cov_x + spde
    )

    lik_ic_nl <- bru_obs(
      formula = fml_ic_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_ic_ap_nl <- bru(
        components = cmp_ic, z_lik, lik_ic_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })

    save(fit_ic_ap_nl, file = "misc/RData/fit_ic_ap_nl.RData")
  }
}

# Two stage NL---------------------------------------------------------------
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

### prep for incomplete pts obs ---------------------------------------------------------
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
    dplyr::select(cov_postsd01, geometry)))
}
# @



## VP for incomplete pts obs --------------------------------------------------------
### valueplugin_pts_nl --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_icov2_pts_nl_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_pts_nl_sd01.RData"))
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

    lik_icov2_pts_nl <- bru_obs(
      formula = fml_icov2_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov2_pts_nl_sd01 <- bru(
        components = cmp_icov2_sd01, lik_icov2_pts_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    save(fit_icov2_pts_nl_sd01, file = here::here("misc", "RData", "fit_icov2_pts_nl_sd01.RData"))
  }
}


### valueplugin_count_nl --------------------------
if (rerun) {
  # sd01
  if (file.exists(here::here("misc", "RData", "fit_icov2_nl_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_nl_sd01.RData"))
  } else {
    cmp_icov2_sd01 <- ~ Intercept(1) +
      cov(nepal_rast_covsd01$cov_postsd01,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(
        block = .block,
        weights = weight
      ),
      state = Intercept - cov + spde
    )

    lik_icov2_nl <- bru_obs(
      formula = fml_icov2_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icov2_nl_sd01 <- bru(
        components = cmp_icov2_sd01, lik_icov2_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    #     iinla: Iteration 4 [max:10]
    #    user  system elapsed
    # 128.110   7.680  82.472

    save(fit_icov2_nl_sd01, file = here::here("misc", "RData", "fit_icov2_nl_sd01.RData"))
  }
}




## VP for incomplete PolyAgg obs --------------------------------------------------------

### valueplugin_pts_ap_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov2_ap_pts_nl.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_ap_pts_nl.RData"))
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

    lik_icov2_pts_nl <- bru_obs(
      formula = fml_icov2_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov2_ap_pts_nl <- bru(
        components = cmp_icov2_ap, lik_icov2_pts_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    save(fit_icov2_ap_pts_nl, file = here::here("misc", "RData", "fit_icov2_ap_pts_nl.RData"))
  }
}


### valueplugin_count_ap_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov2_nl_ap.RData"))) {
    load(here::here("misc", "RData", "fit_icov2_nl_ap.RData"))
  } else {
    cmp_icov2_ap <- ~ Intercept(1) +
      cov(nepal_rast_icap$icap,
        model = "linear",
        mean.linear = 0,
        prec.linear = 1,
        marginal = bru_mapper_marginal(qexp, rate = 1)
      ) +
      spde(main = geometry, model = maternc)

    fml_icov2_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(
        block = .block,
        weights = weight
      ),
      state = Intercept - cov + spde
    )

    lik_icov2_nl <- bru_obs(
      formula = fml_icov2_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icov2_nl_ap <- bru(
        components = cmp_icov2_ap, lik_icov2_nl,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    #     iinla: Iteration 4 [max:10]
    #    user  system elapsed
    # 128.110   7.680  82.472

    save(fit_icov2_nl_ap, file = here::here("misc", "RData", "fit_icov2_nl_ap.RData"))
  }
}



## UP for incomplete pts obs ------------------------------------------------------------




### unplugin_pts_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov3_pts_nl_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov3_pts_nl_sd01.RData"))
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

    fml_icov3_pts <- geometry ~ Intercept + spde -
      beta_x * (cov_pt_est + cov_uncertainty)

    lik_icov3_pts_nl <- bru_obs(
      formula = fml_icov3_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icov3_pts_nl_sd01 <- bru(
        components = cmp_icov3sd01, lik_icov3_pts_nl,
        options = list(
          bru_initial = list(beta_x = 0), # exp(0) = 1
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })
    # iinla: Iteration 4 [max:10]
    #    user  system elapsed
    # 250.809   9.988  73.011
    save(fit_icov3_pts_nl_sd01, file = here::here("misc", "RData", "fit_icov3_pts_nl_sd01.RData"))
  }
}

### unplugin_count_nl --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icov3_nl_sd01.RData"))) {
    load(here::here("misc", "RData", "fit_icov3_nl_sd01.RData"))
  } else {
    # sd01
    fit_qsd01 <- fit_icov_sd01$misc$configs$config[[1]]$Q
    fit_qqsd01 <- make_sym(fit_qsd01)
    cmp_icov3_sd01 <- ~ Intercept(1) +
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

    fml_icov3_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(block = .block, weights = weight),
      state = Intercept +
        spde -
        beta_x * (cov_pt_est + cov_uncertainty)
    )

    lik_icov3_nl <- bru_obs(
      formula = fml_icov3_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icov3_nl_sd01 <- bru(
        components = cmp_icov3_sd01, lik_icov3_nl,
        options = list(
          bru_initial = list(beta_x = 0),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 30
        )
      )
    })

    save(fit_icov3_nl_sd01, file = here::here("misc", "RData", "fit_icov3_nl_sd01.RData"))
  }
}



## UP for incomplete PolyAgg obs ----------------------------------------------------------------
### unplugin_pts_aggpoly_nl --------------------------

if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icovap_pts_nl.RData"))) {
    load(here::here("misc", "RData", "fit_icovap_pts_nl.RData"))
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

    lik_icovap_pts_nl <- bru_obs(
      formula = fml_icovap_pts,
      family = "cp",
      data = pts_samp_nl_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )

    system.time({
      fit_icovap_pts_nl <- bru(
        components = cmp_icovap, lik_icovap_pts_nl,
        options = list(
          bru_initial = list(beta_x = -1),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    save(fit_icovap_pts_nl, file = here::here("misc", "RData", "fit_icovap_pts_nl.RData"))
  }
}


### unplugin_count_aggpoly_nl --------------------------

if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_icovap_nl.RData"))) {
    load(here::here("misc", "RData", "fit_icovap_nl.RData"))
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

    fml_icovap_nl <- pts_count_nl ~ ibm_eval(agg_nc,
      input = list(
        block = .block,
        weights = weight
      ),
      state = Intercept +
        spde -
        beta_x * (cov_pt_est + cov_uncertainty)
    )

    lik_icovap_nl <- bru_obs(
      formula = fml_icovap_nl,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )

    system.time({
      fit_icovap_nl <- bru(
        components = cmp_icovap, lik_icovap_nl,
        options = list(
          bru_initial = list(beta_x = -1),
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    save(fit_icovap_nl, file = here::here("misc", "RData", "fit_icovap_nl.RData"))
  }
}
