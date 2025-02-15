
# description -------------------------------------------------------------

# This script fits the Observation Plugin models. 
# See scenarios 1-3 in Figure 2 in the main text. 

## Pts ---------------------------------------------------------------------
inla.setOption(num.threads = ncore)

### pts_rast --------------------------
system.time({
  if (rerun) {
    if (file.exists(here("misc", "RData", "fit_pts.RData"))) {
      load(here("misc", "RData", "fit_pts.RData"))
    } else {
      cmp_pts <- ~ Intercept(1) + cov(nepal_rast$z, model = "linear") +
        spde(main = geometry, model = matern)
      
      fml_pts <- geometry ~ Intercept + cov + spde
      
      lik_pts <- bru_obs(
        formula = fml_pts,
        family = "cp",
        data = pts_samp_sf,
        domain = list(geometry = mesh_fm),
        samplers = bnd
      )
      
      system.time({
        fit_pts <- bru(
          components = cmp_pts, lik_pts,
          options = list(
            control.inla = list(int.strategy = "eb"),
            bru_verbose = 3, bru_max_iter = 20
          )
        )
      })
      
      save(fit_pts, file = here("misc", "RData", "fit_pts.RData"))
    }
  }
})



### pts_agg --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_pts_agg.RData"))) {
    load(here::here("misc", "RData", "fit_pts_agg.RData"))
  } else {
    cmp_pts_agg <- ~ Intercept(1) + cov(nepal_rast_agg$z, model = "linear") +
      spde(main = geometry, model = matern)
    
    fml_pts_agg <- geometry ~ Intercept + cov + spde
    
    lik_pts_agg <- bru_obs(
      formula = fml_pts_agg,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )
    
    system.time({
      fit_pts_agg <- bru(
        components = cmp_pts_agg, lik_pts_agg,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    
    save(fit_pts_agg, file = here("misc", "RData", "fit_pts_agg.RData"))
  }
}


### pts_poly --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_pts_poly.RData"))) {
    load(here::here("misc", "RData", "fit_pts_poly.RData"))
  } else {
    cmp_pts_poly <- ~ Intercept(1) + cov(nepal_poly["cov"], model = "linear") +
      spde(main = geometry, model = matern)
    
    fml_pts_poly <- geometry ~ Intercept + cov + spde
    
    lik_pts_poly <- bru_obs(
      formula = fml_pts_poly,
      family = "cp",
      data = pts_samp_sf,
      domain = list(geometry = mesh_fm),
      samplers = bnd
    )
    system.time({
      fit_pts_poly <- bru(
        components = cmp_pts_poly, lik_pts_poly,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 50
        )
      )
    })
    
    save(fit_pts_poly, file = here::here("misc", "RData", "fit_pts_poly.RData"))
  }
}



## Count -------------------------------------------------------------------

### count_rast --------------------------
system.time({
  if (rerun) {
    if (file.exists(here::here("misc", "RData", "fit_count.RData"))) {
      load(here::here("misc", "RData", "fit_count.RData"))
    } else {
      cmp_count <- ~ Intercept(1) + cov(nepal_rast$z, model = "linear") +
        spde(main = geometry, model = matern)
      
      fml_count <- pts_count ~ ibm_eval(agg_nc,
                                        input = list(
                                          block = .block,
                                          weights = weight
                                        ),
                                        state = Intercept + cov + spde
      )
      
      lik_count <- bru_obs(
        formula = fml_count,
        family = "poisson",
        data = ips_nc,
        response_data = nepal_nc
      )
      
      system.time({
        fit_count <- bru(
          components = cmp_count, lik_count,
          options = list(
            control.inla = list(int.strategy = "eb"),
            bru_verbose = 3, bru_max_iter = 20
          )
        )
      })
      save(fit_count, file = here("misc", "RData", "fit_count.RData"))
    }
  }
})


### count_agg --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fit_agg.RData"))) {
    load(here::here("misc", "RData", "fit_agg.RData"))
  } else {
    cmp_agg <- ~ Intercept(1) + cov(nepal_rast_agg$z, model = "linear") +
      spde(main = geometry, model = matern)
    
    fml_agg <- pts_count ~ ibm_eval(agg_nc,
                                    input = list(
                                      block = .block,
                                      weights = weight
                                    ),
                                    state = Intercept + cov + spde
    )
    
    lik_agg <- bru_obs(
      formula = fml_agg,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )
    
    system.time({
      fit_agg <- bru(
        components = cmp_agg, lik_agg,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    #     iinla: Iteration 3 [max:5]
    #    user  system elapsed
    # 247.525   6.522 180.780
    save(fit_agg, file = here("misc", "RData", "fit_agg.RData"))
  }
}



### count_poly --------------------------
if (rerun) {
  if (file.exists(here::here("misc", "RData", "fitly.RData"))) {
    load(here::here("misc", "RData", "fitly.RData"))
  } else {
    cmpoly <- pts_count ~ Intercept(1) + cov(nepal_poly["cov"], model = "linear") +
      spde(main = geometry, model = matern)
    
    fmly <- pts_count ~ ibm_eval(agg_nc,
                                 input = list(block = .block, weights = weight),
                                 state = Intercept + cov + spde
    )
    
    likly <- bru_obs(
      formula = fmly,
      family = "poisson",
      data = ips_nc,
      response_data = nepal_nc
    )
    system.time({
      fitly <- bru(
        components = cmpoly, likly,
        options = list(
          control.inla = list(int.strategy = "eb"),
          bru_verbose = 3, bru_max_iter = 20
        )
      )
    })
    # iinla: Iteration 3 [max:5]
    #    user  system elapsed
    # 263.937   6.946 197.258
    save(fitly, file = here("misc", "RData", "fitly.RData"))
  }
}



