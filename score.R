# Score -------------------------------------------------------------------
pkgs <- c(
  "INLA", "inlabru", "fmesher", "sf", # core
  "ggplot2", "here", "dplyr", "patchwork", # general
  "terra", "stars", "tidyterra", # raster
  "purrr"
)

for (i in 1:length(pkgs)) {
  suppressMessages(library(pkgs[i], character.only = TRUE))
  # Perform other operations here
}
## complete_obs --------------------------
if (rerun) {
  mod_nm <- c(
    "fit_pts", "fit_pts_agg", "fit_pts_poly",
    "fit_count", "fit_agg", "fitly"
  )

  mod_t <- list(
    fit_pts = "Point RastFull OP",
    fit_pts_agg = "Point RastAgg OP",
    fit_pts_poly = "Point PolyAgg OP",
    fit_count = "Count RastFull OP",
    fit_agg = "Count RastAgg OP",
    fitly = "Count PolyAgg OP"
  )

  ls_dse <- sapply(mod_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse.RData"))) {
    load(here::here("misc", "RData", "ls_dse.RData"))
  } else {
    for (i in mod_nm) {
      ls_dse[[i]] <- dse_score(get(i),
        formula = lambda ~ exp(Intercept + cov + spde),
        obs = "lambda",
        newdata = pts_int, seed = rep(seed, 2)
      )
    }

    save(ls_dse, file = here::here("misc", "RData", "ls_dse.RData"))
  }


  # SE
  # lims <- range(lapply(ls_dse, function(x) as.numeric(unlist(x$pred["SE_score"])))))
  se_sum <- summary(sapply(mod_nm, function(x) {
    unlist(ls_dse[[x]][["score"]]["SE_score"])
  }))

  se_sc <- scale_fill_viridis_c(
    limits = c(0, 1e-05),
    name = "SE"
  )
  # p_pt <-  gg_subhd("Point Observation", orientation = "left-rotated", x =.05, y=.68)
  # p_ra <-  gg_subhd("RastAgg")
  # p_pa <-  gg_subhd("PolyAgg")
  # ls_p_se <- lapply(mod_nm, function(x) {
  #   if(stringr::str_detect(mod_t[[x]], "Point")){
  #   ggplot() +
  #     gg(data = ls_dse[[x]][["pred"]]["SE_score"], geom = "tile") +
  #     se_sc +
  #     labs(title=stringr::word(mod_t[[x]],-1), subtitle = mod_t[[x]])
  #   } else {
  #     ggplot() +
  #       gg(data = ls_dse[[x]][["pred"]]["SE_score"], geom = "tile") +
  #       se_sc +
  #       labs(subtitle = mod_t[[x]])
  #   }
  # })
  # p_subhd <- p_rf+p_ra+p_pa
  # p_se <- wrap_plots(ls_p_se, guides = "collect", byrow = TRUE, ncol = 3) & theme(legend.position = "right")
  # p_subhd + p_se + plot_layout(heights = c(.1, 1))
  ls_p_se <- lapply(mod_nm, function(x) {
      ggplot() +
        gg(data = ls_dse[[x]][["pred"]]["SE_score"], geom = "tile") +
        se_sc +
        ggtitle(mod_t[[x]])
  })
  wrap_plots(ls_p_se, guides = "collect", byrow = TRUE, ncol = 3) & theme(legend.position = "right")
  # & theme(legend.position = "right")& theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  # ggsave("paper/figure/SE_.png", width = tw, height = tw / 2.5, dpi = 300)
  ggsave("paper/figure/SE.png", width = tw, height = tw / 2.5, dpi = 300)
  ggsave("paper/figure/SE.pdf", width = tw, height = tw / 2.5)

  # DS
  ds_sum <- summary(sapply(mod_nm, function(x) {
    unlist(ls_dse[[x]][["score"]]["DS_score"])
  }))

  ds_sum_min <- min(sapply(mod_nm, function(x) {
    unlist(ls_dse[[x]][["score"]]["DS_score"])
  }))
  # map(ds_sum,  1)
  # min(ds_sum[1,])

  ds_sc <- scale_fill_viridis_c(
    limits = c(ds_sum_min, 0),
    name = "DS"
  )
  ls_p_ds <- lapply(mod_nm, function(x) {
    ggplot() +
      gg(data = ls_dse[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(mod_t[[x]])
  })

  wrap_plots(ls_p_ds, guides = "collect", byrow = TRUE, ncol = 3) & theme(legend.position = "right")
  ggsave("paper/figure/DS.png", width = tw, height = tw / 2.5, dpi = 300)
  ggsave("paper/figure/DS.pdf", width = tw, height = tw / 2.5)

  print(
    xtable::xtable(t(data.frame(
      point_model = as.numeric(ls_dse$fit_pts$mean_score),
      point_agg_model = as.numeric(ls_dse$fit_pts_agg$mean_score),
      point_poly_model = as.numeric(ls_dse$fit_pts_poly$mean_score),
      count_model = as.numeric(ls_dse$fit_count$mean_score),
      agg_model = as.numeric(ls_dse$fit_agg$mean_score),
      poly_model = as.numeric(ls_dse$fitly$mean_score)
    )), math.style.exponents = TRUE, digits = -3),
    file = "paper/table/agg_dse_tab.txt"
  )
}
# @

## cov_score --------------------------
if (rerun) {
  cov_nm <- c(
    "fit_ic_pts_sd01",
    "fit_ic_sd01",
    "fit_icov_sd01",
    "fit_ic_pts_sd01_nl",
    "fit_ic_sd01_nl"
  )

  cov_t <- list(
    fit_ic_pts_sd01 = "JU Point PointVal",
    fit_ic_sd01 = "JU Count PointVal",
    fit_icov_sd01 = "VP & UP Point & Count & NL PointVal",
    fit_ic_pts_sd01_nl = "JU Point NL PointVal",
    fit_ic_sd01_nl = "JU Count NL PointVal"
  )

  ls_dse_cov <- sapply(cov_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_cov.RData"))) {
    load(here::here("misc", "RData", "ls_dse_cov.RData"))
  } else {
    for (i in cov_nm) {
      ls_dse_cov[[i]] <- dse_score(get(i),
        formula = ~cov_x,
        obs = "cov",
        newdata = pts_int, seed = rep(seed, 2)
      )
    }
    save(ls_dse_cov, file = here::here("misc", "RData", "ls_dse_cov.RData"))
  }

  se_sum_cov <- summary(sapply(cov_nm, function(x) {
    unlist(ls_dse_cov[[x]][["score"]]["SE_score"])
  }))



  se_sc <- scale_fill_viridis_c(
    limits = c(0, 1e-02),
    name = "SE"
  )

  ls_p_se_cov <- lapply(cov_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_cov[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(cov_t[[x]])
  })

  wrap_plots(ls_p_se_cov, ncol = 2, guides = "collect") & theme(legend.position = "right")
  ggsave("paper/figure/SE_cov.png", width = tw, height = tw, dpi = 300)
  ggsave("paper/figure/SE_cov.pdf", width = tw, height = tw)

  # ds
  ds_sum_cov <- summary(sapply(cov_nm, function(x) {
    unlist(ls_dse_cov[[x]][["score"]]["DS_score"])
  }))
  ds_sum_cov_min <- min(sapply(cov_nm, function(x) {
    unlist(ls_dse_cov[[x]][["score"]]["DS_score"])
  }))

  ds_sc <- scale_fill_viridis_c(
    limits = c(ds_sum_cov_min, 0),
    name = "DS"
  )

  ls_p_ds_cov <- lapply(cov_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_cov[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(cov_t[[x]])
  })

  wrap_plots(ls_p_ds_cov, ncol = 2, guides = "collect") & theme(legend.position = "right")
  ggsave("paper/figure/DS_cov.png", width = tw, height = tw, dpi = 300)
  ggsave("paper/figure/DS_cov.pdf", width = tw, height = tw)

  # print(
  #   xtable::xtable(t(data.frame(
  #     point_model = as.numeric(ls_dse$fit_pts$mean_score),
  #     point_agg_model = as.numeric(ls_dse$fit_pts_agg$mean_score),
  #     point_poly_model = as.numeric(ls_dse$fit_pts_poly$mean_score),
  #     count_model = as.numeric(ls_dse$fit_count$mean_score),
  #     agg_model = as.numeric(ls_dse$fit_agg$mean_score),
  #     poly_model = as.numeric(ls_dse$fitly$mean_score)
  #   )), math.style.exponents = TRUE, digits = -3),
  #   file = "paper/table/agg_dse_tab.txt"
  # )

  print(
    xtable::xtable(
      # t(
      data.frame(
        map_dfr(
          sapply(cov_nm, function(x) {
            ls_dse_cov[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
        # )
      ),
      digits = -3
    ),
    file = "paper/table/cov_dse_tab.txt"
  )
}
# @


## cov_score_ap --------------------------
if (rerun) {
  cov_ap_nm <- c(
    "fit_ap",
    "fit_ic_ap_pts",
    "fit_ic_ap",
    "fit_ic_ap_nl",
    "fit_ic_ap_pts_nl"
  )

  cov_ap_t <- list(
    fit_ap = "VP & UP PolyAgg",
    fit_ic_ap_pts = "JU Point PolyAgg",
    fit_ic_ap = "JU Count PolyAgg",
    fit_ic_ap_pts_nl = "JU Point NL PolyAgg",
    fit_ic_ap_nl = "JU Count NL PolyAgg"
  )

  ls_dse_cov_ap <- sapply(cov_ap_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_cov_ap.RData"))) {
    load(here::here("misc", "RData", "ls_dse_cov_ap.RData"))
  } else {
    for (i in cov_ap_nm) {
      ls_dse_cov_ap[[i]] <- dse_score(get(i),
        formula = ~cov_x,
        obs = "cov",
        newdata = pts_int, seed = rep(seed, 2)
      )
    }
    save(ls_dse_cov_ap, file = here::here("misc", "RData", "ls_dse_cov_ap.RData"))
  }
  se_sum_cov_ap <- summary(sapply(cov_ap_nm, function(x) {
    unlist(ls_dse_cov_ap[[x]][["score"]]["SE_score"])
  }))



  se_sc <- scale_fill_viridis_c(
    limits = c(0, 1e-04),
    name = "SE"
  )

  ls_p_se_cov_ap <- lapply(cov_ap_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_cov_ap[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(cov_ap_t[[x]])
  })

  wrap_plots(ls_p_se_cov_ap, ncol = 2, guides = "collect") & theme(legend.position = "right")
  ggsave("paper/figure/SE_cov_ap.png", width = tw, height = tw, dpi = 300)
  ggsave("paper/figure/SE_cov_ap.pdf", width = tw, height = tw)

  # ds
  ds_sum_cov_ap <- summary(sapply(cov_ap_nm, function(x) {
    unlist(ls_dse_cov_ap[[x]][["score"]]["DS_score"])
  }))
  ds_sum_cov_ap_min <- min(sapply(cov_ap_nm, function(x) {
    unlist(ls_dse_cov_ap[[x]][["score"]]["DS_score"])
  }))

  ds_sc <- scale_fill_viridis_c(
    limits = c(ds_sum_cov_ap_min, -2),
    name = "DS"
  )

  ls_p_ds_cov_ap <- lapply(cov_ap_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_cov_ap[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(cov_ap_t[[x]])
  })

  wrap_plots(ls_p_ds_cov_ap, ncol = 2, guides = "collect") & theme(legend.position = "right")
  ggsave("paper/figure/DS_cov_ap.png", width = tw, height = tw, dpi = 300)
  ggsave("paper/figure/DS_cov_ap.pdf", width = tw, height = tw)

  # print(
  #   xtable::xtable(t(data.frame(
  #     point_model = as.numeric(ls_dse$fit_pts$mean_score),
  #     point_agg_model = as.numeric(ls_dse$fit_pts_agg$mean_score),
  #     point_poly_model = as.numeric(ls_dse$fit_pts_poly$mean_score),
  #     count_model = as.numeric(ls_dse$fit_count$mean_score),
  #     agg_model = as.numeric(ls_dse$fit_agg$mean_score),
  #     poly_model = as.numeric(ls_dse$fitly$mean_score)
  #   )), math.style.exponents = TRUE, digits = -3),
  #   file = "paper/table/agg_dse_tab.txt"
  # )

  print(
    xtable::xtable(
      # t(
      data.frame(
        map_dfr(
          sapply(cov_ap_nm, function(x) {
            ls_dse_cov_ap[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
        # )
      ),
      digits = -3
    ),
    file = "paper/table/cov_dse_tab_ap.txt"
  )
}
# @


### ju_score --------------------------

if (rerun) {
  mod_ju_nm <- c(
    "fit_ic_pts_sd01",
    "fit_ic_sd01",
    "fit_ic_ap_pts",
    "fit_ic_ap"
  )

  ju_t <- list(
    fit_ic_pts_sd01 = "JU Point PointVal",
    fit_ic_sd01 = "JU Count PointVal",
    fit_ic_ap_pts = "JU Point PolyAgg",
    fit_ic_ap = "JU Count PolyAgg"
  )

  ls_dse_ju <- sapply(mod_ju_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_ju.RData"))) {
    load(here::here("misc", "RData", "ls_dse_ju.RData"))
  } else {
    for (i in mod_ju_nm) {
      ls_dse_ju[[i]] <- dse_score(get(i),
        formula = lambda ~ exp(Intercept - beta_x * cov_x + spde),
        obs = "lambda",
        newdata = pts_int, seed = rep(seed, 2)
      )
      save(ls_dse_ju, file = here::here("misc", "RData", "ls_dse_ju.RData"))
    }
  }
  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_ju_nm, function(x) {
            ls_dse_ju[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/ju_dse_tab.txt"
  )
}


### ju_nl_score --------------------------

if (rerun) {
  mod_ju_nl_nm <- c(
    "fit_ic_pts_sd01_nl",
    "fit_ic_sd01_nl",
    "fit_ic_ap_pts_nl",
    "fit_ic_ap_nl"
  )

  ju_nl_t <- list(
    fit_ic_pts_sd01_nl = "JU Point NL PointVal",
    fit_ic_sd01_nl = "JU Count NL PointVal",
    fit_ic_ap_pts_nl = "JU Point NL PolyAgg",
    fit_ic_ap_nl = "JU Count NL PolyAgg"
  )

  ls_dse_ju_nl <- sapply(mod_ju_nl_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_ju_nl.RData"))) {
    load(here::here("misc", "RData", "ls_dse_ju_nl.RData"))
  } else {
    for (i in mod_ju_nl_nm) {
      ls_dse_ju_nl[[i]] <- dse_score(get(i),
        formula = lambda_nl ~ exp(Intercept - beta_x * cov_x + spde),
        obs = "lambda_nl",
        newdata = pts_int, seed = rep(seed, 2)
      )
      save(ls_dse_ju_nl, file = here::here("misc", "RData", "ls_dse_ju_nl.RData"))
    }
  }

  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_ju_nl_nm, function(x) {
            ls_dse_ju_nl[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/ju_nl_dse_tab.txt"
  )
}


### vp_score ----------------------------------------------------------------
if (rerun) {
  mod_vp_nm <- c(
    "fit_icov2_pts_sd01",
    "fit_icov2_sd01",
    "fit_icov2_ap_pts",
    "fit_icov2_ap"
  )

  vp_t <- list(
    fit_icov2_pts_sd01 = "VP Point PointVal",
    fit_icov2_sd01 = "VP Count PointVal",
    fit_icov2_ap_pts = "VP Point PolyAgg",
    fit_icov2_ap = "VP Count PolyAgg"
  )

  ls_dse_vp <- sapply(mod_vp_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_vp.RData"))) {
    load(here::here("misc", "RData", "ls_dse_vp.RData"))
  } else {
    for (i in mod_vp_nm) {
      ls_dse_vp[[i]] <- dse_score(get(i),
        formula = lambda ~ exp(Intercept - cov + spde),
        obs = "lambda",
        newdata = pts_int, seed = rep(seed, 2)
      )
      save(ls_dse_vp, file = here::here("misc", "RData", "ls_dse_vp.RData"))
    }
  }

  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_vp_nm, function(x) {
            ls_dse_vp[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/vp_dse_tab.txt"
  )
}
# @

### vp_nl_score ----------------------------------------------------------------
if (rerun) {
  mod_vp_nl_nm <- c(
    "fit_icov2_pts_nl_sd01",
    "fit_icov2_nl_sd01",
    "fit_icov2_ap_pts_nl",
    "fit_icov2_nl_ap"
  )

  vp_nl_t <- list(
    fit_icov2_pts_nl_sd01 = "VP Point NL PointVal",
    fit_icov2_nl_sd01 = "VP Count NL PointVal",
    fit_icov2_ap_pts_nl = "VP Point NL PolyAgg",
    fit_icov2_nl_ap = "VP Count NL PolyAgg"
  )

  ls_dse_vp_nl <- sapply(mod_vp_nl_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_vp_nl.RData"))) {
    load(here::here("misc", "RData", "ls_dse_vp_nl.RData"))
  } else {
    for (i in mod_vp_nl_nm) {
      ls_dse_vp_nl[[i]] <- dse_score(get(i),
        formula = lambda_nl ~ exp(Intercept - cov + spde),
        obs = "lambda_nl",
        newdata = pts_int, seed = rep(seed, 2)
      )
      #   save(ls_dse_vp_nl, file = here::here("misc", "RData", "ls_dse_vp_nl.RData"))
    }
  }
  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_vp_nl_nm, function(x) {
            ls_dse_vp_nl[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/vp_nl_dse_tab.txt"
  )
}
# @

# up_score ----------------------------------------------------------------
if (rerun) {
  mod_up_nm <- c(
    "fit_icov3_pts_sd01",
    "fit_icov3_sd01",
    "fit_icovap_pts",
    "fit_icovap"
  )

  up_t <- list(
    fit_icov3_pts_sd01 = "UP Point PointVal",
    fit_icov3_sd01 = "UP Count PointVal",
    fit_icovap_pts = "UP Point PolyAgg",
    fit_icovap = "UP Count PolyAgg"
  )

  ls_dse_up <- sapply(mod_up_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_up.RData"))) {
    load(here::here("misc", "RData", "ls_dse_up.RData"))
  } else {
    for (i in mod_up_nm) {
      ls_dse_up[[i]] <- dse_score(get(i),
        formula = lambda ~ exp(Intercept +
          spde -
          beta_x * (cov_pt_est + cov_uncertainty)),
        obs = "lambda",
        newdata = pts_int, seed = rep(seed, 2)
      )
      save(ls_dse_up, file = here::here("misc", "RData", "ls_dse_up.RData"))
    }
  }
  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_up_nm, function(x) {
            ls_dse_up[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/up_dse_tab.txt"
  )
}


# up_nl_score ----------------------------------------------------------------
if (rerun) {
  mod_up_nl_nm <- c(
    "fit_icov3_pts_nl_sd01",
    "fit_icov3_nl_sd01",
    "fit_icovap_pts_nl",
    "fit_icovap_nl"
  )

  up_nl_t <- list(
    fit_icov3_pts_nl_sd01 = "UP Point NL PointVal",
    fit_icov3_nl_sd01 = "UP Count NL PointVal",
    fit_icovap_pts_nl = "UP Point NL PolyAgg",
    fit_icovap_nl = "UP Count NL PolyAgg "
  )

  ls_dse_up_nl <- sapply(mod_up_nl_nm, function(x) NULL)

  if (file.exists(here::here("misc", "RData", "ls_dse_up_nl.RData"))) {
    load(here::here("misc", "RData", "ls_dse_up_nl.RData"))
  } else {
    for (i in mod_up_nl_nm) {
      ls_dse_up_nl[[i]] <- dse_score(get(i),
        formula = lambda_nl ~ exp(Intercept +
          spde -
          beta_x * (cov_pt_est + cov_uncertainty)),
        obs = "lambda_nl",
        newdata = pts_int, seed = rep(seed, 2)
      )
      save(ls_dse_up_nl, file = here::here("misc", "RData", "ls_dse_up_nl.RData"))
    }
  }

  print(
    xtable::xtable(
      data.frame(
        map_dfr(
          sapply(mod_up_nl_nm, function(x) {
            ls_dse_up_nl[[x]]["mean_score"]
          }), bind_rows,
          .id = "tib"
        )
      ),
      digits = -3
    ),
    file = "paper/table/up_nl_dse_tab.txt"
  )
}


# JUUPVP plot -------------------------------------------------------------


if (rerun) {
  ds_sum_up_min <- min(sapply(mod_up_nm, function(x) {
    unlist(ls_dse_up[[x]][["score"]]["DS_score"])
  }))

  ds_sum_vp_min <- min(sapply(mod_vp_nm, function(x) {
    unlist(ls_dse_vp[[x]][["score"]]["DS_score"])
  }))

  ds_sum_ju_min <- min(sapply(mod_ju_nm, function(x) {
    unlist(ls_dse_ju[[x]][["score"]]["DS_score"])
  }))

  # scale
  se_sc <- scale_fill_viridis_c(
    limits = c(0, 5e-06),
    name = "SE"
  )

  ds_sc <- scale_fill_viridis_c(
    limits = c(min(ds_sum_ju_min, ds_sum_vp_min, ds_sum_up_min), -10),
    name = "DS"
  )

  # JU SE
  se_sum_ju <- summary(sapply(mod_ju_nm, function(x) {
    unlist(ls_dse_ju[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_ju <- lapply(mod_ju_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_ju[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(ju_t[x]))
  })

  # DS
  ds_sum_ju <- summary(sapply(mod_ju_nm, function(x) {
    unlist(ls_dse_ju[[x]][["score"]]["DS_score"])
  }))



  ls_p_ds_ju <- lapply(mod_ju_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_ju[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(ju_t[x]))
  })

  # VP SE
  se_sum_vp <- summary(sapply(mod_vp_nm, function(x) {
    unlist(ls_dse_vp[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_vp <- lapply(mod_vp_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_vp[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(vp_t[x]))
  })

  # DS
  ds_sum_vp <- summary(sapply(mod_vp_nm, function(x) {
    unlist(ls_dse_vp[[x]][["score"]]["DS_score"])
  }))


  ls_p_ds_vp <- lapply(mod_vp_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_vp[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(vp_t[x]))
  })


  # UP SE
  se_sum_up <- summary(sapply(mod_up_nm, function(x) {
    unlist(ls_dse_up[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_up <- lapply(mod_up_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_up[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(up_t[x]))
  })

  # UP DS
  ds_sum_up <- summary(sapply(mod_up_nm, function(x) {
    unlist(ls_dse_up[[x]][["score"]]["DS_score"])
  }))

  ls_p_ds_up <- lapply(mod_up_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_up[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(up_t[x]))
  })

  # col_label_1 <- wrap_elements(panel = textGrob('Joint Uncertainty'))
  # col_label_2 <- wrap_elements(panel = textGrob('Value Plugin'))
  # col_label_3 <- wrap_elements(panel = textGrob('Uncertainty Plugin'))
  # (col_label_1/ls_p_se_ju) + plot_layout(widths = c(.1,1,1))
  # wrap_plots(c(ls_p_se_ju, ls_p_se_vp, ls_p_se_up),
  #   byrow = FALSE,
  #   ncol = 3, guides = "collect", tag_level = "new"
  # ) & theme(legend.position = "right")
  # ggsave("paper/figure/SE_juupvp2.png", width = tw, height = .8*tw, dpi = 300)

  wrap_plots(c(ls_p_se_ju, ls_p_se_vp, ls_p_se_up),
    byrow = FALSE,
    ncol = 3, guides = "collect"
  ) & theme(legend.position = "right")
  ggsave("paper/figure/SE_juupvp.png", width = tw, height = .8 * tw, dpi = 300)
  ggsave("paper/figure/SE_juupvp.pdf", width = tw, height = .8 * tw)

  wrap_plots(c(ls_p_ds_ju, ls_p_ds_vp, ls_p_ds_up),
    byrow = FALSE,
    ncol = 3, guides = "collect"
  ) & theme(legend.position = "right")
  ggsave("paper/figure/DS_juupvp.png", width = tw, height = .8 * tw, dpi = 300)
  ggsave("paper/figure/DS_juupvp.pdf", width = tw, height = .8 * tw)
}

# JUUPVP_nl plot -------------------------------------------------------------


if (rerun) {
  ds_sum_ju_nl_min <- min(sapply(mod_ju_nl_nm, function(x) {
    unlist(ls_dse_ju_nl[[x]][["score"]]["DS_score"])
  }))
  ds_sum_vp_nl_min <- min(sapply(mod_vp_nl_nm, function(x) {
    unlist(ls_dse_vp_nl[[x]][["score"]]["DS_score"])
  }))
  ds_sum_up_nl_min <- min(sapply(mod_up_nl_nm, function(x) {
    unlist(ls_dse_up_nl[[x]][["score"]]["DS_score"])
  }))

  # scale
  se_sc <- scale_fill_viridis_c(
    limits = c(0, 5e-05),
    name = "SE"
  )

  ds_sc <- scale_fill_viridis_c(
    limits = c(min(ds_sum_ju_nl_min, ds_sum_vp_nl_min, ds_sum_up_nl_min), 40),
    name = "DS"
  )

  # JU SE
  se_sum_ju_nl <- summary(sapply(mod_ju_nl_nm, function(x) {
    unlist(ls_dse_ju_nl[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_ju_nl <- lapply(mod_ju_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_ju_nl[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(ju_nl_t[x]))
  })

  # DS
  ds_sum_ju_nl <- summary(sapply(mod_ju_nl_nm, function(x) {
    unlist(ls_dse_ju_nl[[x]][["score"]]["DS_score"])
  }))


  ls_p_ds_ju_nl <- lapply(mod_ju_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_ju_nl[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(ju_nl_t[x]))
  })

  # VP SE
  se_sum_vp_nl <- summary(sapply(mod_vp_nl_nm, function(x) {
    unlist(ls_dse_vp_nl[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_vp_nl <- lapply(mod_vp_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_vp_nl[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(vp_nl_t[x]))
  })

  # DS
  ds_sum_vp_nl <- summary(sapply(mod_vp_nl_nm, function(x) {
    unlist(ls_dse_vp_nl[[x]][["score"]]["DS_score"])
  }))

  ls_p_ds_vp_nl <- lapply(mod_vp_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_vp_nl[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(vp_nl_t[x]))
  })


  # UP SE
  se_sum_up_nl <- summary(sapply(mod_up_nl_nm, function(x) {
    unlist(ls_dse_up_nl[[x]][["score"]]["SE_score"])
  }))


  ls_p_se_up_nl <- lapply(mod_up_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_up_nl[[x]][["pred"]]["SE_score"], geom = "tile") +
      se_sc +
      ggtitle(unlist(up_nl_t[x]))
  })

  # UP DS
  ds_sum_up_nl <- summary(sapply(mod_up_nl_nm, function(x) {
    unlist(ls_dse_up_nl[[x]][["score"]]["DS_score"])
  }))

  ls_p_ds_up_nl <- lapply(mod_up_nl_nm, function(x) {
    ggplot() +
      gg(data = ls_dse_up_nl[[x]][["pred"]]["DS_score"], geom = "tile") +
      ds_sc +
      ggtitle(unlist(up_nl_t[x]))
  })

  wrap_plots(c(ls_p_se_ju_nl, ls_p_se_vp_nl, ls_p_se_up_nl),
    byrow = FALSE,
    ncol = 3, guides = "collect"
  ) & theme(legend.position = "right")
  ggsave("paper/figure/SE_juupvp_nl.png", width = tw, height = .8 * tw, dpi = 300)
  ggsave("paper/figure/SE_juupvp_nl.pdf", width = tw, height = .8 * tw)

  wrap_plots(c(ls_p_ds_ju_nl, ls_p_ds_vp_nl, ls_p_ds_up_nl),
    byrow = FALSE,
    ncol = 3, guides = "collect"
  ) & theme(legend.position = "right")
  ggsave("paper/figure/DS_juupvp_nl.png", width = tw, height = .8 * tw, dpi = 300)
  ggsave("paper/figure/DS_juupvp_nl.pdf", width = tw, height = .8 * tw)
}

### subdivide --------------------------
# TODO might consider presenting the subdivide
# submesh <- fmesher:::fm_subdivide(mesh = mesh_fm,
#                              n = 1)
# submesh <- fmesher:::fm_subdivide(mesh = submesh,
#                              n = 1)
# submesh_sf <- stelfi::mesh_2_sf(submesh) %>% st_set_crs(4326)

# mesh_inner_sf <- st_crop(mesh_sf, xmin = 83, xmax = 84, ymin = 28.5, ymax = 29)
# mesh_area <- sf::st_area(mesh_inner_sf)
# range(units::set_units(mesh_area, km^2))
# mean(units::set_units(mesh_area, km^2))
# st_crs(bnd) <- NA
# ggplot() + geom_sf(data = mesh_inner_sf, fill = "white", color = "black")
### @


# eta ---------------------------------------------------------------------


if (rerun) {
  nepal_rast_int <- crop(nepal_rast, vect(bnd)[1, ], mask = TRUE)
  minmax_rast <- c(min(minmax(nepal_rast_int)), max(minmax(nepal_rast_int)))
  p_cov_rast <- ggplot() +
    gg(data = nepal_rast_int$z, na.rm = TRUE) +
    ggtitle("Covariate field") +
    scale_fill_viridis_c(
      name = "cov",
      limits = minmax_rast,
      na.value = "transparent"
    ) +
    coord_sf()

  p_cov_rast_nl <- ggplot() +
    gg(data = nepal_rast_int$znl, na.rm = TRUE) +
    ggtitle("Covariate field NL") +
    scale_fill_viridis_c(
      name = "cov",
      limits = minmax_rast,
      na.value = "transparent"
    ) +
    coord_sf()

  p_nu <- ggplot() +
    gg(data = pts_int, aes(fill = u), geom = "tile") +
    ggtitle("Matern field") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_fill_viridis_c(
      name = "Matern",
      na.value = "transparent"
    ) +
    coord_sf()

  minmax_ll <- c(
    min(min(pts_int$loglambda), min(pts_int$loglambda_nl)),
    max(max(pts_int$loglambda), max(pts_int$loglambda_nl))
  )

  p_loglambda <- ggplot() +
    gg(data = pts_int, aes(fill = loglambda), geom = "tile") +
    ggtitle("Log lambda") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_fill_viridis_c(
      name = "log_lambda",
      limits = minmax_ll,
      na.value = "transparent"
    ) +
    coord_sf()

  p_loglambda_nl <- ggplot() +
    gg(data = pts_int, aes(fill = loglambda_nl), geom = "tile") +
    ggtitle("Log lambda NL") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_fill_viridis_c(
      name = "log_lambda",
      limits = minmax_ll,
      na.value = "transparent"
    ) +
    coord_sf()

  (p_cov_rast / p_cov_rast_nl / p_nu / p_loglambda / p_loglambda_nl) + plot_layout(guides = "collect") & theme(legend.position = "right")
  ggsave(here("paper", "figure", "eta.png"),
    width = tw / 2,
    height = tw,
    dpi = 300
  )


  # (p_cov_rast / p_cov_rast_nl + plot_layout(guides = "collect")& theme(legend.position = "bottom")) |
  #   (p_nu & theme(legend.position = "bottom")) |
  #   (p_loglambda / p_loglambda_nl + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
  #
  # (p_cov_rast | p_cov_rast_nl & plot_layout(guides = "collect") & theme(legend.position = "bottom")) /
  #   (p_nu & theme(legend.position = "right")) /
  #   (p_loglambda | p_loglambda_nl + plot_layout(guides = "collect") & theme(legend.position = "right"))

  (p_cov_rast | p_cov_rast_nl) + plot_layout(guides = "collect") & theme(legend.position = "right")
  ggsave(here("paper", "figure", "eta1a.png"),
    width = tw,
    height = tw / 4,
    dpi = 300
  )
  p_nu & theme(legend.position = "right")
  ggsave(here("paper", "figure", "eta1b.png"),
    width = tw / 2,
    height = tw / 4,
    dpi = 300
  )
  (p_loglambda | p_loglambda_nl) + plot_layout(guides = "collect") & theme(legend.position = "right")
  ggsave(here("paper", "figure", "eta1c.png"),
    width = tw,
    height = tw / 4,
    dpi = 300
  )
}


# pts_obs -----------------------------------------------------------------


if (rerun) {
  p_aggcnt <- ggplot() +
    gg(
      data = nepal_nc["pts_count"],
      mapping = aes(fill = pts_count)
    ) +
    gg(pts_samp_sf, col = "red", size = .3) +
    scale_fill_continuous(name = "count") +
    ggtitle("Aggregated Count under aggregation and incomplete covariate data scenarios") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank()
    )

  p_aggcnt_nl <- ggplot() +
    gg(
      data = nepal_nc["pts_count_nl"],
      mapping = aes(fill = pts_count_nl)
    ) +
    gg(pts_samp_nl_sf, col = "red", size = .3) +
    scale_fill_continuous(name = "count_nl") +
    ggtitle("Aggregated Count NL") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank()
    )
  (p_aggcnt / p_aggcnt_nl) + plot_layout(guides = "collect")

  ggsave(
    file = here("paper", "figure", "pts_obs.png"),
    width = tw / 2,
    height = tw / 2,
    dpi = 300
  )
  ggsave(
    file = here("paper", "figure", "pts_obs.pdf"),
    width = tw / 2,
    height = tw / 2
  )
}


# ctime -------------------------------------------------------------------

if (rerun) {
  all_mod_nm <- data.frame(
    fit_pts = "Point RastFull",
    fit_pts_agg = "Point RastAgg",
    fit_pts_poly = "Point PolyAgg",
    fit_count = "Count RastFull",
    fit_agg = "Count RastAgg",
    fitly = "Count PolyAgg",
    fit_icov_sd01 = "VP and UP point+count NL PointVal",
    fit_ap = "VP and UP PolyAgg",
    fit_ic_pts_sd01 = "JU Point PointVal",
    fit_ic_sd01 = "JU Count PointVal",
    fit_ic_ap_pts = "JU Point PolyAgg",
    fit_ic_ap = "JU Count PolyAgg",
    fit_ic_pts_sd01_nl = "JU Point NL PointVal",
    fit_ic_sd01_nl = "JU Count NL PointVal",
    fit_ic_ap_pts_nl = "JU Point NL PolyAgg",
    fit_ic_ap_nl = "JU Count NL PolyAgg",
    fit_icov2_pts_sd01 = "VP Point PointVal",
    fit_icov2_sd01 = "VP Count PointVal",
    fit_icov2_ap_pts = "VP Point PolyAgg",
    fit_icov2_ap = "VP Count PolyAgg",
    fit_icov2_pts_nl_sd01 = "VP Point NL PointVal",
    fit_icov2_nl_sd01 = "VP Count NL PointVal",
    fit_icov2_ap_pts_nl = "VP Point NL PolyAgg",
    fit_icov2_nl_ap = "VP Count NL PolyAgg",
    fit_icov3_pts_sd01 = "UP Point PointVal",
    fit_icov3_sd01 = "UP Count PointVal",
    fit_icovap_pts = "UP Point PolyAgg",
    fit_icovap = "UP Count PolyAgg",
    fit_icov3_pts_nl_sd01 = "UP Point NL PointVal",
    fit_icov3_nl_sd01 = "UP Count NL PointVal",
    fit_icovap_nl = "UP Point NL PolyAgg",
    fit_icovap_pts_nl = "UP Count NL PolyAgg"
  )

  all_mod_nm_df <- data.table::transpose(all_mod_nm)
  ctime_df <- data.frame(
    sapply(names(all_mod_nm), function(x) {
      sum(get(x)[["bru_timings"]][["Time"]])
    })
  )

  ctime_df_ <- data.frame(all_mod_nm_df, ctime_df)

  names(ctime_df_) <- NULL

  print(
    xtable::xtable(
      ctime_df_
    ),
    file = "paper/table/ctime_tab.txt"
  )


  ls_ctime_id <- lapply(
    names(all_mod_nm),
    function(x) get(x)[["bru_timings"]]
  )
  names(ls_ctime_id) <- names(all_mod_nm)
  ls_ctime_id_ <- ls_ctime_id %>%
    bind_rows(.id = "type")


  # arranging order
  ls_ctime_id_$type <- factor(ls_ctime_id_$type,
    levels = unique(ls_ctime_id_$type)
  )

  ggplot(
    data = ls_ctime_id_,
    aes(Time, type, fill = Task)
  ) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Time (sec)") +
    ylab("Model")
  ggsave(here("paper", "figure", "comp_time_type.pdf"),
    width = tw, height = tw / 2
  )

  rownames(ls_ctime_id_) <- NULL
  print(
    xtable::xtable(
      ls_ctime_id_
    ),
    include.rownames = FALSE,
    file = "paper/table/ls_ctime_tab.txt"
  )
}


# z_locs ------------------------------------------------------------------

p_zlocs <- ggplot() +
  gg(data = pts_int["cov_pts"], mapping = aes(fill = cov_pts), geom = "tile") +
  gg(nepal_nc, fill = NA, col = "black") +
  gg(z_locs, col = "red", size = .3) +
  ggtitle("Sampled covariate point locations") +
  scale_fill_viridis_c(
    na.value = "transparent", name = "cov"
  )
p_zlocs
ggsave(here("paper", "figure", "z_locs.png"),
  width = tw, height = tw / 2,
  dpi = 300
)

fit_mod <- ls(pattern = "fit")
 fit_mod_cls <-  sapply(sapply(fit_mod, get), class)
 fit_mod_ls <- fit_mod[sapply(fit_mod_cls,"[[",1) == "bru"]
# convergence plot --------------------------------------------------------
for (i in fit_mod_ls) {
  bru_convergence_plot(get(i)) +
    plot_annotation(title = paste0("Convergence plots for", i))
  ggsave(
    paste0("paper/figure/bru_plot_", i, ".pdf"),
    width = 5 / 3 * tw,
    height = 5 / 3 * 3 / 4 * tw
  )
}
