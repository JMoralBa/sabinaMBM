#' Fast standardization of SpatRaster stack
#' @noRd
.standardize_rasters <- function(sp_rast, var_names = NULL, n_cores = 1) {
  
  if(is.null(var_names)) var_names <- names(sp_rast)
  if(length(var_names) == 0) return(list(rast = sp_rast, params = list()))
  
  # Parallel if n_cores > 1
  if(n_cores > 1) {
    furrr::plan(furrr::multisession, workers = n_cores, quiet = TRUE)
    on.exit(furrr::plan(furrr::sequential), add = TRUE)
    
    stats_list <- furrr::future_lapply(
      var_names,
      function(v) {
        m <- terra::global(sp_rast[[v]], "mean", na.rm = TRUE)[1, 1]
        s <- terra::global(sp_rast[[v]], "sd", na.rm = TRUE)[1, 1]
        list(mean = m, sd = s)
      },
      .options = furrr::furrr_options(seed = NULL)
    )
  } else {
    stats_list <- lapply(var_names, function(v) {
      m <- terra::global(sp_rast[[v]], "mean", na.rm = TRUE)[1, 1]
      s <- terra::global(sp_rast[[v]], "sd", na.rm = TRUE)[1, 1]
      list(mean = m, sd = s)
    })
  }
  
  names(stats_list) <- var_names
  
  # Standardize
  for(v in var_names) {
    sp_rast[[v]] <- (sp_rast[[v]] - stats_list[[v]]$mean) / stats_list[[v]]$sd
  }
  
  list(rast = sp_rast, params = stats_list)
}


# -----------------------------


#' Build inlabru likelihood objects
#' @noRd
.build_likelihoods <- function(fam, lnk, rhs_glo, rhs_reg,
                                pp_glo, pp_reg, bdy_glo, bdy_reg, dom,
                                coupling.intercept,
                                w_glo = NULL, w_reg = NULL) {
  lik_glo <- NULL
  if(fam == "cp") {
    pres_glo <- if(!is.null(coupling.intercept)) pp_glo[pp_glo$resp != 0L, ] else NULL
    pres_reg <- pp_reg[pp_reg$resp != 0L, ]
    if(!is.null(coupling.intercept)) {
      lik_glo <- inlabru::like(
        family = "cp",
        formula = as.formula(paste0("geometry ~ ", rhs_glo)),
        data = pres_glo, samplers = bdy_glo, domain = dom)
    }
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ ", rhs_reg)),
      data = pres_reg, samplers = bdy_reg, domain = dom)
  } else {
    if(!is.null(coupling.intercept)) {
      lik_glo <- inlabru::like(
        family = fam,
        formula = as.formula(paste0("resp ~ ", rhs_glo)),
        data = pp_glo, samplers = bdy_glo, domain = dom,
        weights = w_glo,
        control.family = list(link = lnk))
    }
    lik_reg <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ ", rhs_reg)),
      data = pp_reg, samplers = bdy_reg, domain = dom,
      weights = w_reg,
      control.family = list(link = lnk))
  }
  list(lik_glo = lik_glo, lik_reg = lik_reg)
}


# -----------------------------


#' Components GL/RE for covariates
#' @noRd
.fcov <- function(obj, 
                 spobjglo, 
                 spobjreg,
                 sp_covglo,
                 sp_covreg,
                 covariate.effects = NULL,
                 coupling.predictors = NULL,
                 pp_glo_sf = NULL,
                 pp_reg_sf= NULL) {

  # selected vars
  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional

  #cmp1 <- paste0(unique(c(vg, vr)), "(1)", collapse = " + ")

  # rw2: knots to enforce min relative spacing (INLA check 1e-3) 
  thin_knots <- function(x, min_ratio = 1e-3) {
    x <- sort(unique(as.numeric(x)))
    if(length(x) <= 2) return(x)
    r <- diff(range(x))
    if(!is.finite(r) || r == 0) return(unique(x))
    out <- x[1]
    for(xi in x[-1]) {
      if((xi - out[length(out)]) / r >= min_ratio) { 
        out <- c(out, xi)
      }
    }
    if(length(out) < 3L) {
      out <- seq(min(x), max(x), length.out = min(max(5L, length(x)), 50L))
    }
    unique(out)
  }

  #rw2 defaults (inla recommends quantile grouping and K = 150-300 (balance stability and flexibility)
  rw2_K <- 300L              #@@@JMB pensar si dejamos esto por defecto
  rw2_method <- "quantile" 

  .build_rw2_values <- function(rast_layer, coords, K = rw2_K, method = rw2_method) {
    vals <- suppressWarnings(as.numeric(terra::extract(rast_layer, coords)[, 1]))
    vals <- vals[is.finite(vals)] 
    rng <- range(vals)
    if(diff(rng) == 0) {
      stop("❌  RW2 requires variability (covariate range = 0).\n",
           "   Use 'linear' instead for this covariate or check the raster values.\n\n")
    }  
    g <- INLA::inla.group(vals, n = K, method = method)
    v <- sort(unique(as.numeric(levels(g))))
    v <- thin_knots(v)
    if(length(v) < 3L) {
      stop("❌  RW2 requires ≥ 3 distinct support values.\n", #@@@JMB aquí también se podría ajustar K, pero demasiados args en mi opinión
           "   Consider reducing smoothing or check covariate variability.\n\n")
    }
    v
  }

  # detect rw2 needs
  need_rw2_global <- FALSE
  need_rw2_regional <- FALSE
  if(is.list(covariate.effects)) {
    if(length(vg) > 0) {
      need_rw2_global <- any(vapply(vg, function(x) 
        .resolve_covariate_effects(x, "global",  covariate.effects)$model == "rw2", logical(1)))
    }
    if(length(vr) > 0) {
      need_rw2_regional <- any(vapply(vr, function(x) 
        .resolve_covariate_effects(x, "regional", covariate.effects)$model == "rw2", logical(1)))
    }
  }

  coords_all <- if(need_rw2_global) {
    rbind(sf::st_coordinates(pp_glo_sf),
          sf::st_coordinates(pp_reg_sf))
  } else NULL
  coords_r <- if(need_rw2_regional) {
    sf::st_coordinates(pp_reg_sf)
  } else NULL

  # global components 
  cmpglobal <- character(0)
  fglobal <- character(0)

  fmt_vals <- function(v) {
    paste0("c(", paste(format(v, digits = 7), collapse = ","), ")")
  }

  for(X in vg) {
    cp_mode <- .resolve_coupling_predictor(X, coupling.predictors, vg)
    if(cp_mode == "NULL") next
    specX <- .resolve_covariate_effects(X, "global", covariate.effects)
    if(specX$model == "drop") next
    if(specX$model == "linear") {
      cmpglobal <- c(cmpglobal,
        paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'linear')"))
    } else if(specX$model == "rw2") {
      vals <- .build_rw2_values(sp_covglo[[X]], coords_all)
      cmpglobal <- c(cmpglobal,
        paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'rw2', scale.model = TRUE, ",
               "hyper = list(prec = list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))), ",
               "values = ", fmt_vals(vals),")"))
    }
    fglobal <- c(fglobal, paste0(X, "GL"))
  }

  # regional components 
  cmpregional <- character(0)
  fregional <- character(0)

  for(X in vr) {
    #cp_mode <- .resolve_coupling_predictor(X, coupling.predictors, vg)
    is_shared <- X %in% vg
    cp_mode <- if(is_shared) .resolve_coupling_predictor(X, coupling.predictors, vg) else "unpooled"

    specX <- .resolve_covariate_effects(X, "regional", covariate.effects)
    if(specX$model == "drop") next

    ## scale_decomposed: global (XGL_reg) + regional anomaly (XRE_delta)
    if(cp_mode == "scale_decomposed") {
      if(specX$model == "linear") {
        cmpregional <- c(cmpregional,
          paste0(X, "GL_ls(main = ", spobjreg, ", main_layer = '", X, "_ls', model = 'linear')"),
          paste0(X, "RE_ss(main = ", spobjreg, ", main_layer = '", X, "_ss', model = 'linear')"))
      } else if(specX$model == "rw2") {
        vals <- .build_rw2_values(sp_covreg[[paste0(X, "_ss")]], coords_r)
        cmpregional <- c(cmpregional,
          paste0(X, "GL_ls(main = ", spobjreg, ", main_layer = '", X, "_ls', model = 'linear')"),
          paste0(X, "RE_ss(main = ", spobjreg, ", main_layer = '", X, "_ss', model='rw2', scale.model=TRUE, values = ", fmt_vals(vals), ", hyper = list(prec=list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))))"))
      }
      fregional <- c(fregional, paste0(X, "GL_ls"), paste0(X, "RE_ss"))
    }

    # ordered_hierarchical: soft constraint beta_RE ~ N(1, 0.5^2)
    # mean.linear=1 with standardized covariates encourages regional effect to mirror global scale.
    #@@@JMB PENDIENTE consultar con Virgilio: La version A con prior fijo N(1, 0.5^2) funciona pero es soft constraint????, no jerarquía real. 
       # La versión B hace copy sobre componente lineal global es jerarquía real beta_RE|beta_GL~N(beta_GL,tau) pero se rompe. Pendiente verificar si es por strategy eb o por linear effects o q????
    else if(cp_mode == "ordered_hierarchical") {
      # version A
      cmpregional <- c(cmpregional,
        paste0(X, "RE_oh(main = ", spobjreg, ", main_layer = '", X, "', model = 'linear', mean.linear = 1, prec.linear = 4)"))
      # version B
      # cmpregional <- c(cmpregional,
      #   paste0(X, "RE_oh(main = ", spobjreg, ", main_layer = '", X, "', ",
      #          "copy = '", X, "GL', fixed = FALSE, ",
      #          "hyper = list(beta = list(prior = 'normal', param = c(1, 4))))"))
      fregional <- c(fregional, paste0(X, "RE_oh"))
    }

    # bayesian_feedback: global priors
    else if(cp_mode == "bayesian_feedback") {
      cmpregional <- c(cmpregional,
        paste0(X, "RE(main = ", spobjreg, ", main_layer = '", X, "', model = 'linear', mean.linear = bf_mean_", X, ", prec.linear = bf_prec_", X, ")"))
      fregional <- c(fregional, paste0(X, "RE"))
    }

    # unpooled/NULL: coefs independent
    else {
      if(specX$model == "linear") {
        cmpregional <- c(cmpregional, paste0(X, "RE(main = ", spobjreg, ", main_layer = '", X, "', model = 'linear')"))
      } else if(specX$model == "rw2") {
        vals <- .build_rw2_values(sp_covreg[[X]], coords_r)
        cmpregional <- c(cmpregional,
          paste0(X, "RE(main = ", spobjreg, ", main_layer = '", X, "', model = 'rw2', scale.model = TRUE, ",
                 "hyper = list(prec = list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))), ",
                 "values = ", fmt_vals(vals),")"))
      }
      fregional <- c(fregional, paste0(X, "RE"))
    }
  }

  # out fcov
  list(
    #cmp = paste(c(cmp1, cmpglobal, cmpregional), collapse = " + "),
    cmp = paste(c(cmpglobal, cmpregional), collapse = " + "),
    like = list(fglobal = if(length(fglobal) > 0) paste(fglobal, collapse = " + ") else "",
                fregional = if(length(fregional) > 0) paste(fregional, collapse = " + ") else "")) 
}


# -----------------------------


#' Fit NSBM sequentially or jointly
#' @noRd
.fit_nsbm <- function(cmp, lik_list, coupling.intercept, 
                      coupling.predictors, needs_feedback, 
                      vr = NULL, n.threads = 1, seed = NULL, 
                      int.strategy = "eb") {
  bru_opts <- list(
    control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
    control.inla = list(int.strategy = int.strategy),
    control.mode = list(restart = TRUE),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
    num.threads = n.threads
  )
  
  if(needs_feedback) {
    env_cmp <- environment(cmp)
    
    # P0. dummies for phase 1 to prevent crashes
    assign("bf_mean_int", 0, envir = env_cmp)
    assign("bf_prec_int", 1, envir = env_cmp)
    if(!is.null(vr)) {
      for(v in vr) {
        assign(paste0("bf_mean_", v), 0, envir = env_cmp)
        assign(paste0("bf_prec_", v), 1, envir = env_cmp)
      }
    }

    # P1. fit global model using ONLY the global likelihood
    message("\nℹ️  Sequential bayesian feedback - Fitting global model to extract posteriors...")
    fit_glo <- do.call(inlabru::bru, c(list(components = cmp), list(lik_list[[1]]), list(options = bru_opts)))

    # P2. extract moments and inject them into the current environment
    message("    Updating regional priors by moments and fitting joint model...\n")

    # helper: safe precision with floor to avoid Inf or near-0 variance
    # floor at 1e-4 precision (sd ~ 100) prevents degenerate flat priors
    # ceiling at 1e6 precision (sd ~ 0.001) prevents point-mass priors
    .safe_prec <- function(sd_val, floor_prec = 1e-4, ceil_prec = 1e6) {
      prec <- 1 / (sd_val^2)
      prec <- pmax(floor_prec, pmin(ceil_prec, prec))
      prec
    }

    # helper: warn if global marginal is strongly skewed (moment matching unreliable)
    .check_skewness <- function(marginal, label) {
      if(is.null(marginal)) return(invisible(NULL))
      sm <- try(INLA::inla.smarginal(marginal), silent = TRUE)
      if(inherits(sm, "try-error")) return(invisible(NULL))
      m  <- INLA::inla.emarginal(function(x) x,       marginal)
      m2 <- INLA::inla.emarginal(function(x) (x-m)^2, marginal)
      m3 <- INLA::inla.emarginal(function(x) (x-m)^3, marginal)
      skew <- m3 / (m2^(3/2))
      if(abs(skew) > 1) {
        warning("⚠️ Bayesian feedback: global posterior for '", label,
                "' has skewness = ", round(skew, 2),
                " — moment matching (mean/sd) may be a poor approximation.\n",
                "   Consider using `coupling.intercept = 'ordered_hierarchical'` instead.\n\n")
      }
    }

    if(!is.null(coupling.intercept) && coupling.intercept == "bayesian_feedback") {
      int_random <- fit_glo$summary.random$IGlobal
      if(is.null(int_random) || nrow(int_random) != 1L) {
        stop("❌ bayesian_feedback: could not extract IGlobal posterior (unexpected structure).\n\n")
      }
      bf_mean_int <- int_random$mean
      bf_sd_int   <- int_random$sd
      .check_skewness(fit_glo$marginals.random$IGlobal[[1]], "IGlobal")
      assign("bf_mean_int", bf_mean_int,              envir = env_cmp)
      assign("bf_prec_int", .safe_prec(bf_sd_int),    envir = env_cmp)
    }

    if(!is.null(fit_glo$summary.fixed)) {
      gl_effs <- rownames(fit_glo$summary.fixed)
      gl_effs <- gl_effs[grepl("GL$", gl_effs)]
      for(eff in gl_effs) {
        base_var <- sub("GL$", "", eff)
        bf_sd_v  <- fit_glo$summary.fixed[eff, "sd"]
        .check_skewness(fit_glo$marginals.fixed[[eff]], eff)
        assign(paste0("bf_mean_", base_var), fit_glo$summary.fixed[eff, "mean"], envir = env_cmp)
        assign(paste0("bf_prec_", base_var), .safe_prec(bf_sd_v),                envir = env_cmp)
      }
    }
    
    # final fit using ONLY regional likelihood
    fit <- do.call(inlabru::bru, c(list(components = cmp), list(lik_list[[length(lik_list)]]), list(options = bru_opts)))
    return(fit)
    
  } else { # standard joint fit
    fit <- do.call(inlabru::bru, c(list(components = cmp), lik_list, list(options = bru_opts)))
    return(fit)
  }
}



# -----------------------------


#' interprete covariate_effects
#' @noRd
.resolve_covariate_effects <- function(varname, scale, covariate.effects, default_model = "linear") {

  path_label <- paste0("covariate.effects$", scale, "$", varname)

  # if no covariate.effects, all linear by default
  if(is.null(covariate.effects) || !is.list(covariate.effects)) {
    return(list(model = default_model, u = NA, alpha = NA))
  }

  # global/regional or default
  if(!is.null(covariate.effects[[scale]]) && !is.null(covariate.effects[[scale]][[varname]])) {
    spec <- covariate.effects[[scale]][[varname]]
  } else if(!is.null(covariate.effects$default)) {
    spec <- covariate.effects$default
    path_label <- paste0("covariate.effects$default (applied to ", varname, ")")
  } else {
    return(list(model = "linear", u = NA, alpha = NA))
  }

  # if spec is a string -> "linear"/"drop"
  if(is.character(spec)) {
    if(spec == "linear") { 
      return(list(model = "linear", u = NA, alpha = NA))
    }
    if(spec == "drop") { 
      return(list(model = "drop", u = NA, alpha = NA))
    } 
    if(spec == "rw2") {
      stop("❌  Invalid RW2 specification in `", path_label, "`.\n",
           "   RW2 must be expressed as a list: list(model='rw2', u=..., alpha=...).\n",
           "   For linear effects use 'linear'; to exclude use 'drop'.\n\n")
    } 
    stop("❌ Invalid keyword '", spec, "' in `", path_label, "`.\n",
         "   Valid options: 'linear', 'drop', or list(model='rw2', u=..., alpha=...) (smoothed RW2 effect).\n\n")
  }

  # if spec is a list, only rw2 admited 
  if(is.list(spec)) {
    if(is.null(spec$model)) {
      stop("❌  Invalid list specification in `", path_label, "`.\n",
           "   Lists are only allowed for RW2 effects. Use: list(model='rw2', u=..., alpha=...).\n",
           "   For linear effects use 'linear'; to exclude use 'drop'.\n\n")
    }
    if(identical(spec$model, "rw2")) {
      # RW2 reequires u alpha
      if(is.null(spec$u) || is.null(spec$alpha)) {
        stop("❌  Incomplete RW2 specification in `", path_label, "`.\n",
             "   Provide both `u` and `alpha`..\n\n")
      }
      return(list(model = "rw2", u = spec$u, alpha = spec$alpha))
    }
    stop("❌  Invalid `model` in `", path_label, "`.\n",
         "   When using a list, only model='rw2' is permitted.\n",
         "   For linear effects use 'linear'; to exclude use 'drop'.\n\n")
  }

  #unsupported type
  stop("❌  Unsupported type in `", path_label, "`.\n",
       "   Must be string ('linear'|'drop') or list(model='rw2', u=..., alpha=...).\n\n")
}


# -----------------------------


#' interpret coupling.predictors
#' @noRd
.resolve_coupling_predictor <- function(var, coupling.predictors, vg = NULL) {

  if(is.null(coupling.predictors)) return("unpooled")

  mode_val <- "unpooled"
  is_explicit <- FALSE

  if(is.character(coupling.predictors)) {
    mode_val <- coupling.predictors
  } else if(is.list(coupling.predictors)) {
    if(!is.null(coupling.predictors$variables) && !is.null(coupling.predictors$variables[[var]])) {
      mode_val <- coupling.predictors$variables[[var]]
      is_explicit <- TRUE
    } else if(!is.null(coupling.predictors$default)) {
      mode_val <- coupling.predictors$default
    }
  }
  
  # unpooled for vars only in regional
  if(!is.null(vg) && !(var %in% vg) && !is_explicit) {
    if(mode_val %in% c("ordered_hierarchical", "scale_decomposed", "bayesian_feedback")) {
      mode_val <- "unpooled"
    }
  }
  
  return(mode_val)
}


# -----------------------------


#' from sf to tif
#' @noRd
.pred_as_tif <- function(pred, template, vars_to_export = c("mean", 
                                                           "sd", 
                                                           "q0.025", 
                                                           "q0.5", 
                                                           "q0.975",
                                                           "median",
                                                           "sd.mc_std_err",
                                                           "mean.mc_std_err")) {

  stopifnot(inherits(pred, c("bru_prediction", "sf")))
  stopifnot(inherits(template, "SpatRaster"))

  r_stack <- lapply(vars_to_export, function(var) {
    if(!is.null(pred[[var]])) {
      df <- as.data.frame(cbind(sf::st_coordinates(pred), value = pred[[var]]))
      sf_pts <- sf::st_as_sf(df, coords = c("X", "Y"), crs = sf::st_crs(pred))
      sf_pts <- sf::st_transform(sf_pts, terra::crs(template))
      r <- terra::rasterize(sf_pts, template, field = "value")
      names(r) <- var
      return(r)
    } else {
      return(NULL)
    }
  })

  r_stack <- Filter(function(x) inherits(x, "SpatRaster"), r_stack)
  if(length(r_stack) == 0) return(NULL)

  r_pred <- do.call(c, r_stack)
  return(r_pred)
}


# -----------------------------


#' random distribution of data in k-folds
#' @noRd
.make_stratified_kfolds <- function(y, K) {
  # y: vector de 1/0
  # K: número de folds
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  f <- integer(length(y))
  f[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  f[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  f
}


# -----------------------------


#' diagnostics
#' @noRd
.nsbm_diagnostics <- function(fit,
                             fam,
                             data_used,
                             n_glo = 0L,
                             priors = NULL,
                             pred_local = NULL,
                             pred_shared = NULL,
                             coupling.intercept,
                             scale_params = NULL) {

  has_Sloc <- "Sloc" %in% names(fit$summary.random)
  has_Sshared <- "Sshared" %in% names(fit$summary.random)

  # Model fit metrics
  dic_val <- if(!is.null(fit$dic$dic) && !is.na(fit$dic$dic)) round(fit$dic$dic, 2) else NA_real_
  waic_val <- if(!is.null(fit$waic$waic) && !is.na(fit$waic$waic)) round(fit$waic$waic, 2) else NA_real_
  mlik_val <- if(!is.null(fit$mlik) && !is.na(fit$mlik[1,1])) round(fit$mlik[1, 1], 2) else NA_real_
  mlpd_val <- if(!is.null(fit$cpo) && !is.null(fit$cpo$cpo)) round(mean(log(fit$cpo$cpo), na.rm = TRUE), 3) else NA_real_
  lcpo_val <- if(!is.null(fit$cpo$cpo)) round(sum(log(fit$cpo$cpo)), 2) else NA_real_

  hyper <- fit$summary.hyperpar
  marginals <- fit$marginals.hyperpar

  ci_ratio <- function(param) {
    if(!param %in% rownames(hyper)) return("—")
    ciw <- hyper[param, "0.975quant"] - hyper[param, "0.025quant"]
    med <- hyper[param, "0.5quant"]
    if(!is.finite(ciw) || !is.finite(med) || med == 0) return(NA_real_)
    ciw / abs(med)
  }

  # hyperparameters spde
  range_res_mean <- if(has_Sloc) hyper["Range for Sloc", "mean"] else "—"
  sigma_res_mean <- if(has_Sloc) hyper["Stdev for Sloc", "mean"] else "—"
  range_lat_mean <- if(has_Sshared) hyper["Range for Sshared", "mean"] else "—"
  sigma_lat_mean <- if(has_Sshared) hyper["Stdev for Sshared", "mean"] else "—"
  range_res_sd <- if(has_Sloc) hyper["Range for Sloc", "sd"] else "—"
  sigma_res_sd <- if(has_Sloc) hyper["Stdev for Sloc", "sd"] else "—"
  range_lat_sd <- if(has_Sshared) hyper["Range for Sshared", "sd"] else "—"
  sigma_lat_sd <- if(has_Sshared) hyper["Stdev for Sshared", "sd"] else "—"

  range_res <- if(has_Sloc) paste0(round(range_res_mean, 2), " ± ", round(range_res_sd, 2)) else "—"
  sigma_res <- if(has_Sloc) paste0(round(sigma_res_mean, 2), " ± ", round(sigma_res_sd, 2)) else "—"
  range_lat <- if(has_Sshared) paste0(round(range_lat_mean, 2), " ± ", round(range_lat_sd, 2)) else "—"
  sigma_lat <- if(has_Sshared) paste0(round(sigma_lat_mean, 2), " ± ", round(sigma_lat_sd, 2)) else "—"
  format_hyper <- function(param) {
   if(!param %in% rownames(hyper)) return("—")
      paste0(round(hyper[param, "mean"], 2), " ± ", round(hyper[param, "sd"], 2))
  }
  prec_IGlobal <- if("Precision for IGlobal" %in% rownames(hyper)) format_hyper("Precision for IGlobal") else "—"
  beta_IRegional <- if("Beta for IRegional" %in% rownames(hyper)) format_hyper("Beta for IRegional") else "—"
  # betas from ordered_hierarchical predictor copies: "Beta for bio1RE_oh", etc.
  beta_pred_oh <- if(!is.null(hyper)) {
    oh_rows <- rownames(hyper)[grepl("^Beta for .+RE_oh$", rownames(hyper))]
    if(length(oh_rows) > 0) {
      stats::setNames(
        lapply(oh_rows, format_hyper),
        oh_rows
      )
    } else NULL
  } else NULL

  # Sloc and Sshared overlap?
  if(is.finite(range_res_mean) && is.finite(range_lat_mean)) {
    ratio <- range_lat_mean / range_res_mean
  }

  # Intercepts
  mean_sd_str <- function(df) {
    if(is.null(df) || nrow(df) == 0) return("—")
    paste0(round(df$mean[1], 5), " ± ", round(df$sd[1], 5), " (", 
           round(df$`0.025quant`[1], 5), ", ", round(df$`0.975quant`[1], 5), ")")
  }
  iGlobal <- mean_sd_str(fit$summary.random$IGlobal)
  iRegional <- mean_sd_str(fit$summary.random$IRegional)

  # Significant vars
  sig_vars <- .signif_vars(fit, scale_params = scale_params)

  if (fam == "cp") {
    auc_full <- "—"; tjur_r2 <- "—"; brier <- "—"; rmse <- "—"; pred_cor <- "—"
    cal_slope <- "—"; ks_pit <- "—"; cov50 <- "—"; cov95 <- "—"; pit_reg <- NULL; moran_I <- "—"
    rs <- NULL

  } else {
    y_obs <- data_used$resp
    idx_reg <- seq(n_glo + 1L, n_glo + nrow(data_used))
    y_pred <- fit$summary.fitted.values$mean[idx_reg]

    # auc_full and Tjur R2 (if binary)
    if(all(y_obs %in% c(0,1))) {
      auc_full <- suppressMessages(as.numeric(pROC::auc(y_obs, y_pred)))
      tjur_r2  <- mean(y_pred[y_obs == 1L], na.rm = TRUE) - mean(y_pred[y_obs == 0L], na.rm = TRUE)
    } else {
      auc_full <- "—"
      tjur_r2  <- "—"
    }

    # brier, pred correlation, RMSE
    brier <- mean((y_pred - y_obs)^2, na.rm = TRUE)
    pred_cor <- stats::cor(y_obs, y_pred, use = "complete.obs")
    rmse <- sqrt(mean((y_obs - y_pred)^2, na.rm = TRUE))

    # pit
    pit_vals <- fit$cpo$pit
    if(!is.null(pit_vals)) {
      pit_reg <- pit_vals[idx_reg]
      pit_reg <- pit_reg[is.finite(pit_reg)]
      ks_pit <- if(length(pit_reg) > 3) suppressWarnings(ks.test(pit_reg, "punif")$p.value) else "—"
    } else {
      ks_pit <- "—"
    }

    # Moran’s I residual autocorrelation
    moran_I <- NA_real_
    rs <- y_obs - y_pred
    xy <- as.matrix(data_used[, c("x", "y")])
    maxdist <- if(has_Sloc) {
      range_res_mean
    } else if(has_Sshared) {
      range_lat_mean
    } else {               #@@@JMB sin Sloc ni Sshared usa 1/4 de la diagonal???
      bb <- apply(xy, 2, range, na.rm = TRUE)
      sqrt(sum((bb[2,] - bb[1,])^2)) / 4
    }   
    if(is.finite(maxdist) && maxdist > 0) {
      nb <- spdep::dnearneigh(xy, 0, maxdist, longlat = FALSE)
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      mi <- spdep::moran(rs, lw, n = length(rs), S0 = spdep::Szero(lw))
      moran_I <- as.numeric(mi$I)
    }

    # calibration & coverage (disable for binomial, enable for continuous)
    if (fam == "binomial") {
      cal_slope <- "—"
      ks_pit <- "—"
      cov50 <- "—"
      cov95 <- "—"
    } else { # for continuous families (gaussian,...)
      if(all(y_pred > 0 & y_pred < 1)) {
        df_cal <- data.frame(
          logit_p = qlogis(y_pred),
          y = y_obs
        )
        cal_mod <- lm(logit_p ~ y, data = df_cal)
        cal_slope <- coef(cal_mod)[2]
      } else {
        cal_slope <- "—"
      }

      # pit ks-test
      if(!is.null(pit_reg) && length(pit_reg) > 3) {
        ks_pit <- suppressWarnings(ks.test(pit_reg, "punif")$p.value)
      } else {
        ks_pit <- "—"
      }

      # coverage 95%
      idx_reg <- seq(n_glo + 1L, n_glo + nrow(data_used))
      if(all(c("0.025quant", "0.975quant") %in% colnames(fit$summary.fitted.values))) {
        low95 <- fit$summary.fitted.values[idx_reg, "0.025quant"]
        up95  <- fit$summary.fitted.values[idx_reg, "0.975quant"]
        cov95 <- mean(y_obs >= low95 & y_obs <= up95, na.rm = TRUE)
      } else {
        cov95 <- "—"
      }

      # coverage 50%
      if(all(c("0.25quant", "0.75quant") %in% colnames(fit$summary.fitted.values))) {
        low50 <- fit$summary.fitted.values[idx_reg, "0.25quant"]
        up50  <- fit$summary.fitted.values[idx_reg, "0.75quant"]
        cov50 <- mean(y_obs >= low50 & y_obs <= up50, na.rm = TRUE)
      } else {
        cov50 <- "—"
      }
    }
  }

  # CI/median ratios
  ci_ratios <- c(range_Sloc = ci_ratio("Range for Sloc"),
                 range_Sshared = ci_ratio("Range for Sshared"),
                 sigma_Sloc = ci_ratio("Stdev for Sloc"),
                 sigma_Sshared = ci_ratio("Stdev for Sshared"))

  # credibility interval ratios (CI/median)
  # High values indicate weak identifiability or non-informative priors.
  max_CIratio <- if(any(is.finite(ci_ratios))) max(ci_ratios, na.rm = TRUE) else NA_real_

  # posterior ≈ prior
  close_rel <- function(post_med, prior_u, tol = 0.15) {  #@@@JMB rev threshold = 15% difference en Bakka et al. 2018)
    if(!is.finite(post_med) || is.null(prior_u) || length(prior_u) < 1 ||
       !is.finite(prior_u[1]) || prior_u[1] == 0) return(FALSE)
    abs(post_med - prior_u[1]) / abs(prior_u[1]) < tol
  }

  prior_close <- list(
    range_Sloc = if(!is.null(priors$local.pcprior.range))
    close_rel(range_res_mean, priors$local.pcprior.range) else FALSE,
  sigma_Sloc = if(!is.null(priors$local.pcprior.sigma))
    close_rel(sigma_res_mean, priors$local.pcprior.sigma) else FALSE,
  range_Sshared = if(!is.null(priors$shared.pcprior.range))
    close_rel(range_lat_mean, priors$shared.pcprior.range) else FALSE,
  sigma_Sshared = if(!is.null(priors$shared.pcprior.sigma))
    close_rel(sigma_lat_mean, priors$shared.pcprior.sigma) else FALSE
  )

  # variance and range ratios
  # sigma Sshared mean / sigma Sloc mean > 1.5 ==> Sshared field dominates (Blangiardo & Cameletti 2015)
  sigma_ratio <- if(has_Sloc && has_Sshared && is.finite(sigma_lat_mean) && is.finite(sigma_res_mean) && sigma_res_mean > 0)
    sigma_lat_mean / sigma_res_mean else "—"
  # 0.5 < range_Sshared / range_Sloc < 2 ==> poor scale separation (Bakka et al. 2018)
  range_ratio <- if(has_Sloc && has_Sshared && is.finite(range_lat_mean) && is.finite(range_res_mean) && range_res_mean > 0)
    range_lat_mean / range_res_mean else "—"

var_explained_sloc <- if(has_Sloc && has_Sshared) {
  var_sloc_contrib <- INLA::inla.emarginal(function(x) x^2, marginals$`Stdev for Sloc`)
  var_sshared_contrib <- INLA::inla.emarginal(function(x) x^2, marginals$`Stdev for Sshared`)
  var_sloc_contrib / (var_sloc_contrib + var_sshared_contrib)
} else "—"
ssi <- if(is.finite(range_ratio) && is.finite(sigma_ratio)) {
  # SSI ideal: range_ratio > 3 AND sigma_ratio ≈ 1–2
  range_penalty <- if(range_ratio > 3) 1 else range_ratio / 3
  sigma_penalty <- min(sigma_ratio, 1/sigma_ratio)  # simétrico
  range_penalty * sigma_penalty
} else NA_real_
sri <- if(has_Sloc && has_Sshared) {
  resid_sloc <- fit$summary.random$Sloc$mean
  resid_sshared <- fit$summary.random$Sshared$mean
  n <- min(length(resid_sloc), length(resid_sshared))
  cor_fields <- cor(resid_sloc[1:n], resid_sshared[1:n], use = "complete.obs")
  cor_fields^2  # R² de correlación
} else "—"

  # correlation Sloc–Sshared fields
  # r > 0.7 00> high correlation
  field_correlation <- "—"
  if(has_Sloc && has_Sshared) {
    f_sp <- fit$summary.random$Sloc$mean
    f_lat <- fit$summary.random$Sshared$mean
    n <- min(length(f_sp), length(f_lat))
    field_correlation <- stats::cor(f_sp[seq_len(n)], f_lat[seq_len(n)], use = "pairwise.complete.obs")
  } 

  # warnings
  warns <- character()
  # overlap Sshared vs Sloc
  if(has_Sloc && has_Sshared && is.finite(range_ratio) && range_ratio > 0.5 && range_ratio < 3) {
    warns <- c(warns,
      paste("⚠️  Insufficient Sloc scale separation: Sshared/Sloc range ratio = ",round(range_ratio, 2), ").\n",
             "   Scales are overlapping (ideal: ratio > 3–5).\n",
             "   Recomended actions: (1) Tighten local.pcprior.range[1] (e.g., 2 → 1.5),\n",
             "   (2) Increase shared.pcprior.range[1] (e.g., 50 → 100),\n",
             "   (3) Run prior sensitivity analysis,\n",
             "   (4) Evaluate if study domain truly supports two scales.\n")
    )
  }  
  # weak identifiability (CI/median ratio)
  if(is.finite(max_CIratio) && max_CIratio > 25) {
    warns <- c(warns, paste0(
      "⚠️ Weak hyperparameter identifiability detected (CI/median > 25).\n",
      "   Recommended: use stronger PC priors."
    ))
  }
  # Sshared field dominating variance
  if(has_Sshared && has_Sloc && is.finite(sigma_ratio) && sigma_ratio > 1.5) {
    warns <- c(warns, paste0(
      "⚠️ Variance imbalance between Sshared and Sloc SPDE. Variance ratio (sigma Sshared / sigma Sloc) = ", round(sigma_ratio, 2), ".\n",
      "   The Sshared field dominates the Sloc variability, making the Sloc SPDE redundant.\n",
      "   Recommended: reduce shared.pcprior.sigma or increase local.pcprior.sigma."
    ))
  }
  # high correlation between Sshared and Sloc
  if(has_Sloc && has_Sshared && is.finite(field_correlation) && field_correlation > 0.7) {
    warns <- c(warns,
      paste0("⚠️ High correlation between Sshared and Sloc fields (r = ", round(field_correlation, 2), ").\n",
             "   Possible redundancy: both SPDE components capture the same spatial pattern.\n",
             "   Recommended: strengthen priors to separate scales, or remove one SPDE fields.")
    )
  }
  # residual autocorrelation
  if(is.finite(moran_I) && moran_I > 0.10) {
    warns <- c(warns,
      paste0("⚠️ Residual spatial autocorrelation detected (Moran’s I ≈ ", round(moran_I, 2), ").\n",
             "   Model missing local spatial structure.\n",
             "   Recommended: add a local SPDE component, refine mesh resolucion (smaller max.edges), or include missing covariates.")  #@@@JMB no estoy segura
    )
  }
  # posterior ≈ prior (weak data information)
  if(any(unlist(prior_close))) {
    warns <- c(warns, paste0(
      "⚠️ Posterior close to PC-prior mode: weak data information relative to prior strength.\n",
      "   Recommended: relax priors or increase data resolution."    #@@@JMB rev recommendation??
    ))
  }    
  # cpo
  if(!is.null(fit$cpo$cpo)) {
    cpo_failures <- sum(fit$cpo$failure, na.rm = TRUE)
    total_obs <- length(fit$cpo$cpo)
    if(cpo_failures > 0) {
      perc_failures <- (cpo_failures / total_obs) * 100
      if(perc_failures > 1) {     #@@@JMB 1% of observations????
        warns <- c(warns,
          paste0("⚠️ CPO failures detected (", round(perc_failures, 2), "% of observations).\n",
          "   Model severely struggles to predict these points (CPO ≈ 0).\n",
          "   Recommended: check for outliers or review model specification/priors."))
      }
    }
  }

  # plots
  # pA Hyperparameters: posterior marginals + PC priors
  post_df <- NA_real_
  if(!is.null(fit$marginals.hyperpar)) {
    for(nm in names(fit$marginals.hyperpar)) {
      sm <- INLA::inla.smarginal(fit$marginals.hyperpar[[nm]])
      post_df <- rbind(post_df, data.frame(x = sm$x, y = sm$y, par = nm))
      post_df <- na.omit(post_df)
    }
  }
  # Create prior reference ticks for vertical lines
  prior_ticks <- do.call(rbind, Filter(Negate(is.null), list(
    if(!is.null(priors$local.pcprior.range)) data.frame(x = priors$local.pcprior.range[1], par = "Range for Sloc"),
    if(!is.null(priors$local.pcprior.sigma)) data.frame(x = priors$local.pcprior.sigma[1], par = "Stdev for Sloc"),
    if(!is.null(priors$shared.pcprior.range)) data.frame(x = priors$shared.pcprior.range[1], par = "Range for Sshared"),
    if(!is.null(priors$shared.pcprior.sigma)) data.frame(x = priors$shared.pcprior.sigma[1], par = "Stdev for Sshared")
  )))
  if(!is.null(prior_ticks) && nrow(prior_ticks) > 0 && !is.null(post_df)) {
    ymax_fac <- aggregate(y ~ par, post_df, function(z) max(z, na.rm = TRUE))
    prior_ticks$y <- ymax_fac$y[match(prior_ticks$par, ymax_fac$par)] * 0.95
  }

  pA <- ggplot2::ggplot(post_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(linewidth = 0.6, color = "#1a5276") +
    ggplot2::facet_wrap(~par, scales = "free", ncol = 2) +
    ggplot2::geom_vline(data = prior_ticks, ggplot2::aes(xintercept = x),
                        linetype = "dashed", linewidth = 0.5, color = "#c0392b", alpha = 0.7) +
    ggplot2::geom_text(data = prior_ticks,
                       ggplot2::aes(x = x, y = y, 
                       label = paste0("PC prior (u) = ", round(x, 2))),
                       vjust = -0.4, hjust = 1, size = 2.8,
                       color = "#c0392b", angle = 90) +
    ggplot2::labs(
      title = "A) Hyperparameters: posterior marginals",
      subtitle = "Red dashed lines = PC prior (u)",
      y = "Density", x = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
      strip.text = ggplot2::element_text(face = "bold", size = 9, color = "#2c3e50"),
      axis.title = ggplot2::element_text(size = 9, color = "#2c3e50"),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#d5d8dc")
    )

  # pB Sloc vs Sshared fields
  ras_to_df <- function(r, nm) {
    rr <- terra::unwrap(r)[["mean"]]
    df <- terra::as.data.frame(rr, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", "mean")
    df$which <- nm
    df
  }
  maps_df <- data.frame()
  if(!is.null(pred_local)) maps_df <- rbind(maps_df, ras_to_df(pred_local, "Sloc field"))
  if(!is.null(pred_shared)) maps_df <- rbind(maps_df, ras_to_df(pred_shared, "Sshared field"))

  ordered_hierarchical <- any(grepl("copy", rownames(fit$summary.hyperpar), ignore.case = TRUE)) ||
                   "IGlobal" %in% names(fit$summary.random)
  # Define midpoint dynamically
  has_copy_structure <- any(grepl("copy", rownames(fit$summary.hyperpar), ignore.case = TRUE)) ||
                   "IGlobal" %in% names(fit$summary.random)
  if(has_copy_structure) {
    midpoint_val <- mean(maps_df$mean, na.rm = TRUE)
  } else {
    midpoint_val <- 0
  }

  if(nrow(maps_df) > 0) {
    zlim <- range(maps_df$mean, na.rm = TRUE)
    pB <- ggplot2::ggplot(maps_df, ggplot2::aes(x = x, y = y, fill = mean)) +
      ggplot2::geom_raster(na.rm = TRUE) +
      #ggplot2::scale_fill_distiller(palette = "YlGnBu", limits = zlim, na.value = "white") +
      ggplot2::scale_fill_gradient2(
        low = "#c0392b", mid = "white", high = "#1a5276",
        midpoint = midpoint_val,
        limits = zlim, na.value = "white",
        oob = scales::squish
      ) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::facet_wrap(~which, ncol = 2, scales = "fixed") +
      ggplot2::labs(
        title = "B) Sloc fields (posterior mean)",
        subtitle = if(has_copy_structure)
          "Red = below global mean, Blue = above global mean (centered at model mean)"
        else
          "Red = below average, Blue = above average (centered at zero)",
        fill = "Posterior\nmean"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
        plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
        strip.text = ggplot2::element_text(face = "bold", size = 9, color = "#2c3e50"),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.key.height = ggplot2::unit(0.3, "cm"),
        legend.key.width = ggplot2::unit(1.2, "cm"),
        legend.title = ggplot2::element_text(size = 9, color = "#2c3e50"),
        legend.text = ggplot2::element_text(size = 8)
      )
  } else {
    pB <- ggplot2::ggplot() +
      ggplot2::labs(
        title = "B) Sloc fields (posterior mean)",
        subtitle = "No Sloc or Sshared fields present in the model"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
        plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e")
      )
  }

  if (!is.null(rs)) {
  # pC Residual correlogram
  cor_df <- data.frame()
  coords <- as.matrix(data_used[, c("x", "y")])
  if(is.finite(moran_I)) {
    dmat <- as.matrix(stats::dist(coords))
    dvec <- dmat[upper.tri(dmat)]
    nbins <- 15L
    brks <- seq(0, max(dvec, na.rm = TRUE), length.out = nbins + 1)
    mid <- 0.5 * (brks[-1] + brks[-length(brks)])
    rho <- rep(NA_real_, nbins)
    ut_r <- row(dmat)[upper.tri(dmat)]
    ut_c <- col(dmat)[upper.tri(dmat)]
    for(i in seq_len(nbins)) {
      sel <- dvec >= brks[i] & dvec < brks[i + 1]
      if(sum(sel) > 30) {
        idx <- which(sel)
        r1 <- rs[ut_r[idx]]; r2 <- rs[ut_c[idx]]
        rho[i] <- stats::cor(r1, r2, use = "pairwise.complete.obs")
      }
    }
    cor_df <- data.frame(dist_mid = mid, rho = rho)
  }
  range_eff <- if(is.finite(range_res_mean)) range_res_mean else if(is.finite(range_lat_mean)) range_lat_mean else NA_real_
  pC <- ggplot2::ggplot(cor_df, ggplot2::aes(x = dist_mid, y = rho)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(na.rm = TRUE, size = 1.2, color = "#1a5276") +
    ggplot2::geom_line(na.rm = TRUE, color = "#1a5276", linewidth = 0.6) +
    ggplot2::labs(
      title = "C) Residual correlogram (residuals: obs − fitted mean)",
      subtitle = if(is.finite(range_eff))
        paste0("Model range ≈ ", round(range_eff, 3), " (map units)\n(distance where correlation vanishes)")
      else
        "No Sloc/Sshared field: full extent shown",
      x = "Distance (map units)",
      y = "Residual correlation"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
      strip.text = ggplot2::element_text(face = "bold", size = 9, color = "#2c3e50"),
      axis.title.x = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, color = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#d5d8dc"),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  # Vertical line for model range
  if(is.finite(range_eff)) {
    pC <- pC +
      ggplot2::geom_vline(xintercept = range_eff, linetype = "dashed", color = "#c0392b", alpha = 0.7) +
      ggplot2::annotate(
        "text",
        x = range_eff, y = max(cor_df$rho, na.rm = TRUE),
        label = "Model range", angle = 90, vjust = -0.8, hjust = 0.9,
        color = "#c0392b", size = 3
      )
  }

  # pD Residual histogram + QQ plot
  res_mean <- mean(rs, na.rm = TRUE)
  df_r <- data.frame(resid = rs)
  center_label <- if(ordered_hierarchical) {
    paste0("Mean residual (≈ ", round(res_mean, 3), ")")
  } else {
    "Model mean = 0"
  }
  line_x <- if(ordered_hierarchical) res_mean else 0

  pD1 <- ggplot2::ggplot(df_r, ggplot2::aes(x = resid)) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = "#1a5276",
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#c0392b", linewidth = 0.4) +
    ggplot2::annotate("text", x = line_x, y = Inf,
    label = center_label, angle = 90, vjust = -0.8,
    hjust = 1.2, size = 2.8, color = "#c0392b"
    ) +
    ggplot2::labs(
      title = "D) Residual histogram & QQ-plot (residuals: obs − fitted mean)",
      subtitle = paste0("Distribution of residuals\n(dashed line = ", center_label, ")"),
      x = "Residuals",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, color = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#d5d8dc"),
      plot.background = ggplot2::element_blank()
    )
  qq <- qqnorm(rs, plot.it = FALSE)
  df_qq <- data.frame(theoretical = qq$x, sample = qq$y)
  pD2 <- ggplot2::ggplot(df_qq, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point(color = "#1a5276", size = 1.3, alpha = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#c0392b", linewidth = 0.4) +
    ggplot2::labs(
      title = " ",
      subtitle = "Residuals vs. theoretical quantiles\n(dashed = normal expectation)",
      x = "Theoretical quantiles (Normal)",
      y = "Sample quantiles (residuals)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, color = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#d5d8dc"),
      plot.background = ggplot2::element_blank()
    )

  # pE semivariaogram
  coords <- as.matrix(data_used[, c("x", "y")])
  if(has_Sloc) {
    sp_vals <- fit$summary.random$Sloc$mean
    sp_vals <- sp_vals[seq_len(nrow(data_used))]
  } else {
    sp_vals <- 0
  }
  if(has_Sshared) {
    lat_vals <- fit$summary.random$Sshared$mean
    lat_vals <- lat_vals[seq_len(nrow(data_used))]
  } else {
    lat_vals <- 0
  }
  rs_net <- y_obs - y_pred - sp_vals - lat_vals
  dmat <- as.matrix(stats::dist(coords))
  dvec <- dmat[upper.tri(dmat)]
  r1 <- rs_net[row(dmat)[upper.tri(dmat)]]
  r2 <- rs_net[col(dmat)[upper.tri(dmat)]]
  delta <- r1 - r2
  nbins <- 15L
  brks <- seq(0, max(dvec, na.rm = TRUE), length.out = nbins + 1)
  mid  <- 0.5 * (brks[-1] + brks[-length(brks)])
  gamma <- numeric(nbins)
  for(i in seq_len(nbins)) {
    sel <- dvec >= brks[i] & dvec < brks[i+1]
    if(sum(sel) > 20) {
      gamma[i] <- var(delta[sel], na.rm = TRUE) / 2
    } else {
      gamma[i] <- NA_real_
    }
  }
  sv_df <- data.frame(
    dist = mid,
    gamma = gamma
  )

  pE <- ggplot2::ggplot(sv_df, ggplot2::aes(x = dist, y = gamma)) +
    ggplot2::geom_point(color = "#1a5276", size = 1.5, alpha = 0.8) +
    ggplot2::geom_line(color = "#1a5276", linewidth = 0.6, alpha = 0.8) +
    ggplot2::labs(
      title = "E) Empirical semivariogram (residuals: obs − fitted mean)",
      subtitle = paste0(" "),
      x = "Distance (map units)",
      y = "Semivariance γ(h)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, color = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, color = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#d5d8dc")
    )
  if(has_Sloc) {
    pE <- pE + ggplot2::geom_vline(xintercept = range_res_mean, linetype = "dashed", color = "#c0392b", linewidth = 0.5) +
      ggplot2::geom_text(
        data = data.frame(x = range_res_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Sloc range"),
        ggplot2::aes(x = x, y = y, label = label), 
        color = "#c0392b", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }
  if(has_Sshared) {
    pE <- pE + ggplot2::geom_vline(xintercept = range_lat_mean, linetype = "dashed", color = "#2980b9", linewidth = 0.5) +
      ggplot2::geom_text(
      data = data.frame(x = range_lat_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Sshared range"),
      ggplot2::aes(x = x, y = y, label = label),
      color = "#2980b9", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }

  } else {
    # if family cp
    pC <- NULL
    pD1 <- NULL
    pD2 <- NULL
    pE <- NULL
  }  

  # pH covariates importance
  #pH <- .nsbm_vars_importance(fit)


  # 
  out <- list(
    bayes_fit = list(dic_val = dic_val,
                      waic_val = waic_val,
                      mlik_val = mlik_val, 
                      mlpd_val = mlpd_val,
                      lcpo_val = lcpo_val),
    hiperpars = list(range_res = range_res,
                      sigma_res = sigma_res,
                      range_lat = range_lat,
                      sigma_lat = sigma_lat,
                      prec_IGlobal = prec_IGlobal),
    intercepts = list(iGlobal = iGlobal,
                       iRegional = iRegional,
                       copy_beta = beta_IRegional,
                       copy_beta_predictors = beta_pred_oh),
    fixed_covariates = sig_vars,
    predictive = list(auc_full = auc_full,
                       tjur_r2 = tjur_r2,
                       brier = brier,
                       rmse = rmse, 
                       corr_obs_pred = pred_cor),
    calibration = list(slope = cal_slope,
                       ks_pit = ks_pit,
                       pit_values = if(exists("pit_reg")) pit_reg else NULL),
    coverage = list(cov50 = cov50,
                    cov95 = cov95),
    diagnostics = list(moran_I = moran_I,
                       max_CIratio = max_CIratio, # relative uncertainty
                       range_ratio = range_ratio, # scale separation ratio
                       sigma_ratio = sigma_ratio, # variance ratio
                       field_correlation = field_correlation), # prior influence
    warnings = warns,
    plots = list(hyperparams = pA,
                 Slocfields = pB,
                 correlogram = pC,
                 hist = pD1,
                 qq = pD2,
                 #vars_importance = pF,
                 semivariogram = pE)
  )

  return(out)
}


# -----------------------------


#' prepare significant covariates
#' @noRd
.signif_vars <- function(fit, scale_params = NULL) {
  sf <- fit$summary.fixed

  # back-transform coef to original scale if vars were standardized
  if(!is.null(scale_params) && nrow(sf) > 0) {
    for(i in seq_len(nrow(sf))) {
      coef_name <- rownames(sf)[i]
      # extract base variable name (remove GL, RE, RE_oh, GL_ls, RE_ss suffixes)
      base_var <- gsub("GL$|RE$|RE_oh$|GL_ls$|RE_ss$", "", coef_name)
      if(base_var %in% names(scale_params) && scale_params[[base_var]]$sd > 0) {
        sd_x <- scale_params[[base_var]]$sd
        sf[i, "mean"]       <- sf[i, "mean"]       / sd_x
        sf[i, "sd"]         <- sf[i, "sd"]          / sd_x
        sf[i, "0.025quant"] <- sf[i, "0.025quant"]  / sd_x
        sf[i, "0.975quant"] <- sf[i, "0.975quant"]  / sd_x
      }
    }
  }
  
  # p-value aprox with marginal posterior
  tailprob <- function(term){
    if (!is.null(fit$marginals.fixed) && !is.null(fit$marginals.fixed[[term]])) {
      p0 <- INLA::inla.pmarginal(0, fit$marginals.fixed[[term]])
      2 * min(p0, 1 - p0)
    } else {
      # fallback Normal
      2 * stats::pnorm(-abs(sf[term,"mean"] / sf[term,"sd"]))
    }
  }

  # 95% CI no cruza 0
  signif95 <- sf[,"0.025quant"] * sf[,"0.975quant"] > 0
  
  stars <- ifelse(signif95, "***", "")
  
  out <- data.frame(
    coef = rownames(sf),
    estimate = sf$mean,
    sd = sf$sd,
    `2.5%` = sf[,"0.025quant"],
    `97.5%` = sf[,"0.975quant"],
    signif = stars,
    tail_p = vapply(rownames(sf), tailprob, numeric(1)),
    row.names = NULL,
    check.names = FALSE
  )
  
  out <- out[order(-abs(out$estimate)), ]

  out$estimate <- round(out$estimate, 6)
  out$sd <- round(out$sd, 6)
  out$`2.5%` <- round(out$`2.5%`, 6)
  out$`97.5%` <- round(out$`97.5%`, 6)
  
  return(out)
}


# -----------------------------


#' prepare summary
#' @noRd
.nsbm_generate_summary <- function(fit, species, fam, lnk, coupling.intercept, coupling.predictors, diag_block, cv_res=NULL, vg=NULL, vr=NULL, scale_params=NULL, has_spatial=FALSE) {
  
  fmt_val <- function(x, digits = 3) {
    if(is.null(x) || length(x) == 0) return("—")
    if(is.character(x)) return(x)
    if(!is.finite(x)) return("—")
    round(x, digits)
  }

  fmt_pm <- function(mean, sd, digits = 3) {
    if(!is.finite(mean) || !is.finite(sd)) return(NULL)
    paste0(round(mean, digits), " ± ", round(sd, digits))
  }

  has_Sloc <- "Sloc" %in% names(fit$summary.random)
  has_Sshared <- "Sshared" %in% names(fit$summary.random)
  has_covariates <- nrow(fit$summary.fixed) > 0

  # metadata
  species_name <- gsub("\\.", " ", species)
  base_name <- if (!has_spatial) "Non-spatial baseline" else "NSBM"
  suf <- c()
  if(has_Sloc) suf <- c(suf, "Sloc")
  if(has_Sshared) suf <- c(suf, "Sshared")
  if(has_covariates) suf <- c(suf, "covariates")

  model_type <- if(length(suf) > 0) {
    paste0(base_name, ": ", paste(suf, collapse = " + "))
  } else {
    base_name
  }
  #model_type <- paste0("NSBM",
  # if(has_Sloc) " + Sloc" else "",
  # if(has_Sshared) " + Sshared" else "",
  # if(has_covariates) " + covariates" else ""
  #)

  cp_int <- if(is.null(coupling.intercept)) "NULL" else coupling.intercept
  cp_pred <- if(is.list(coupling.predictors)) "custom list" else if(is.null(coupling.predictors)) "NULL" else coupling.predictors
  
  tbl_metadata <- data.frame(
    Field = c("Species name:", "Model type:", "Family | Link:", 
              "Coupling (Intercept):", "Coupling (Predictors):"),
    Value = c(species_name,
              model_type,
              paste0(fam, " | ", ifelse(is.null(lnk), "—", lnk)),
              cp_int,
              cp_pred),
    stringsAsFactors = FALSE
  )

  # bayesian model fit
  tbl_fit <- data.frame(
    Metric = c(
      "DIC",
      "WAIC",
      "Marginal log-likelihood (log ML)",
      "LCPO (sum log-CPO)",
      "MLPD (mean log predictive density)"
    ),
    Value = c(diag_block$bayes_fit$dic_val, 
              diag_block$bayes_fit$waic_val,
              diag_block$bayes_fit$mlik_val, 
              diag_block$bayes_fit$lcpo_val,
              diag_block$bayes_fit$mlpd_val),
    stringsAsFactors = FALSE
  )

  # hyperparametres
  hyp_block <- diag_block$hiperpars
  params <- character(0)
  values <- character(0)

  if(has_Sloc) {
    params <- c(
      params,
      "Sloc local field – Range (posterior mean ± SD)",
      "Sloc local field – Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_res,
      hyp_block$sigma_res
    )
  }
  if(has_Sshared) {
    params <- c(
      params,
      "Sshared global field – Range (posterior mean ± SD)",
      "Sshared global field – Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_lat,
      hyp_block$sigma_lat
    )
  }
  params <- c(params, "Precision for IGlobal (mean ± SD)")
  values <- c(values, hyp_block$prec_IGlobal)

  tbl_hyper <- data.frame(
    Parameter = params,
    Value = values,
    stringsAsFactors = FALSE
  )

  # intercepts
  int_block <- diag_block$intercepts
  rows <- list()
  # IGlobal exists only for coupling.intercept = "additive"/"hierarchical"
  if(!is.null(int_block$iGlobal)) {
    rows[["IGlobal (mean ± SD, CI95%)"]] <- int_block$iGlobal
  }
  # IRegional always exists
  if(!is.null(int_block$iRegional)) {
    if(!is.null(coupling.intercept) && coupling.intercept == "ordered_hierarchical") {
      rows[["IRegional (copy from IGlobal, mean ± SD, CI95%)"]] <- int_block$iRegional
    } else {
      rows[["IRegional (mean ± SD, CI95%)"]] <- int_block$iRegional
    }
  }
  # copy_beta intercept: only exists when coupling.intercept = "ordered_hierarchical"
  if(!is.null(int_block$copy_beta) && int_block$copy_beta != "—") {
    rows[["Copy β (IGlobal -> IRegional) (mean ± SD)"]] <- int_block$copy_beta
  }
  # copy_beta predictors: one per variable with ordered_hierarchical coupling
  if(!is.null(int_block$copy_beta_predictors)) {
    for(nm in names(int_block$copy_beta_predictors)) {
      var_name <- sub("^Beta for (.+)RE_oh$", "\\1", nm)
      rows[[paste0("Copy β (", var_name, "GL -> ", var_name, "RE_oh) (mean ± SD)")]] <- int_block$copy_beta_predictors[[nm]]
    }
  }

  tbl_intercepts <- data.frame(
    Term = names(rows),
    Value = unname(unlist(rows)),
    stringsAsFactors = FALSE
  )

  # fixed covariates
  tbl_fixed <- diag_block$fixed_covariates
  if(is.null(tbl_fixed) || !nrow(tbl_fixed)) {
    tbl_fixed <- data.frame(Term = "No significant covariates detected")
  } else {
    # var name
    base_vars <- gsub("GL$|RE$", "", tbl_fixed$coef)
    
    # which coupling?
    tbl_fixed$Coupling <- vapply(base_vars, function(v) {
      if(v %in% c(vg, vr)) {
        .resolve_coupling_predictor(v, coupling.predictors, vg)
      } else {
        "—"
      }
    }, character(1))
    
    cols_order <- c("coef", "coupling", "estimate", "sd", "2.5%", "97.5%", "signif", "tail_p")
    cols_order <- intersect(cols_order, names(tbl_fixed))
    tbl_fixed <- tbl_fixed[, cols_order]
  }

  # predictive performance
  tbl_pred <- data.frame(
    Metric = c("AUC (full model)", 
               "Brier score", 
               "RMSE", 
               "Observed-predicted correlation (r)"),
    Value  = c(diag_block$predictive$auc_full,
               fmt_val(diag_block$predictive$brier),
               fmt_val(diag_block$predictive$rmse),
               fmt_val(diag_block$predictive$corr_obs_pred)),
    stringsAsFactors = FALSE
  )

  # calibration & coverage
  if (fam %in% c("cp", "binomial")) {
    tbl_cal <- NULL
  } else {
    cal_block <- diag_block$calibration
    cov_block <- diag_block$coverage

    tbl_cal <- data.frame(
      Metric = c("Calibration slope",
                 "PIT KS p-value",
                 "Coverage (central 50%)",
                 "Coverage (central 95%)"),
      Value = c(fmt_val(cal_block$slope),
                fmt_val(cal_block$ks_pit),
                fmt_val(cov_block$cov50),
                fmt_val(cov_block$cov95)),
      stringsAsFactors = FALSE
    )
  }

  # diagnostics
  tbl_diag <- data.frame(
    Metric = c(
      "Residual Moran's I",
      "Max CI/median ratio (hyperparameters)",
      "Scale-separation ratio (range_Sshared / range_Sloc)",
      "Variance ratio (sigma_Sshared / sigma_Sloc)",
      "Sshared–Sloc field correlation (r)",
      "Variance explained by Sloc (%)",
      "Scale separation Index (SSI) [0–1]",
      "Spatial redundancy Index (SRI) [0–1]"
    ),
    Value = c(
      fmt_val(diag_block$diagnostics$moran_I),
      fmt_val(diag_block$diagnostics$max_CIratio),
      fmt_val(diag_block$diagnostics$range_ratio),
      fmt_val(diag_block$diagnostics$sigma_ratio),
      fmt_val(diag_block$diagnostics$field_correlation),
      fmt_val(var_explained_sloc * 100),
      fmt_val(ssi),
      fmt_val(sri) 
    ),
    stringsAsFactors = FALSE
  )


  # cv
  tbl_cv <- NULL
  if(!is.null(cv_res)) {
    metric_label <- paste0("CV ", cv_res$metric_name, " (mean ± sd)")
    tbl_cv <- data.frame(
      Metric = c("CV folds", metric_label),
      Value = c(cv_res$cv.folds,
                sprintf("%.3f ± %.3f", cv_res$metric_mean, cv_res$metric_sd)),
      stringsAsFactors = FALSE
    )
  }


  drop_null_rows <- function(df) {
    ok <- !sapply(df$Value, function(x) is.null(x) || identical(x, "—"))
    df[ok, , drop = FALSE]
  }

  tbl_fit <- drop_null_rows(tbl_fit)
  tbl_hyper <- drop_null_rows(tbl_hyper)
  tbl_intercepts <- drop_null_rows(tbl_intercepts)
  tbl_pred <- drop_null_rows(tbl_pred)
  if (!is.null(tbl_cal)) tbl_cal <- drop_null_rows(tbl_cal)
  tbl_diag <- drop_null_rows(tbl_diag)

  # list of tables
  out <- list(
    Metadata = tbl_metadata,
    `Model fit` = tbl_fit,
    Hyperparameters = tbl_hyper,
    Intercepts = tbl_intercepts,
    `Fixed effects` = tbl_fixed,
    `Predictive performance` = tbl_pred,
    Diagnostics = tbl_diag
  )
  if(!is.null(tbl_cal) && nrow(tbl_cal) > 0) out[["Calibration & coverage"]] <- tbl_cal
  if(!is.null(tbl_cv)) out[["Cross-validation"]] <- tbl_cv

  return(out)
}


# -----------------------------


#' plot covariates importance
#' @noRd
.nsbm_vars_importance <- function(fit) {

  df <- fit$summary.fixed
  df$var <- rownames(df)

  # exclude intercepts and Sshared
  drop <- c("IGlobal", "IRegional", "beta_GL")
  df <- df[!(df$var %in% drop), , drop = FALSE]

  df$sign <- ifelse(df$mean > 0, "Positive", "Negative")

  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = reorder(var, abs(mean)),
                 y = abs(mean),
                 fill = sign)) +    ggplot2::geom_col(alpha = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Positive" = "#1a5276","Negative" = "#c0392b")) +
    ggplot2::labs(
      title = "F) Variable importance (|β|)",
      subtitle = "Higher bars indicate stronger effect magnitude",
      x = "Covariate",
      y = "|Effect size|"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold",size = 11,color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9,color = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9,margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9,margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8,color = "#2c3e50"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2,color = "#d5d8dc")
    )


  return(p)
}


###----------------###


