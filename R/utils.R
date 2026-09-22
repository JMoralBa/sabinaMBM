
#' @importFrom stats aggregate as.formula binomial coef dist ks.test lm predict qlogis qqnorm reorder sd var
#' @importFrom utils write.csv
#' @importFrom sf st_sfc
NULL
# global vars check
utils::globalVariables(c(
  "x", "y", "dist_mid", "rho", "resid", "intercept", 
  "significant", "density", "theoretical", "label", "type"
))


# -----------------------------


#' Logs
#' @noRd
.info <- function(msg, verbose = TRUE) if(verbose) message("\u2139 ", msg)
.check <- function(msg, verbose = TRUE) if(verbose) message("  ✓ ", msg)
.item <- function(msg, verbose = TRUE) if(verbose) message("    - ", msg)
.warn <- function(msg, verbose = TRUE) warning("⚠️  ", msg, call. = FALSE)
.stop <- function(msg, verbose = TRUE) stop("❌ ", msg, call. = FALSE)


# -----------------------------


#' Fast standardization of SpatRaster stack
#' @noRd
.standardize_rasters <- function(sp_rast, var_names = NULL, n_cores = 1) {

  if(is.null(sp_rast)) return(list(rast = NULL, params = list()))
  if(!inherits(sp_rast, "SpatRaster")) return(list(rast = sp_rast, params = list()))

  if(is.null(var_names)) var_names <- names(sp_rast)
  if(length(var_names) == 0) return(list(rast = sp_rast, params = list()))

  stats_all <- terra::global(sp_rast[[var_names]], fun = c("mean", "sd"), na.rm = TRUE)
  
  stats_list <- list()
  for(v in var_names) {
    m_val <- stats_all[v, "mean"]
    s_val <- stats_all[v, "sd"]
    
    stats_list[[v]] <- list(mean = m_val, sd = s_val)
    
    # standarization
    if(!is.na(s_val) && s_val > 1e-10) {
      sp_rast[[v]] <- (sp_rast[[v]] - m_val) / s_val
    } else {
      .warn(paste0("Variable '", v, "' has near-zero variance. Centering only."))
      sp_rast[[v]] <- sp_rast[[v]] - m_val
    }
  }

  list(rast = sp_rast, params = stats_list)
}


# -----------------------------


#' Multiscale Z-standardization with COS correction
#' @noRd
.standardize_multiscale <- function(sp_covglo, sp_covreg, unified_vars = character(0),
                                    standardize_global = TRUE) {

  scale_params_glo <- list()
  scale_params_reg <- list()

  vars_glo <- if (!is.null(sp_covglo)) names(sp_covglo) else character(0)
  vars_reg <- if (!is.null(sp_covreg)) names(sp_covreg) else character(0)

  # unified variables (shared + coupled with scale_decomposed or bayesian_feedback)
  # Use global statistics for both rasters. X_GL and X_RE in the same Z-space, eliminating COS artifacts in cross-scale coefficient transfer.
  for (v in unified_vars) {
    if (!(v %in% vars_glo)) {
      .warn(paste0("Unified standardization requested for '", v,
                   "' but variable is missing in global raster. Falling back to regional stats."))
      next
    }
    m_glo <- terra::global(sp_covglo[[v]], "mean", na.rm = TRUE)[1, 1]
    s_glo <- terra::global(sp_covglo[[v]], "sd", na.rm = TRUE)[1, 1]


    if (is.na(s_glo) || s_glo <= 1e-10) {
      .warn(paste0("Variable '", v, "' has near-zero global variance. Centering only."))
      if (standardize_global) sp_covglo[[v]] <- sp_covglo[[v]] - m_glo
      if (v %in% vars_reg) sp_covreg[[v]] <- sp_covreg[[v]] - m_glo
      scale_params_glo[[v]] <- list(mean = m_glo, sd = NA_real_)
      scale_params_reg[[v]] <- list(mean = m_glo, sd = NA_real_, unified = TRUE)
    } else {
      if (standardize_global) sp_covglo[[v]] <- (sp_covglo[[v]] - m_glo) / s_glo
      if (v %in% vars_reg)    sp_covreg[[v]] <- (sp_covreg[[v]] - m_glo) / s_glo
      scale_params_glo[[v]] <- list(mean = m_glo, sd = s_glo)
      scale_params_reg[[v]] <- list(mean = m_glo, sd = s_glo, unified = TRUE)
    }
  }

  # non-unified variables. Independent standardization
  if (standardize_global) {
    non_unif_glo <- setdiff(vars_glo, unified_vars)
    if (length(non_unif_glo) > 0) {
      res_glo <- .standardize_rasters(sp_covglo, var_names = non_unif_glo)
      if (!is.null(res_glo$rast)) {
        sp_covglo <- res_glo$rast
        for (v in names(res_glo$params)) {
          scale_params_glo[[v]] <- res_glo$params[[v]]
        }
      }
    }
  }

  non_unif_reg <- setdiff(vars_reg, unified_vars)
  if (length(non_unif_reg) > 0) {
    res_reg <- .standardize_rasters(sp_covreg, var_names = non_unif_reg)
    if (!is.null(res_reg$rast)) {
      sp_covreg <- res_reg$rast
      for (v in names(res_reg$params)) {
        scale_params_reg[[v]] <- list(mean = res_reg$params[[v]]$mean, 
                                      sd = res_reg$params[[v]]$sd, 
                                      unified = FALSE)
      }
    }
  }
  
  list(sp_covglo = sp_covglo,
       sp_covreg = sp_covreg,
       scale_params_glo = scale_params_glo,
       scale_params_reg = scale_params_reg)
}


# -----------------------------


#' Build inlabru likelihood objects
#' @noRd
.build_likelihoods <- function(fam, lnk, rhs_glo, rhs_reg,
                               pp_glo, pp_reg, bdy_glo, bdy_reg, dom,
                               coupling.intercept) {
  lik_glo <- NULL
  if(!is.null(coupling.intercept)) {
    lik_glo <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ ", rhs_glo)),
      data = pp_glo, samplers = bdy_glo, domain = dom,
      control.family = list(link = lnk))
  }
  lik_reg <- inlabru::like(
    family = fam,
    formula = as.formula(paste0("resp ~ ", rhs_reg)),
    data = pp_reg, samplers = bdy_reg, domain = dom,
    control.family = list(link = lnk))
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
                 coupling.covariates = NULL,
                 slope_delta_prior = NULL,
                 pp_glo_sf = NULL,
                 pp_reg_sf = NULL) {

  # selected vars
  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional

  # Build rw2 mesh using fmesher::fm_mesh_1d
  .build_rw2_values <- function(rast_layer, coords, n_knots = 100L) {
    vals <- suppressWarnings(as.numeric(terra::extract(rast_layer, coords)[, 1]))
    vals <- vals[is.finite(vals)]
    rng <- range(vals)
    if(diff(rng) == 0) {
      .stop("RW2 requires variability (covariate range = 0).\n   Use 'linear' instead for this covariate or check the raster values.")
    }
    fmesher::fm_mesh_1d(
      loc      = seq(rng[1], rng[2], length.out = n_knots),
      degree   = 2,
      boundary = c("neumann", "free")
    )
  }

 
  spec_glo <- list()
  if(length(vg) > 0) {
    spec_glo <- lapply(stats::setNames(vg, vg), function(x) .resolve_covariate_effects(x, "global", covariate.effects))
  }
  
  spec_reg <- list()
  if(length(vr) > 0) {
    spec_reg <- lapply(stats::setNames(vr, vr), function(x) .resolve_covariate_effects(x, "regional", covariate.effects))
  }

  # detect rw2 needs
  need_rw2_global <- if(length(spec_glo) > 0) any(vapply(spec_glo, function(s) s$model == "rw2", logical(1))) else FALSE
    need_rw2_regional <- if(length(spec_reg) > 0) any(vapply(spec_reg, function(s) s$model == "rw2", logical(1))) else FALSE

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

  for(X in vg) {
    cp_mode <- .resolve_coupling_predictor(X, coupling.covariates, vg, vr)
    if(cp_mode == "NULL") next
    specX <- spec_glo[[X]]
    if(specX$model == "drop") next
    if(specX$model == "linear") {
      if(cp_mode %in% c("nested_shrinkage", "ordered_hierarchical")) {
           cmpglobal <- c(cmpglobal,
                           paste0(X, "GL(main = rep(1L, nrow(.data.)), model = 'iid', weights = as.numeric(terra::extract(", spobjglo,
                              ", sf::st_transform(.data., terra::crs(", spobjglo,
                              ")), method = 'bilinear', ID = FALSE)[['", X, "']]), hyper = list(prec = list(initial = -10, fixed = TRUE)))"))
        } else {
          cmpglobal <- c(cmpglobal,
                         paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'linear')"))
      }
    } else if(specX$model == "rw2") {
    mesh_1d <- .build_rw2_values(sp_covglo[[X]], coords_all)
    cmpglobal <- c(cmpglobal,
                   paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'rw2', scale.model = TRUE, ",
                   "mapper = mesh_1d, ",
                   "hyper = list(prec = list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))))"))
    }
    fglobal <- c(fglobal, paste0(X, "GL"))
  }

  # regional components 
  cmpregional <- character(0)
  fregional <- character(0)

  for(X in vr) {
    is_shared <- X %in% vg
    cp_mode <- if(is_shared) .resolve_coupling_predictor(X, coupling.covariates, vg, vr) else "unpooled"

    specX <- spec_reg[[X]]
    if(specX$model == "drop") next

    ## scale_decomposed: global (XGL_reg) + regional anomaly (XRE_delta)
    if(cp_mode == "scale_decomposed") {
      if(specX$model == "linear") {
        cmpregional <- c(cmpregional,
          paste0(X, "GL_glo_res(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "_glo_res']]), model = 'linear')"),
          paste0(X, "RE_reg_anom(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "_reg_anom']]), model = 'linear')"))
      } else if(specX$model == "rw2") {
        mesh_1d <- .build_rw2_values(sp_covreg[[paste0(X, "_reg_anom")]], coords_r)
        cmpregional <- c(cmpregional,
          paste0(X, "GL_glo_res(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "_glo_res']]), model = 'linear')"),
          paste0(X, "RE_reg_anom(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "_reg_anom']]), model='rw2', scale.model=TRUE, mapper = mesh_1d, hyper = list(prec=list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))))"))
      }
      fregional <- c(fregional, paste0(X, "GL_glo_res"), paste0(X, "RE_reg_anom"))
    }

    ## nested_shrinkage: beta_RE = beta_GL (shared via copy= real, fixed=TRUE) + delta_RE, delta_RE ~ N(0, sigma_delta^2).
    else if(cp_mode == "nested_shrinkage") {
      cmpregional <- c(cmpregional,
        paste0(X, "GL_copy_reg(main = rep(1L, nrow(.data.)), model = 'iid', weights = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), copy = '", X, "GL', hyper = list(beta = list(fixed = TRUE, initial = 1)))"),
        paste0(X, "_delta(main = rep(1L, nrow(.data.)), model = 'iid', weights = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), ", slope_delta_prior, ")"))
      fregional <- c(fregional, paste0(X, "GL_copy_reg"), paste0(X, "_delta"))
    }

    ## ordered_hierarchical: beta_RE = beta_copy * beta_GL via copy= real, beta_copy free (Krainski et al. 2018; Knorr-Held & Best 2001).
    else if(cp_mode == "ordered_hierarchical") {
      cmpregional <- c(cmpregional,
        paste0(X, "RE_oh(main = rep(1L, nrow(.data.)), model = 'iid', weights = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), copy = '", X, "GL', hyper = list(beta = list(prior = 'normal', param = c(1, 4))))"))
      fregional <- c(fregional, paste0(X, "RE_oh"))
    }

    # bayesian_feedback: global priors
    else if(cp_mode == "bayesian_feedback") {
      cmpregional <- c(cmpregional,
        paste0(X, "RE(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), model = 'linear', mean.linear = bf_mean_", X, ", prec.linear = bf_prec_", X, ")"))
      fregional <- c(fregional, paste0(X, "RE"))
    }

    # unpooled/NULL: coefs independent
    else {
      if(specX$model == "linear") {
        cmpregional <- c(cmpregional, paste0(X, "RE(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), model = 'linear')"))
      } else if(specX$model == "rw2") {
        mesh_1d <- .build_rw2_values(sp_covreg[[X]], coords_r)
        cmpregional <- c(cmpregional,
          paste0(X, "RE(main = as.numeric(terra::extract(", spobjreg, ", .data., ID = FALSE)[['", X, "']]), model = 'rw2', scale.model = TRUE, ",
                 "mapper = mesh_1d, ",
                 "hyper = list(prec = list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))))"))
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


#' Bayesian-feedback intercept
#' Global and regional background ratios differ. Fix before share the global intercept.
#'  Only for background points. Real absences offset 0.
#' @noRd
.bf_intercept_offset <- function(resp_glo, resp_reg, uses_background) {
  if(!isTRUE(uses_background)) return(0)

  n1_glo <- sum(resp_glo == 1L); n0_glo <- sum(resp_glo == 0L)
  n1_reg <- sum(resp_reg == 1L); n0_reg <- sum(resp_reg == 0L)
  if(n1_glo == 0 || n0_glo == 0 || n1_reg == 0 || n0_reg == 0) return(0)

  log((n1_reg / n0_reg) / (n1_glo / n0_glo))
}


# -----------------------------


#' Fit MBM
#' @noRd
.fit_jmbm <- function(cmp, lik_list, coupling.intercept,
                      coupling.covariates, needs_feedback,
                      vr = NULL, n.threads = 1, seed = NULL,
                      int.strategy = "eb", bf_delta_int = 0,
                      has_Sshared = FALSE, spde.mesh = NULL,
                      verbose = TRUE) {
  bru_opts <- list(
    control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
    control.inla = list(int.strategy = int.strategy),
    control.mode = list(restart = TRUE),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
    num.threads = n.threads
  )
  
  if(needs_feedback) {
    env_cmp <- environment(cmp)
    
    assign("bf_mean_int", 0, envir = env_cmp)
    assign("bf_prec_int", 1, envir = env_cmp)
    if(!is.null(vr)) {
      for(v in vr) {
        assign(paste0("bf_mean_", v), 0, envir = env_cmp)
        assign(paste0("bf_prec_", v), 1, envir = env_cmp)
      }
    }

    # fit global model using only the global likelihood (Phase 1)
    .info("Sequential bayesian feedback", verbose = verbose)
    .check("Phase 1: Fitting global model to extract posteriors...", verbose = verbose)
    fit_glo <- do.call(inlabru::bru, c(list(components = cmp), list(lik_list[[1]]), list(options = bru_opts)))

    # extract moments and inject them into the current environment

    # safe precision with floor to avoid Inf or near-0 variance
    .safe_prec <- function(sd_val, floor_prec = 1e-4, ceil_prec = 1e6) {
      prec <- 1 / (sd_val^2)
      prec <- pmax(floor_prec, pmin(ceil_prec, prec))
      prec
    }

    # warn if global marginal is strongly skewed (moment matching unreliable)
    .check_skewness <- function(marginal, label) {
      if(is.null(marginal)) return(invisible(NULL))
      sm <- try(INLA::inla.smarginal(marginal), silent = TRUE)
      if(inherits(sm, "try-error")) return(invisible(NULL))
      m  <- INLA::inla.emarginal(function(x) x, marginal)
      m2 <- INLA::inla.emarginal(function(x) (x-m)^2, marginal)
      m3 <- INLA::inla.emarginal(function(x) (x-m)^3, marginal)
      skew <- m3 / (m2^(3/2))
      if(abs(skew) > 1) {
        .warn(paste0("Bayesian feedback: global posterior for '", label,
                     "' has skewness = ", round(skew, 2),
                     " — moment matching (mean/sd) may be a poor approximation.",
                     " Consider using `coupling.intercept = 'ordered_hierarchical'` instead."))
      }
    }

    if(!is.null(coupling.intercept) && coupling.intercept == "bayesian_feedback") {
      int_random <- fit_glo$summary.random$IGlobal
      if(is.null(int_random) || nrow(int_random) != 1L) {
        .stop("bayesian_feedback: could not extract IGlobal posterior (unexpected structure).")
      }
      bf_mean_int <- int_random$mean + bf_delta_int
      bf_sd_int <- int_random$sd
      .check_skewness(fit_glo$marginals.random$IGlobal[[1]], "IGlobal")
      assign("bf_mean_int", bf_mean_int, envir = env_cmp)
      assign("bf_prec_int", .safe_prec(bf_sd_int), envir = env_cmp)
      .item(sprintf("Intercept mean = %.3f, sd = %.3f", bf_mean_int, bf_sd_int), verbose = verbose)
    }

    if(!is.null(fit_glo$summary.fixed)) {
      gl_effs <- rownames(fit_glo$summary.fixed)
      gl_effs <- gl_effs[grepl("GL$", gl_effs)]
      for(eff in gl_effs) {
        base_var <- sub("GL$", "", eff)
        bf_sd_v <- fit_glo$summary.fixed[eff, "sd"]
        .check_skewness(fit_glo$marginals.fixed[[eff]], eff)
        assign(paste0("bf_mean_", base_var), fit_glo$summary.fixed[eff, "mean"], envir = env_cmp)
        assign(paste0("bf_prec_", base_var), .safe_prec(bf_sd_v),                envir = env_cmp)
        .item(sprintf("%s mean = %.3f, sd = %.3f", eff, fit_glo$summary.fixed[eff, "mean"], bf_sd_v), verbose = verbose)
      }
    }

    # transfer prams Sshared
    cmp_stage2 <- cmp
    if(needs_feedback && isTRUE(has_Sshared)) {
      bf_hyper_sshared <- fit_glo$summary.hyperpar
      bf_mean_range <- bf_hyper_sshared["Range for Sshared", "mean"]
      bf_mean_sigma <- bf_hyper_sshared["Stdev for Sshared", "mean"]

      matern_shared_bf <- INLA::inla.spde2.pcmatern(
        mesh = spde.mesh,
        prior.range = c(bf_mean_range, 0.5),
        prior.sigma = c(bf_mean_sigma, 0.5)
      )
      assign("matern_shared_bf", matern_shared_bf, envir = env_cmp)

      cmp_str <- paste(deparse(cmp), collapse = "")
      cmp_str2 <- sub(
        'Sshared\\(main *= *geometry, *model *= *matern_shared\\)',
        "Sshared(main = geometry, model = matern_shared_bf)",
        cmp_str
      )
      cmp_stage2 <- stats::as.formula(cmp_str2, env = env_cmp)

      .item(sprintf("Sshared range = %.3f, sigma = %.3f", bf_mean_range, bf_mean_sigma),
            verbose = verbose)
    }

    # fit using only regional likelihood (Phase 2)
    .check("Phase 2: Fitting regional model using Phase 1 posteriors as priors...", verbose = verbose)
    fit <- do.call(inlabru::bru, c(list(components = cmp_stage2), list(lik_list[[length(lik_list)]]), list(options = bru_opts)))
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
      # defaults (u=5, alpha=0.01)
      return(list(model = "rw2", u = 0.5, alpha = 0.01))
    }
    .stop(paste0("Invalid keyword '", spec, "' in `", path_label, "`.\n",
                 "   Valid options: 'linear', 'drop', or list(model='rw2', u=..., alpha=...)."))
  }

  # if spec is a list, only rw2 admited 
  if(is.list(spec)) {
    if(is.null(spec$model)) {
      .stop(paste0("Invalid list specification in `", path_label, "`.\n",
                   "   Lists are only allowed for RW2 effects. Use: list(model='rw2', u=..., alpha=...).\n",
                   "   For linear effects use 'linear'; to exclude use 'drop'."))
    }
    if(identical(spec$model, "rw2")) {
      u_val  <- if(is.null(spec$u)) 5 else spec$u
      alpha_val <- if(is.null(spec$alpha)) 0.01 else spec$alpha
      return(list(model = "rw2", u = u_val, alpha = alpha_val))
    }
    .stop(paste0("Invalid `model` in `", path_label, "`.\n",
                 "   When using a list, only model='rw2' is permitted.\n",
                 "   For linear effects use 'linear'; to exclude use 'drop'."))
  }

  # unsupported type
  .stop(paste0("Unsupported type in `", path_label, "`.\n",
               "   Must be string ('linear'|'drop') or list(model='rw2', u=..., alpha=...)."))
}


# -----------------------------


#' interpret coupling.covariates
#' @noRd
.resolve_coupling_predictor <- function(var, coupling.covariates, vg = NULL, vr = NULL) {

  if(is.null(coupling.covariates)) return("unpooled")

  mode_val <- "unpooled"
  is_explicit <- FALSE

  if(is.character(coupling.covariates)) {
    mode_val <- coupling.covariates
  } else if(is.list(coupling.covariates)) {
    if(!is.null(coupling.covariates$variables) && !is.null(coupling.covariates$variables[[var]])) {
      mode_val <- coupling.covariates$variables[[var]]
      is_explicit <- TRUE
    } else if(!is.null(coupling.covariates$default)) {
      mode_val <- coupling.covariates$default
    }
  }
  
  # unpooled for vars only in regional (no global counterpart to couple with)
  if(!is.null(vg) && !(var %in% vg) && !is_explicit) {
    if(mode_val %in% c("ordered_hierarchical", "scale_decomposed", "bayesian_feedback", "nested_shrinkage")) {
      mode_val <- "unpooled"
    }
  }

  # unpooled for vars only in global (no regional counterpart to couple with)
  if(!is.null(vr) && !(var %in% vr) && !is_explicit) {
    if(mode_val %in% c("ordered_hierarchical", "scale_decomposed", "bayesian_feedback", "nested_shrinkage")) {
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
  
  vars_avail <- intersect(vars_to_export, names(pred))
  if(length(vars_avail) == 0) return(NULL)
  
  if(sf::st_crs(pred) != terra::crs(template)) {
    pred <- sf::st_transform(pred, terra::crs(template))
  }
  
  vect_pred <- terra::vect(pred)
  r_pred <- terra::rasterize(vect_pred, template, field = vars_avail)
  
  rm(vect_pred, pred)
  gc(verbose = FALSE)
  
  names(r_pred) <- vars_avail
  return(r_pred)
}


.pred_as_tif2 <- function(pred, template, vars_to_export = c("mean", 
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
      sv_pts <- terra::vect(df, geom = c("X", "Y"), crs = terra::crs(sf::st_crs(pred)$wkt))
      sv_pts <- terra::project(sv_pts, terra::crs(template))
      r <- terra::rasterize(sv_pts, template, field = "value")
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
  # K= nomber of folds
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  f <- integer(length(y))
  f[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  f[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  f
}


# -----------------------------
#' format mean, sd and 95% CI strings
#' @noRd
.mean_sd_str <- function(df) {
  if(is.null(df) || nrow(df) == 0) return("—")
  paste0(round(df$mean[1], 5), " ± ", round(df$sd[1], 5), " (", 
         round(df$`0.025quant`[1], 5), ", ", round(df$`0.975quant`[1], 5), ")")
}


# -----------------------------


#' diagnostics
#' @noRd
.jmbm_diagnostics <- function(fit,
                             fam,
                             data_used,
                             n_glo = 0L,
                             priors = NULL,
                             pred_Sre = NULL,
                             pred_Sshared = NULL,
                             coupling.intercept,
                             scale_params_glo = NULL,
                             scale_params_reg = NULL) {

  has_Sre <- "Sre" %in% names(fit$summary.random)
  has_Sshared <- "Sshared" %in% names(fit$summary.random)

  # model fit metrics
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
  range_res_mean <- if(has_Sre) hyper["Range for Sre", "mean"] else "—"
  sigma_res_mean <- if(has_Sre) hyper["Stdev for Sre", "mean"] else "—"
  range_lat_mean <- if(has_Sshared) hyper["Range for Sshared", "mean"] else "—"
  sigma_lat_mean <- if(has_Sshared) hyper["Stdev for Sshared", "mean"] else "—"
  range_res_sd <- if(has_Sre) hyper["Range for Sre", "sd"] else "—"
  sigma_res_sd <- if(has_Sre) hyper["Stdev for Sre", "sd"] else "—"
  range_lat_sd <- if(has_Sshared) hyper["Range for Sshared", "sd"] else "—"
  sigma_lat_sd <- if(has_Sshared) hyper["Stdev for Sshared", "sd"] else "—"

  range_res <- if(has_Sre) paste0(round(range_res_mean, 2), " ± ", round(range_res_sd, 2)) else "—"
  sigma_res <- if(has_Sre) paste0(round(sigma_res_mean, 2), " ± ", round(sigma_res_sd, 2)) else "—"
  range_lat <- if(has_Sshared && is.finite(range_lat_mean)) paste0(round(range_lat_mean, 2), " ± ", round(range_lat_sd, 2)) else "—"
  sigma_lat <- if(has_Sshared && is.finite(sigma_lat_mean)) paste0(round(sigma_lat_mean, 2), " ± ", round(sigma_lat_sd, 2)) else "—"
  format_hyper <- function(param) {
    if(!param %in% rownames(hyper)) return("—")
      paste0(round(hyper[param, "mean"], 2), " ± ", round(hyper[param, "sd"], 2))
    }
  prec_IGlobal <- if("Precision for IGlobal" %in% rownames(hyper)) format_hyper("Precision for IGlobal") else "—"
  beta_IRegional <- if("Beta for IRegional" %in% rownames(hyper)) format_hyper("Beta for IRegional") else "—"
  beta_pred_oh <- if(!is.null(hyper)) {
    oh_rows <- rownames(hyper)[grepl("^Beta for .+RE_oh$", rownames(hyper))]
    if(length(oh_rows) > 0) {
      stats::setNames(lapply(oh_rows, format_hyper),oh_rows)
    } else NULL
  } else NULL

  # nested_shrinkage
  rs_delta_names <- names(fit$summary.random)[grepl("^.+_delta$", names(fit$summary.random))]
  beta_pred_rs <- if(length(rs_delta_names) > 0) {
    stats::setNames(lapply(rs_delta_names, function(nm) .mean_sd_str(fit$summary.random[[nm]])), rs_delta_names)
  } else NULL
  prec_pred_rs <- if(length(rs_delta_names) > 0 && !is.null(hyper)) {
    stats::setNames(lapply(rs_delta_names, function(nm) {
      hp_name <- paste0("Precision for ", nm)
      if(hp_name %in% rownames(hyper)) format_hyper(hp_name) else "—"
    }), rs_delta_names)
  } else NULL

  # Sre and Sshared overlap?
  if(is.finite(range_res_mean) && is.finite(range_lat_mean)) {
    ratio <- range_lat_mean / range_res_mean
  }

  # Intercepts
  iGlobal <- .mean_sd_str(fit$summary.random$IGlobal)
  iRegional <- .mean_sd_str(fit$summary.random$IRegional)

  # significant vars
  sig_vars <- .signif_vars(fit,  scale_params_glo = scale_params_glo,  scale_params_reg = scale_params_reg)

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

  # Brier, pred correlation, RMSE, Brier Skill Score
  brier <- mean((y_pred - y_obs)^2, na.rm = TRUE)
  pred_cor <- stats::cor(y_obs, y_pred, use = "complete.obs")
  rmse <- sqrt(mean((y_obs - y_pred)^2, na.rm = TRUE))
  prevalence <- mean(y_obs, na.rm = TRUE)
  brier_ref <- mean((prevalence - y_obs)^2, na.rm = TRUE)
  bss <- if(is.finite(brier_ref) && brier_ref > 0) 1 - (brier / brier_ref) else NA_real_

  # pit
  pit_vals <- fit$cpo$pit
  if(!is.null(pit_vals)) {
    pit_reg <- pit_vals[idx_reg]
    pit_reg <- pit_reg[is.finite(pit_reg)]
    ks_pit <- if(length(pit_reg) > 3) suppressWarnings(ks.test(pit_reg, "punif")$p.value) else "—"
  } else {
    ks_pit <- "—"
  }

  # Moran's I residual autocorrelation
  moran_I <- NA_real_
  rs <- y_obs - y_pred
  valid_idx <- is.finite(rs)
  rs <- rs[valid_idx]
  xy <- as.matrix(data_used[valid_idx, c("x", "y")])
  maxdist <- if(has_Sre) {
    range_res_mean
  } else if(has_Sshared) {
    range_lat_mean
  } else {
  # range unknown: default Moran's I threshold to 1/4 of bounding box diagonal
    bb <- apply(xy, 2, range, na.rm = TRUE)
    sqrt(sum((bb[2,] - bb[1,])^2)) / 4
  }
  if(is.finite(maxdist) && maxdist > 0) {
    nb <- spdep::dnearneigh(xy, 0, maxdist, longlat = FALSE)
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    mi <- spdep::moran(rs, lw, n = length(rs), S0 = spdep::Szero(lw))
    moran_I <- as.numeric(mi$I)
  }

  # calibration & coverage (disabled: only binomial supported)
  cal_slope <- "—"
  ks_pit <- "—"
  cov50 <- "—"
  cov95 <- "—"

  # CI/median ratios
  ci_ratios <- c(range_Sre = ci_ratio("Range for Sre"),
                 range_Sshared = ci_ratio("Range for Sshared"),
                 sigma_Sre = ci_ratio("Stdev for Sre"),
                 sigma_Sshared = ci_ratio("Stdev for Sshared"))

  # credibility interval ratios (CI/median)
  max_CIratio <- if(any(is.finite(ci_ratios))) max(ci_ratios, na.rm = TRUE) else NA_real_

  # posterior ≈ prior
  # prior sensitivity: contraction = 1 - (SD_post / SD_prior). 
  # Near 0 = prior dominance; near 1 = strong data update (Gelman et al. 2008).
  prior_sd_fixed <- 1 / sqrt(0.001)  # INLA default fixed effect prior precision
  prior_sd_linear <- 1 / sqrt(0.001) # 31.62: inlabru model='linear' default
  prior_close <- list(
    fixed_effects = if(!is.null(fit$summary.fixed) && nrow(fit$summary.fixed) > 0) {
    post_sds    <- fit$summary.fixed[["sd"]]
    contraction <- 1 - (post_sds / prior_sd_linear)
    uninformative <- rownames(fit$summary.fixed)[
      is.finite(contraction) & contraction < 0.1]
    length(uninformative) > 0
  } else FALSE,
  random_effects = FALSE
  )

  # variance and range ratios
  # sigma Sshared mean / sigma Sre mean > 1.5 --> Sshared field dominates (Blangiardo & Cameletti 2015)
  sigma_ratio <- if(has_Sre && has_Sshared && is.finite(sigma_lat_mean) && is.finite(sigma_res_mean) && sigma_res_mean > 0)
    sigma_lat_mean / sigma_res_mean else "—"
  # 0.5 < range_Sshared / range_Sre < 2 --> poor scale separation (Bakka et al. 2018)
  range_ratio <- if(has_Sre && has_Sshared && is.finite(range_lat_mean) && is.finite(range_res_mean) && range_res_mean > 0)
    range_lat_mean / range_res_mean else "—"

  var_explained_Sre <- "—"
  ssi <- NA_real_
  r2_fields <- "—"

  if(has_Sre && has_Sshared) {
    
    if(!is.null(marginals$`Stdev for Sre`) && !is.null(marginals$`Stdev for Sshared`)) {
      var_Sre_contrib <- INLA::inla.emarginal(function(x) x^2, marginals$`Stdev for Sre`)
      var_Sshared_contrib <- INLA::inla.emarginal(function(x) x^2, marginals$`Stdev for Sshared`)
      var_explained_Sre <- var_Sre_contrib / (var_Sre_contrib + var_Sshared_contrib)
    }

    # SSI (scale separation index)
    if(is.numeric(range_ratio) && is.numeric(sigma_ratio) && is.finite(range_ratio) && is.finite(sigma_ratio)) {
      range_penalty <- if(range_ratio > 3) 1 else range_ratio / 3
      sigma_penalty <- min(sigma_ratio, 1/sigma_ratio)
      ssi <- range_penalty * sigma_penalty
    }

    # Posterior field redundancy r^2(S_RE, S_shared). 
    # r^2 > 0.5 indicates both fields capture the same spatial pattern (redundancy)
    resid_Sre <- fit$summary.random$Sre$mean
    resid_Sshared <- fit$summary.random$Sshared$mean
    if(!is.null(resid_Sre) && !is.null(resid_Sshared)) {
      n_min <- min(length(resid_Sre), length(resid_Sshared))
      cor_fields <- stats::cor(resid_Sre[1:n_min], resid_Sshared[1:n_min], use = "complete.obs")
      r2_fields <- cor_fields^2
    }
  }

  # correlation Sre–Sshared fields
  # r > 0.7 00> high correlation
  field_correlation <- "—"
  if(has_Sre && has_Sshared) {
    f_sp <- fit$summary.random$Sre$mean
    f_lat <- fit$summary.random$Sshared$mean
    n <- min(length(f_sp), length(f_lat))
    field_correlation <- stats::cor(f_sp[seq_len(n)], f_lat[seq_len(n)], use = "pairwise.complete.obs")
  } 

  # warnings
  # overlap Sshared vs Sre
  if(has_Sre && has_Sshared && is.finite(range_ratio) && range_ratio > 0.5 && range_ratio < 3) {
    .warn(paste0(
      "Insufficient Sre scale separation: Sshared/Sre range ratio = ",round(range_ratio, 2), ").\n",
      "   Scales are overlapping (ideal: ratio > 3–5).\n",
      "   Recomended actions: (1) Tighten regional.pcprior.range[1] (e.g., 2 → 1.5),\n",
      "   (2) Increase shared.pcprior.range[1] (e.g., 50 → 100),\n",
      "   (3) Run prior sensitivity analysis,\n",
      "   (4) Evaluate if study domain truly supports two scales.\n"
    ))
  }  
  # weak identifiability (CI/median ratio)
  if(is.finite(max_CIratio) && max_CIratio > 25) {
    .warn(paste0(
      "Weak hyperparameter identifiability detected (CI/median > 25).\n",
      "   Recommended: use stronger PC priors."
    ))
  }
  # Sshared field dominating variance
  if(has_Sshared && has_Sre && is.finite(sigma_ratio) && sigma_ratio > 1.5) {
    .warn(paste0(
      "Variance imbalance between Sshared and Sre SPDE. Variance ratio (sigma Sshared / sigma Sre) = ", round(sigma_ratio, 2), ".\n",
      "   The Sshared field dominates the Sre variability, making the Sre SPDE redundant.\n",
      "   Recommended: reduce shared.pcprior.sigma or increase regional.pcprior.sigma."
    ))
  }
  # high correlation between Sshared and Sre
  if(has_Sre && has_Sshared && is.finite(field_correlation) && field_correlation > 0.7) {
    .warn(paste0(
      "High correlation between Sshared and Sre fields (r = ", round(field_correlation, 2), ").\n",
      "   Possible redundancy: both SPDE components capture the same spatial pattern.\n",
      "   Recommended: strengthen priors to separate scales, or remove one SPDE fields."
    ))
  }
  # residual autocorrelation
  if(is.finite(moran_I) && moran_I > 0.10) {
    .warn(paste0(
      "Residual spatial autocorrelation detected (Moran's I ≈ ", round(moran_I, 2), ").\n",
      "   Model missing local spatial structure.\n",
      "   Recommended: add a local SPDE component, refine mesh resolucion (smaller max.edges), or include missing covariates."
    ))
  }
  # posterior ≈ prior (weak data information)
  if(any(unlist(prior_close))) {
    .warn(paste0(
      "Posterior matches the uninformative prior: data provided little new information.\n",
      "   Posterior SD is close to prior SD (default N(0, sd=31.6) for linear effects).\n",
      "   This usually indicates near-constant covariates, very small sample size, or model misspecification.\n",
      "   Recommended: check covariate variability, increase sample size, or simplify the model."
    ))
  }    
  # cpo
  if(!is.null(fit$cpo$cpo)) {
    cpo_failures <- sum(fit$cpo$failure, na.rm = TRUE)
    total_obs <- length(fit$cpo$cpo)
    if(cpo_failures > 0) {
      perc_failures <- (cpo_failures / total_obs) * 100
      if(perc_failures > 1) {  # 1% threshold  (Held et al. 2010)
        .warn(paste0(
          "CPO failures detected (", round(perc_failures, 2), "% of observations).\n",
          "   Model severely struggles to predict these points (CPO ≈ 0).\n",
          "   Recommended: check for outliers or review model specification/priors."
        ))
      }
    }
  }

  # plots
  # Hyperparameters: posterior marginals + PC priors
  post_df <- data.frame()
  if(!is.null(fit$marginals.hyperpar)) {
    for(nm in names(fit$marginals.hyperpar)) {
      sm <- tryCatch(INLA::inla.smarginal(fit$marginals.hyperpar[[nm]]), error = function(e) NULL)
      if(!is.null(sm)) {
        post_df <- rbind(post_df, data.frame(x = sm$x, y = sm$y, par = nm))
      } else {
        .warn(paste0("Degenerate posterior marginal for '", nm, "'."))
      }
    }
  }
  if(nrow(post_df) == 0) post_df <- NULL
  prior_ticks <- do.call(rbind, Filter(Negate(is.null), list(
    if(!is.null(priors$regional.pcprior.range)) data.frame(x = priors$regional.pcprior.range[1], par = "Range for Sre"),
    if(!is.null(priors$regional.pcprior.sigma)) data.frame(x = priors$regional.pcprior.sigma[1], par = "Stdev for Sre"),
    if(!is.null(priors$shared.pcprior.range)) data.frame(x = priors$shared.pcprior.range[1], par = "Range for Sshared"),
    if(!is.null(priors$shared.pcprior.sigma)) data.frame(x = priors$shared.pcprior.sigma[1], par = "Stdev for Sshared")
  )))
  if(!is.null(prior_ticks) && nrow(prior_ticks) > 0 && !is.null(post_df)) {
    ymax_fac <- aggregate(y ~ par, post_df, function(z) max(z, na.rm = TRUE))
    prior_ticks$y <- ymax_fac$y[match(prior_ticks$par, ymax_fac$par)] * 0.95
  }

  ordered_hierarchical <- any(grepl("copy", rownames(fit$summary.hyperpar), ignore.case = TRUE)) ||
                   "IGlobal" %in% names(fit$summary.random)

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
  attr(cor_df, "range_eff") <- range_eff

  # Residual histogram + QQ plot
  res_mean <- mean(rs, na.rm = TRUE)
  df_r <- data.frame(resid = rs)
  center_label <- if(ordered_hierarchical) {
    paste0("Mean residual (≈ ", round(res_mean, 3), ")")
  } else {
    "Model mean = 0"
  }
  line_x <- if(ordered_hierarchical) res_mean else 0

  attr(df_r, "center_label") <- center_label
  attr(df_r, "line_x") <- line_x
  qq <- qqnorm(rs, plot.it = FALSE)
  df_qq <- data.frame(theoretical = qq$x, sample = qq$y)

  # semivariaogram
  coords <- as.matrix(data_used[, c("x", "y")])
  if(has_Sre) {
    sp_vals <- fit$summary.random$Sre$mean
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
  attr(sv_df, "range_res_mean") <- if(has_Sre) range_res_mean else NA_real_
  attr(sv_df, "range_lat_mean") <- if(has_Sshared) range_lat_mean else NA_real_

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
                       copy_beta_predictors = beta_pred_oh,
                       nested_shrinkage_delta = beta_pred_rs,
                       nested_shrinkage_delta_precision = prec_pred_rs),
    fixed_covariates = sig_vars,
    predictive = list(auc_full = auc_full,
                       tjur_r2 = tjur_r2,
                       brier = brier,
                       bss = bss,
                       rmse = rmse,
                       corr_obs_pred = pred_cor),
    calibration = list(slope = cal_slope,
                       ks_pit = ks_pit,
                       pit_values = if(exists("pit_reg")) pit_reg else NULL),
    coverage = list(cov50 = cov50,
                    cov95 = cov95),
    diagnostics = list(moran_I = moran_I,
                       max_CIratio = max_CIratio,
                       range_ratio = range_ratio,
                       sigma_ratio = sigma_ratio,
                       field_correlation = field_correlation,
                       var_explained_Sre = var_explained_Sre,
                       ssi = ssi,
                       r2_fields = r2_fields), 
    diagnostic_data = list(
      hyperparams = list(posteriors = post_df, prior_ticks = prior_ticks),
      correlogram = cor_df,
      hist = df_r,
      qq = df_qq,
      semivariogram = sv_df
    )
  )

  return(out)
}


# -----------------------------


#' prepare significant covariates
#' @noRd
.signif_vars <- function(fit, scale_params_glo = NULL, scale_params_reg = NULL) {
  sf <- fit$summary.fixed

  # recover GL-side and ordered_hierarchical regional-side iid covariates into Fixed effects
  gl_iid_names <- names(fit$summary.random)[grepl("GL$|RE_oh$", names(fit$summary.random))]
  if(length(gl_iid_names) > 0) {
    gl_iid_rows <- do.call(rbind, lapply(gl_iid_names, function(nm) {
      row <- fit$summary.random[[nm]]
      cols <- intersect(colnames(sf), colnames(row))
      row[1, cols, drop = FALSE]
    }))
    rownames(gl_iid_rows) <- gl_iid_names
    sf <- rbind(sf, gl_iid_rows)
  }

  has_params <- !is.null(scale_params_glo) || !is.null(scale_params_reg)

  # back-transform coef to original scale if vars were standardized
  if(has_params && nrow(sf) > 0) {
    for(i in seq_len(nrow(sf))) {
      coef_name <- rownames(sf)[i]
      # extract base variable name (remove GL, RE, RE_oh, GL_glo_res, RE_reg_anom suffixes)
      base_var <- gsub("GL$|RE$|RE_oh$|GL_glo_res$|RE_reg_anom$", "", coef_name)
               is_global <- grepl("GL$|GL_glo_res$", coef_name)
               params <- if(is_global) scale_params_glo else scale_params_reg

      if(!is.null(params) && base_var %in% names(params) && params[[base_var]]$sd > 0) {
        sd_x <- params[[base_var]]$sd
        sf[i, "mean"] <- sf[i, "mean"] / sd_x
        sf[i, "sd"] <- sf[i, "sd"] / sd_x
        sf[i, "0.025quant"] <- sf[i, "0.025quant"] / sd_x
        sf[i, "0.975quant"] <- sf[i, "0.975quant"] / sd_x
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

  # 95% CI no corss 0
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
.jmbm_generate_summary <- function(fit, species, fam, lnk, coupling.intercept, coupling.covariates, diag_block, cv_res=NULL, vg=NULL, vr=NULL, scale_params=NULL, scale_params_glo=NULL, scale_params_reg=NULL, has_spatial=FALSE) {
  
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

  has_Sre <- "Sre" %in% names(fit$summary.random)
  has_Sshared <- "Sshared" %in% names(fit$summary.random)
  has_covariates <- nrow(fit$summary.fixed) > 0

  # metadata
  species_name <- gsub("\\.", " ", species)
  base_name <- if (!has_spatial) "Non-spatial baseline" else "MBM"
  suf <- c()
  if(has_Sre) suf <- c(suf, "Sre")
  if(has_Sshared) suf <- c(suf, "Sshared")
  if(has_covariates) suf <- c(suf, "covariates")

  model_type <- if(length(suf) > 0) {
    paste0(base_name, ": ", paste(suf, collapse = " + "))
  } else {
    base_name
  }

  cp_int <- if(is.null(coupling.intercept)) "NULL" else coupling.intercept
  cp_pred <- if(is.list(coupling.covariates)) "custom list" else if(is.null(coupling.covariates)) "NULL" else coupling.covariates
  
  tbl_metadata <- data.frame(
    Field = c("Species name:", "Model type:", "Family | Link:", 
              "Coupling (Intercept):", "Coupling (Covariates):"),
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

  if(has_Sre) {
    params <- c(
      params,
      "Sre field: Range (posterior mean ± SD)",
      "Sre field: Sigma (posterior mean ± SD)"
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
      "Sshared field: Range (posterior mean ± SD)",
      "Sshared field: Sigma (posterior mean ± SD)"
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
  #nested_shrinkage predictors: delta_RE (regional deviation from the shared beta_GL) plus its estimated precision, one pair per variable.
  if(!is.null(int_block$nested_shrinkage_delta)) {
    for(nm in names(int_block$nested_shrinkage_delta)) {
      var_name <- sub("_delta$", "", nm)
      rows[[paste0("Random slope δ (", var_name, "RE deviation from ", var_name, "GL, mean ± SD, CI95%)")]] <- int_block$nested_shrinkage_delta[[nm]]
      prec_val <- int_block$nested_shrinkage_delta_precision[[nm]]
      if(!is.null(prec_val) && prec_val != "—") {
        rows[[paste0("Random slope δ precision (", var_name, ", mean ± SD)")]] <- prec_val
      }
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
    base_vars <- gsub("GL$|RE$", "", tbl_fixed$coef)
    
    # which coupling?
    tbl_fixed$Coupling <- vapply(base_vars, function(v) {
      if(v %in% c(vg, vr)) {
        .resolve_coupling_predictor(v, coupling.covariates, vg, vr)
      } else {
        "—"
      }
    }, character(1))
    
    cols_order <- c("coef", "coupling", "estimate", "sd", "2.5%", "97.5%", "signif", "tail_p")
    cols_order <- intersect(cols_order, names(tbl_fixed))
    tbl_fixed <- tbl_fixed[, cols_order]

    is_gl <- grepl("GL$|GL_glo_res$", tbl_fixed$coef)
    base_n <- gsub("GL$|GL_glo_res$|RE$|RE_oh$|RE_reg_anom$", "", tbl_fixed$coef)
    idx_gl <- which( is_gl)[order(base_n[ is_gl])]
    idx_re <- which(!is_gl)[order(base_n[!is_gl])]
      tbl_fixed <- tbl_fixed[c(idx_gl, idx_re), ]
      rownames(tbl_fixed) <- NULL
  }

  # predictive performance
  tbl_pred <- data.frame(
    Metric = c("AUC (full model)",
               "Tjur R\u00b2 (discrimination coefficient)",
               "Brier score",
               "Brier Skill Score",
               "RMSE",
               "Observed-predicted correlation (r)"),
    Value  = c(diag_block$predictive$auc_full,
               fmt_val(diag_block$predictive$tjur_r2),
               fmt_val(diag_block$predictive$brier),
               fmt_val(diag_block$predictive$bss),
               fmt_val(diag_block$predictive$rmse),
               fmt_val(diag_block$predictive$corr_obs_pred)),
    stringsAsFactors = FALSE
  )

  # calibration & coverage (disabled: only binomial supported)
  tbl_cal <- NULL

  # diagnostics
  tbl_diag <- data.frame(
    Metric = c(
      "Residual Moran's I",
      "Max CI/median ratio (hyperparameters)",
      "Scale-separation ratio (range_Sshared / range_Sre)",
      "Variance ratio (sigma_Sshared / sigma_Sre)",
      "Sshared–Sre field correlation (r)",
      "Variance explained by Sre (%)",
      "Scale separation Index (SSI) [0–1]",
      "Posterior field redundancy r\u00b2(S_RE, S_shared) [0–1]"
    ),
    Value = c(
      fmt_val(diag_block$diagnostics$moran_I),
      fmt_val(diag_block$diagnostics$max_CIratio),
      fmt_val(diag_block$diagnostics$range_ratio),
      fmt_val(diag_block$diagnostics$sigma_ratio),
      fmt_val(diag_block$diagnostics$field_correlation),
      if(is.numeric(diag_block$diagnostics$var_explained_Sre)) 
        fmt_val(diag_block$diagnostics$var_explained_Sre * 100) else "—",
      fmt_val(diag_block$diagnostics$ssi),
      fmt_val(diag_block$diagnostics$r2_fields)
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
    `Random effects` = tbl_intercepts,
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
.jmbm_vars_importance <- function(fit) {

  df <- fit$summary.fixed
  df$var <- rownames(df)

  # exclude intercepts and Sshared
  drop <- c("IGlobal", "IRegional", "beta_GL")
  df <- df[!(df$var %in% drop), , drop = FALSE]

  df$sign <- ifelse(df$mean > 0, "Positive", "Negative")

  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = reorder(var, abs(mean)),
                 y = abs(mean),
                 fill = sign)) + 
    ggplot2::geom_col(alpha = 0.9) +
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


