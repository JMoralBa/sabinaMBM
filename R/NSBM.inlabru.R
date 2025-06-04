#' @name NSBM.inlabru
#'
#' @title Nested species distribution modeling ...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru.
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param mesh An INLA mesh object created externally with `create_mesh()`. Required only if \code{spatial = TRUE}. Ignored if \code{spatial = FALSE}.
#' @param output A character `"intensity"` or `"probability"` (default).
#' @param prior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param prior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param spatial Logical. Include spatial latent field (SPDE) in the model (default: TRUE). 
#' @param proj.new.env Logical. Whether to compute predictions under new scenarios (default: TRUE).
#' @param seed Optional integer. If provided, sets a random seed for reproducibility.
#' @param save.output Logical. If TRUE, saves key model outputs (predictions, evaluation, summary ...).
#'
#' @return A named list of class `nsbm.inlabru` with the following elements:
#' \item{Species.Name}{Species name}
#' \item{args}{List of arguments used in the model fitting.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{current.projections}{List with: fitted model (`fit`), prediction (`pred`), spatial field (`pred_sp`).}
#' \item{new.projections}{List of projections to new.env (if `proj.new.env = TRUE`).}
#' \item{Summary}{\code{data.frame}  with key evaluation metrics and significant variables.} #@@@JMB revisar y refinar
#'
#' @export
NSBM.inlabru <- function(nsbm_obj, 
                         mesh = NULL, 
                         output = "probability",
                         spatial = TRUE,
                         prior.range = c(5, 0.01),
                         prior.sigma = c(1, 0.01),
                         proj.new.env = TRUE,
                         seed = NULL,
                         save.output = FALSE) {

  if(!inherits(nsbm_obj, "nsdm.vinput")) {
    stop("The 'nsbm_obj' must be of class 'nsdm.vinput', Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  if(!(output %in% c("probability", "intensity"))) {
    stop("Invalid 'output'. Please, use 'probability' or 'intensity'.")
  }
  if(spatial && is.null(mesh)) {
    stop("If spatial = TRUE, you must provide a mesh object using create_mesh().")
  }
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a numeric value.")
    }
    set.seed(seed)
  }

  sabina <- list()

  # Data preparation
  pp_regional <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Regional, coords = c("x", "y"))
  pp_global   <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  sp_covglo   <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg   <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)

  crs <- sf::st_crs(sp_covglo)
  sf::st_crs(pp_global) <- crs
  sf::st_crs(pp_regional) <- crs
  
  pp_global <- sf::st_transform(pp_global, crs)
  pp_regional <- sf::st_transform(pp_regional, crs)

  pts_reg <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs) 

  geom_comb <- c(sf::st_geometry(pp_global), sf::st_geometry(pts_reg))
  aux <- sf::st_sf(geometry = geom_comb)
  sf::st_crs(aux) <- crs
  
  # Spatial domain definition
  bdy_global <- sf::st_convex_hull(sf::st_union(aux))
  bdy_regional <- sf::st_union(sf::st_make_valid(sf::st_as_sf(raster::rasterToPolygons(raster::raster(sp_covreg)))))
  sf::st_crs(bdy_global) <- crs
  sf::st_crs(bdy_regional) <- crs

  # Description SPDE
  if(spatial) {
    matern <- INLA::inla.spde2.pcmatern(
      mesh,
      prior.range = prior.range,
      prior.sigma = prior.sigma
    )
  }

  # Model components
  cmp_cov <- fcov(nsbm_obj, "sp_covglo", "sp_covreg")
  cmp_formula <- if(spatial) {
    paste0("~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ IGlobal(1) + IRegional(1) + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  if(output == "probability") {
    pp_regional$presence <- 1
    pseudo_regional <- sf::st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
    pseudo_regional$presence <- 0
    sf::st_crs(pseudo_regional) <- crs
    pseudo_regional <- sf::st_transform(pseudo_regional, crs)
    pp_regional <- rbind(pp_regional, pseudo_regional)

    pp_global$presence <- 1
    pseudo_global <- sf::st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
    sf::st_crs(pseudo_global) <- crs
    pseudo_global$presence <- 0
    pseudo_global <- sf::st_transform(pseudo_global, crs) 
    pp_global <- rbind(pp_global, pseudo_global)

    # Likelihoods
    f_spatial <- if(spatial) " + spatial" else ""

    lik_global <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- if(spatial) {
      as.formula(paste0("~ 1 / (1 + exp(-(IRegional + spatial + ", cmp_cov$like$fregional, ")))"))
      } else {
        as.formula(paste0("~ 1 / (1 + exp(-(IRegional + ", cmp_cov$like$fregional, ")))"))
      }

  } else {
    lik_global <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- if(spatial) {
      as.formula(paste0("~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")"))
    } else {
      as.formula(paste0("~ exp(IRegional + ", cmp_cov$like$fregional, ")"))
    }
  }

  # Model fitting
  fit <- inlabru::bru(
    components = cmp,
    lik_global,
    lik_regional,
    options = list(
      control.compute = list(
        cpo  = TRUE,
        waic = TRUE,
        dic  = TRUE
      )
    )
  )

  # log sum of conditional predictive ordinates
  lcpo_val <- if(!is.null(fit$cpo$cpo)) round(sum(log(fit$cpo$cpo)), 2) else "Not computed"

  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)

  # Predictions
  pred <- predict(fit, pred.df, pred_formula)  
  pred <- pred_as_tif(pred, sp_covreg) # from sf to tif

  if(spatial) {
    pred_sp <- predict(fit, pred.df, ~ spatial)
    pred_sp <- pred_as_tif(pred_sp, sp_covreg) # from sf to tif
  } else {
    pred_sp <- NULL
  }

  # New scenarios
  proj_list <- list()
  if (proj.new.env && !is.null(nsbm_obj$Scenarios)) {
    for (sc in names(nsbm_obj$Scenarios)) {
      #sc <- 1
      scen_rast <- terra::unwrap(nsbm_obj$Scenarios[[sc]])
      scen_df <- sf::st_as_sf(as.points(scen_rast))
      sf::st_crs(scen_df) <- crs
      scen_df <- sf::st_transform(scen_df, crs)
      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[paste0("proj_", sc)]] <- pred_as_tif(proj_pred, template = scen_rast)  # from sf to tif
    }
  }

  species <- nsbm_obj$Species.Name

  # save outputs
  if(save.output) {
    base_results <- "Results"
    nsbm_path <- file.path(base_results, "NSBM")
    values_path <- file.path(nsbm_path, "Values")
    projections_path <- file.path(nsbm_path, "Projections")

    # Create directories
    fs::dir_create(values_path, recurse = TRUE)
    fs::dir_create(projections_path, recurse = TRUE)

    # save fixed effects
    if(!is.null(fit$summary.fixed)) {
      write.csv(fit$summary.fixed, file = file.path(values_path, paste0(species, "_fixed_effects.csv")), row.names = TRUE)
    }

    # save spatial random effects
    if(!is.null(fit$summary.random)) {
       for(ran in names(fit$summary.random)) {
        write.csv(fit$summary.random[[ran]], file = file.path(values_path, paste0(species, "_random_", ran, ".csv")), row.names = TRUE)
      }
    }

    # save hypermarams
    if(!is.null(fit$summary.hyperpar)) {
      write.csv(fit$summary.hyperpar, file = file.path(values_path, paste0(species, "_hyperparameters.csv")), row.names = TRUE)
    }

    # save evaluation metrics
    eval_metrics <- data.frame(
      Metric = c("WAIC", "DIC", "MLik", "LCPO (sum log-CPO)"),
      Value = c(fit$waic$waic, fit$dic$dic, fit$mlik[1], lcpo_val)
    )
    write.csv(eval_metrics, file = file.path(values_path, paste0(species, "_evaluation.csv")), row.names = FALSE)

    # Save CPO values (one per observation)   #@@@JMB useful for leave-one-out diagnostics or model comparison??
    write.csv(data.frame(CPO = fit$cpo$cpo), file = file.path(values_path, paste0(species, "_pointwise_CPO.csv")), row.names = FALSE)

    # save full model object (fit)
    saveRDS(fit, file = file.path(values_path, paste0(species, "_model_fit.rds")))

    # save pred current
    if(!is.null(pred)) {
      file_path <- file.path(projections_path, paste0(species, "_Current.tif"))
      terra::writeRaster(terra::unwrap(pred), file_path, overwrite = TRUE)   
    }

    # save pred_sf
    if(!is.null(pred_sp)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_Spatial.tif"))
      terra::writeRaster(terra::unwrap(pred_sp), file_path, overwrite = TRUE)
    }

    # save new scenarios
    if(length(proj_list) > 0 && !is.null(nsbm_obj$Scenarios)) {
      for(i in seq_along(proj_list)) {
        sc_name <- names(proj_list)[i]
        file_path <- file.path(projections_path, paste0(species, "_", sc_name, ".tif"))
        terra::writeRaster(terra::unwrap(proj_list[[i]]), file_path, overwrite = TRUE)
      }
    }

    message("Results saved in the following local folder(s):")
    message(paste(
    "  - Current projection (pred), spatial field (pred_sp), and new scenarios: ", projections_path, "\n",
    " - Fixed/random spatial effects, hyperparameters, evaluation: ", values_path, "\n",
    " - Full model object: ", file.path(values_path, paste0(species, "_model_fit.rds")), "\n"
    ))
  }

  #
  sabina <- list(
      Species.Name = nsbm_obj$SpeciesName,
      args = list(
        spatial = spatial,
        prior.range = prior.range,
        prior.sigma = prior.sigma,
        proj.new.env = proj.new.env
      ),
      Selected.Variables.Global = nsbm_obj$Selected.Variables.Global,
      Selected.Variables.Regional = nsbm_obj$Selected.Variables.Regional,
      current.projections = list(
        fit = fit,
        pred = terra::wrap(pred),
        pred_sp = if (!is.null(pred_sp)) terra::wrap(pred_sp) else NULL
      ),
      new.projections = rapply(proj_list, terra::wrap, how = "list"),
      Summary = generate_summary_nsbm(fit, species, spatial, lcpo_val) #@@@JMB pendiente revisar/completar...
  )
 
  attr(sabina, "class") <- "nsbm.inlabru"
  return(sabina)

}


# formulas
fcov <- function(obj, spobjglo, spobjreg) {
  vars <- unique(c(obj$Selected.Variables.Global, obj$Selected.Variables.Regional))
  cmp1 <- paste(paste(vars, "(1)", sep = ""), collapse = " + ")

  cmpglobal <- paste(sapply(obj$Selected.Variables.Global, function(X) {
    paste0(X, "GL(main = ", spobjglo, ", main_layer = \"", X, "\", model = \"const\")")
  }), collapse = " + ")

  cmpregional <- paste(sapply(obj$Selected.Variables.Regional, function(X) {
    paste0(X, "RE(main = ", spobjreg, ", main_layer = \"", X, "\", model = \"const\")")
  }), collapse = " + ")

  fglobal <- paste(sapply(obj$Selected.Variables.Global, function(X) {
    paste0(X, " * ", X, "GL")
  }), collapse = " + ")

  fregional <- paste(sapply(obj$Selected.Variables.Regional, function(X) {
    paste0(X, " * ", X, "RE")
  }), collapse = " + ")

  list(
    cmp = paste(c(cmp1, cmpglobal, cmpregional), collapse = " + "),
    like = list(fglobal = fglobal, fregional = fregional)
  )
}



# prepare summary
generate_summary_nsbm <- function(fit, species=species, spatial, lcpo_val) {

  # Fixed effects and hyperpar
  summary_fixed <- fit$summary.fixed
  hyper         <- fit$summary.hyperpar

  # Filter valid vars
  valid_vars <- summary_fixed[!is.na(summary_fixed$mean), ]
  valid_vars$type <- ifelse(
    rownames(valid_vars) %in% c("IGlobal", "IRegional"), "Intercept",
    ifelse(grepl("GL$", rownames(valid_vars)), "Global",
           ifelse(grepl("RE$", rownames(valid_vars)), "Regional", "Unclassified"))
  )

  # Filter significant vars (CI does not cross 0 and mean is relevant)
  signif_vars <- valid_vars[
    valid_vars[,"0.025quant"] * valid_vars[,"0.975quant"] > 0 &
    abs(valid_vars[,"mean"]) > 0.05,
  ]

  # Order by abs mean (importance?)
  signif_vars <- signif_vars[order(-abs(signif_vars[,"mean"])), ]

  # Label
  var_labels <- paste0(
    rownames(signif_vars), " ",
    ifelse(signif_vars$mean > 0, "(+)", "(–)")
  )

  # Evaluation metrics
  dic_val  <- if (!is.null(fit$dic$dic) && !is.na(fit$dic$dic)) round(fit$dic$dic, 2) else "Not computed"
  waic_val <- if (!is.null(fit$waic$waic) && !is.na(fit$waic$waic)) round(fit$waic$waic, 2) else "Not computed"
  mlik_val <- if (!is.null(fit$mlik) && !is.na(fit$mlik[1,1])) round(fit$mlik[1, 1], 2) else "Not computed"
  
  # 
  spatial_range <- if (!is.null(hyper) && "Range for spatial" %in% rownames(hyper)) {
    paste0(round(hyper["Range for spatial", "mean"], 2), " ± ",
           round(hyper["Range for spatial", "sd"], 2))
  } else {
    "Not computed"
  }

  #
  summary_df <- data.frame(
    Field = c(
      "Species name:",
      "Model type:",
      "DIC;",
      "WAIC:",
      "Marginal log-likelihood:",
      # log sum of conditional predictive ordinates: a bayesian metric for model validation (lower values indicate better fit)
      "LCPO (sum log-CPO):",   
      "SPDE spatial range (mean ± sd):",
      "Significant variables (ordered):"
    ),
    Value = c(
      gsub("\\.", " ", species),
      if(spatial) {
        "NSBM.inlabru (with SPDE)"
      } else {
        "NSBM.inlabru (no SPDE)"
      },
      dic_val,
      waic_val,
      mlik_val,
      lcpo_val,
      spatial_range,
      if (length(var_labels) > 0) paste(var_labels, collapse = ", ") else "None"
    ),
    stringsAsFactors = FALSE
  )

  return(summary_df)
}


# from sf to tif
pred_as_tif <- function(pred, template, vars_to_export = c("mean", 
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
    if (!is.null(pred[[var]])) {
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
  if (length(r_stack) == 0) return(NULL)

  r_pred <- do.call(c, r_stack)
  return(r_pred)
}

