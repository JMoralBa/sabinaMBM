#' @name NSBM.covariate
#'
#' @title Covariate-based species distribution model....
#'
#' @description bla bla...
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param output Character; `"intensity"` or `"probability"` (default).
#' @param family  Character; `"binomial"` (default), supported `"poisson"`, `"nbinomial"`, or `"cp"`.
#' @param link  Character; link function `"logit"` (default for `family = "binomial"`), otherwise `"log"`.
#' @param prior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param prior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param spatial Logical. Include spatial latent field (SPDE) in the model (default: TRUE). 
#' @param mesh An INLA mesh object created externally with `create_mesh()`. Required if `spatial = TRUE`.
#' @param rm.corr (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} controlling whether environmental covariates correlated with the global model should be removed. The threshold value used for identifying collinearity is the same used with \code{\link{NSDM.SelectCovariates}} function.
#' #@param corcut (\emph{optional, default} \code{0.7}) \cr
#' #A \code{numeric} value for the correlation coefficient threshold used for identifying collinearity.   #@@@JMB ponerlo como argumento en lugar de traerlo de NSDM.SelectCovariates()???
#' @param proj.new.env Logical; whether to compute predictions under new scenarios (default: TRUE).
#' @param cv.folds Number of k-folds for cross-validation (default: 1 = no CV). If >1, returns mean ± sd AUC.
#' @param n.threads Number of threads for analysis
#' @param seed Optional integer to set the random seed for reproducibility.
#' @param save.output Logical; if TRUE, saves key model outputs (predictions, evaluation, summary ...).
#' @param save.independent Logical; if `TRUE`, also saves the individual global and regional model vlaues and predictions.
#'
#'
#' @return A named list of class `nsbm.inlabru, including:
#' \item{Species.Name}{Species name.}
#' \item{args}{List of arguments used in the function.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{Selected.Variables.Covariate}{Names of regional covariates after dropping high correlation.}
#' \item{current.projections}{List with prediction: \code{pred.covariate} and, if requested, separate \code{pred.global}.}
#' \item{new.projections}{List of combined projections under new environmental scenarios (if `proj.new.env = TRUE`).}
#' \item{Summary}{A `data.frame` with key evaluation metrics and significant variables.}
#'
#' @export
NSBM.covariate <- function(nsbm_obj,
                           output = "probability",
                           family = "binomial",
                           link = "logit",
                           prior.range = c(5, 0.01),
                           prior.sigma = c(1, 0.01),
                           spatial = TRUE,
                           mesh = NULL,
                           rm.corr = TRUE,
                           corcut = 0.7,     #@@@JMB podemos obligar a tomar el de nsbm_obj$args$corcut
                           proj.new.env = TRUE,
                           n.threads = 1,
                           cv.folds = 1,
                           seed = NULL,
                           save.output = FALSE,
                          save.independent = FALSE) {

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
  available_cores <- parallel::detectCores(logical = TRUE)
  if(!is.null(n.threads) && n.threads > available_cores) {
    stop(paste0("Requested n.threads = ", n.threads, " exceeds available hardware threads (", available_cores,")."))
  }
  INLA::inla.setOption(num.threads = n.threads)
  if(!is.null(seed)) set.seed(seed)

  sabina <- list()
  current.projections <- list(
    pred.covariate = NULL,
    pred_sp.covariate = NULL,
    pred.global = NULL,
    pred_sp.global = NULL
  )
  if (proj.new.env) {
    scens <- names(nsbm_obj$Scenarios)
    if (is.null(scens) || length(scens) == 0) {
      warning("proj.new.env = TRUE but no 'Scenarios' found; skipping new projections.")
      new.projections <- NULL
    } else {
      new.projections <- setNames(vector("list", length(scens)), paste0("pred.", scens))
    }
  } else {
    new.projections <- NULL
  }

  # Prepare data
  # GLOBAL
  pres_glo <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  abs_glo  <- sf::st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
  sp_covglo <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  crs <- sf::st_crs(sp_covglo)
  sf::st_crs(pres_glo) <- crs
  sf::st_crs(abs_glo) <- crs

  pres_glo <- sf::st_transform(pres_glo, crs)
  pres_glo$presence <- 1L
  abs_glo$presence <- 0L
  pp_glo <- rbind(pres_glo, abs_glo)

  # REGIONAL
  pres_reg <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Regional, coords = c("x", "y"))
  abs_reg <- sf::st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
  sp_covreg <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)
  sf::st_crs(pres_reg) <- crs
  sf::st_crs(abs_reg) <- crs

  pres_reg <- sf::st_transform(pres_reg, crs)
  pres_reg$presence <- 1L
  abs_reg$presence  <- 0L
  pp_reg <- rbind(pres_reg, abs_reg)

  pts_reg <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs) 

  geom_comb <- c(sf::st_geometry(pres_glo), sf::st_geometry(pts_reg))
  aux <- sf::st_sf(geometry = geom_comb)
  sf::st_crs(aux) <- crs

  # Spatial domain definition
  bdy_glo <- sf::st_convex_hull(sf::st_union(aux))
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(raster::rasterToPolygons(raster::raster(sp_covreg)))))
  sf::st_crs(bdy_glo) <- crs
  sf::st_crs(bdy_reg) <- crs

  # SPDE
  if(spatial) {
    matern <- INLA::inla.spde2.pcmatern(mesh, prior.range = prior.range, prior.sigma = prior.sigma)
  }

  # prefilter select terms
  #.....

  # cmp global
  cmp_glo <- indiv_fcov(vars = nsbm_obj$Selected.Variables.Global,
                        spobj = "sp_covglo",
                        tag = "GL")

  # Formula global
  cmp_formula_glo <- if(spatial) {
    paste0("~ IGlobal(1) + spatial(geometry, model = matern) + ", cmp_glo$cmp)
  } else {
    paste0("~ IGlobal(1) + ", cmp_glo$cmp)
  }
  cmp_formula_glo <- as.formula(cmp_formula_glo)

  if(output == "probability") {
    f_spatial <- if(spatial) " + spatial" else ""
    # Likelihoods
    lik_glo <- inlabru::like(
      family = family,
      formula = as.formula(paste0("presence ~ IGlobal ", f_spatial, " + ", cmp_glo$like)),
      data = pp_glo,
      samplers = bdy_glo,
      domain = list(geometry = mesh),
      control.family = list(link = link)
    )
    pred_formula_glo <- as.formula(paste0("~ 1 / (1 + exp(-(IGlobal", f_spatial, " + ", cmp_glo$like,")))"))
  } else {  # output = "intensity"
    lik_glo <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, " + ", cmp_glo$like)),
      data = pres_glo,
      samplers = bdy_glo,
      domain = list(geometry = mesh)
    )
    pred_formula_glo <- as.formula(paste0("~ exp(IGlobal", f_spatial, " + ", cmp_glo$like,")"))
  }

  # Fit global
  fit_glo <- inlabru::bru(
    components = cmp_formula_glo,
    lik_glo,
    options = list(control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE))
  )

  # log sum of conditional predictive ordinates
  lcpo_val_glo <- if(!is.null(fit_glo$cpo$cpo)) round(sum(log(fit_glo$cpo$cpo)), 2) else "Not computed"

  # Prediction global on regional scale layers
  sp_covglo_reg <- terra::unwrap(nsbm_obj$IndVar.Global.Selected.reg)
  pred.df <- sf::st_as_sf(as.points(sp_covglo_reg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)

  old_sp_covglo <- sp_covglo # Backup original global raster
  sp_covglo <- sp_covglo_reg # overwrite for predict() using the fine-scale raster

  pred_glo <- predict(fit_glo, pred.df, pred_formula_glo)

  sp_covglo <- old_sp_covglo # restore

  pred_glo <- pred_as_tif(pred_glo, sp_covglo_reg)       #@@@JMB save uncertainty?????

  if(spatial) {
    pred_sp_glo <- predict(fit_glo, pred.df, ~ spatial)
    pred_sp_glo <- pred_as_tif(pred_sp_glo, sp_covglo_reg)
  } else {
    pred_sp_glo <- NULL
  }

  # Add global model as new covariate of regional model.
  SDM.global <- pred_glo[["mean"]]
  names(SDM.global) <- "SDM.global"
  sp_covreg <- c(sp_covreg, SDM.global)

  # remove correlated covariates
  if(rm.corr) {
    coords <- sf::st_coordinates(pp_reg)
    myExpl.covsel <- terra::extract(sp_covreg, coords)
    # covsel.filteralgo
    sel_df <- covsel::covsel.filteralgo(
      covdata = myExpl.covsel,
      pa = as.vector(pp_reg$presence),
      force = "SDM.global",
      corcut = corcut
    )
    IndVar.Regional.Covariate <- intersect(names(sel_df), names(sp_covreg))
    dropped <- setdiff(names(sp_covreg), names(sel_df))
    sp_covreg <- sp_covreg[[IndVar.Regional.Covariate]]
    dropped_msg <- if(length(dropped) == 0) "None" else paste(dropped, collapse = ", ")
    message("\nRemoved covariates: ", dropped_msg, "\n")
  }

  # cmp regional/covariate
  cmp_cov <- indiv_fcov(vars = names(sp_covreg),
                        spobj= "sp_covreg",
                        tag = "RE")

  # Formula
  cmp_formula_cov <- if(spatial) {
    paste0("~ IRegional(1) + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ IRegional(1) + ", cmp_cov$cmp)
  }
  cmp_formula_cov <- as.formula(cmp_formula_cov)

  if(output == "probability") {
    f_spatial <- if(spatial) " + spatial" else ""
    # Likelihoods
    lik_cov <- inlabru::like(
      family = family,
      formula = as.formula(paste0("presence ~ IRegional ", f_spatial, " + ", cmp_cov$like)),
      data = pp_reg,
      samplers = bdy_reg,
      domain = list(geometry = mesh),
      control.family = list(link = link)
    )
    pred_formula_cov <- as.formula(paste0("~ 1 / (1 + exp(-(IRegional", f_spatial, " + ", cmp_cov$like,")))"))
  } else {  # output = "intensity"
    lik_cov <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", f_spatial, " + ", cmp_cov$like)),
      data = pres_reg,
      samplers = bdy_reg,
      domain = list(geometry = mesh)
    )
    pred_formula_cov <- as.formula(paste0("~ exp(IRegional", f_spatial, " + ", cmp_cov$like,")"))
  }

  # fit model covariate
  fit_cov <- inlabru::bru(
    components = cmp_formula_cov,
    lik_cov,
    options = list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
    )

  # log sum of conditional predictive ordinates
  lcpo_val_cov <- if(!is.null(fit_cov$cpo$cpo)) round(sum(log(fit_cov$cpo$cpo)), 2) else "Not computed"

  # Prediction covariate
  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)

  pred_cov <- predict(fit_cov, pred.df, pred_formula_cov)  
  pred_cov <- pred_as_tif(pred_cov, sp_covreg)

  if(spatial) {
    pred_sp_cov <- predict(fit_cov, pred.df, ~ spatial)
    pred_sp_cov <- pred_as_tif(pred_sp_cov, sp_covreg)
  } else {
    pred_sp_cov <- NULL
  }

  current.projections <- list(
    pred.covariate = terra::wrap(pred_cov),
    pred_sp.covariate = if(!is.null(pred_sp_cov)) terra::wrap(pred_sp_cov) else NULL
  )

  if(save.independent) {
    current.projections$pred.global <- terra::wrap(pred_glo)
    current.projections$pred_sp.global <- if(!is.null(pred_sp_glo)) terra::wrap(pred_sp_glo)
  }

  # New scenarios
  proj_list <- list()
  if(proj.new.env && !is.null(nsbm_obj$Scenarios)) {
    for (sc in names(nsbm_obj$Scenarios)) {
      scen_rast <- terra::unwrap(nsbm_obj$Scenarios[[sc]])
      scen_df <- sf::st_as_sf(as.points(scen_rast))
      sf::st_crs(scen_df) <- crs
      scen_df <- sf::st_transform(scen_df, crs)
      proj_pred <- predict(fit_cov, scen_df, pred_formula_cov)
      proj_list[[paste0("proj_", sc)]] <- pred_as_tif(proj_pred, template = scen_rast)  # from sf to tif
    }						#@@@JMB guardar uncertainty new env?
  }

  # k-fold CV for global and regional
  if (cv.folds > 1) {
    cv_res_glo <- cv_individual_inlabru(
      lik_obj = lik_glo,
      pp_data = pp_glo,
      cmp_formula = cmp_formula_glo,
      pred_formula = pred_formula_glo,
      cv.folds = cv.folds
    )
    cv_res_cov <- cv_individual_inlabru(
      lik_obj = lik_cov,
      pp_data = pp_reg,
      cmp_formula = cmp_formula_cov,
      pred_formula = pred_formula_cov,
      cv.folds = cv.folds
    )
  } else {
    cv_res_glo <- cv_res_cov <- NULL
  }

  # save covariate
  species <- nsbm_obj$Species.Name
  if(save.output) {
    # Create directories
    values_path <- file.path("Results", "NSBM_covariate", "Values")
    projections_path <- file.path("Results", "NSBM_covariate", "Projections")
    fs::dir_create(values_path, recurse = TRUE)
    fs::dir_create(projections_path, recurse = TRUE)

    # save fixed effects
    if(!is.null(fit_cov$summary.fixed)) {
      write.csv(fit_cov$summary.fixed, file = file.path(values_path, paste0(species, "_fixed_effects.csv")), row.names = TRUE)
    }

    # save spatial random effects
    if(!is.null(fit_cov$summary.random)) {
       for(ran in names(fit_cov$summary.random)) {
        write.csv(fit_cov$summary.random[[ran]], file = file.path(values_path, paste0(species, "_random_", ran, ".csv")), row.names = TRUE)
      }
    }

    # save hypermarams
    if(!is.null(fit_cov$summary.hyperpar)) {
      write.csv(fit_cov$summary.hyperpar, file = file.path(values_path, paste0(species, "_hyperparameters.csv")), row.names = TRUE)
    }

    # save evaluation metrics
    eval_metrics <- data.frame(
      Metric = c("WAIC", "DIC", "MLik", "LCPO (sum log-CPO)"),
      Value = c(fit_cov$waic$waic, fit_cov$dic$dic, fit_cov$mlik[1], lcpo_val_cov)
    )
    write.csv(eval_metrics, file = file.path(values_path, paste0(species, "_evaluation.csv")), row.names = FALSE)

    # Save CPO values (one per observation)   #@@@JMB useful for leave-one-out diagnostics or model comparison??
    write.csv(data.frame(CPO = fit_cov$cpo$cpo), file = file.path(values_path, paste0(species, "_pointwise_CPO.csv")), row.names = FALSE)

    # save full model object (fit)
    saveRDS(fit_cov, file = file.path(values_path, paste0(species, "_model_fit.rds")))

    # save pred current
    if(!is.null(pred)) {
      file_path <- file.path(projections_path, paste0(species, "_Current.tif"))
      terra::writeRaster(terra::unwrap(pred_cov), file_path, overwrite = TRUE)   
    }

    # save pred_sf
    if(!is.null(pred_sp_cov)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_Spatial.tif"))
      terra::writeRaster(terra::unwrap(pred_sp_cov), file_path, overwrite = TRUE)
    }

    # save new scenarios
    if(length(proj_list) > 0 && !is.null(nsbm_obj$Scenarios)) {
      for(i in seq_along(proj_list)) {
        sc_name <- names(proj_list)[i]
        file_path <- file.path(projections_path, paste0(species, "_", sc_name, ".tif"))
        terra::writeRaster(terra::unwrap(proj_list[[i]]), file_path, overwrite = TRUE)
      }
    }

    # save global
    if(save.independent) {
      # save fixed effects
      if(!is.null(fit_glo$summary.fixed)) {
        write.csv(fit_glo$summary.fixed, file = file.path(values_path, paste0(species, "_fixed_effects_glo.csv")), row.names = TRUE)
      }

      # save spatial random effects
      if(!is.null(fit_glo$summary.random)) {
        for(ran in names(fit_glo$summary.random)) {
          write.csv(fit_glo$summary.random[[ran]], file = file.path(values_path, paste0(species, "_random_", ran, "_global.csv")), row.names = TRUE)
        }
      }

      # save hypermarams
      if(!is.null(fit_glo$summary.hyperpar)) {
        write.csv(fit_glo$summary.hyperpar, file = file.path(values_path, paste0(species, "_hyperparameters_global.csv")), row.names = TRUE)
      }

      # save evaluation metrics
      eval_metrics <- data.frame(
        Metric = c("WAIC", "DIC", "MLik", "LCPO (sum log-CPO)"),
        Value = c(fit_glo$waic$waic, fit_glo$dic$dic, fit_glo$mlik[1], lcpo_val_glo)
      )
      write.csv(eval_metrics, file = file.path(values_path, paste0(species, "_evaluation_global.csv")), row.names = FALSE)

      # Save CPO values (one per observation)   #@@@JMB useful for leave-one-out diagnostics or model comparison??
      write.csv(data.frame(CPO = fit_glo$cpo$cpo), file = file.path(values_path, paste0(species, "_pointwise_CPO_global.csv")), row.names = FALSE)

      # save full model object (fit)
      saveRDS(fit_glo, file = file.path(values_path, paste0(species, "_model_fit_global.rds")))

      # save pred
      file_path <- file.path(projections_path, paste0(species, "_Current_global.tif")) 
      terra::writeRaster(terra::unwrap(pred_glo), file_path, overwrite = TRUE)

      # save pred_sp
      if(!is.null(pred_sp_glo)) {
        file_path <- file.path(projections_path, paste0(species, "_Current_Spatial_global.tif"))
        terra::writeRaster(terra::unwrap(pred_sp_glo), file_path, overwrite = TRUE)
      }
    }

    message("Results saved in the following local folder(s):")
    message(paste(
    "  - Current projection (pred.covariate), spatial field (pred_sp), and new scenarios: ", projections_path, "\n",
    " - Fixed/random spatial effects, hyperparameters, evaluation: ", values_path, "\n",
    " - Full model object: ", file.path(values_path, paste0(species, "_model_fit.rds")), "\n",
    if(save.independent) {
      paste0(
      " - Individual projections (global): ", projections_path, "\n",
      "  - Individual model fits, fixed/random effects, hyperparameters, evaluation: ", values_path, "\n")
    }))
  }

  # summary
  separator <- data.frame(Field = "--------------------------------", Value = "---------------------------", stringsAsFactors = FALSE)
  summary_glo <- generate_summary_nsbm(fit_glo, species, spatial, lcpo_val_glo, model = "global")
  if(cv.folds > 1) {
    cv_rows <- data.frame(
      Field = c("CV folds:", "AUC mean ± sd:"),
      Value = c(cv.folds, paste0(round(cv_res_glo$auc_mean,3), " ± ", round(cv_res_glo$auc_sd,3))),
      stringsAsFactors = FALSE
    )
    summary_glo <- rbind(summary_glo, cv_rows)
  }
  summary_cov <- generate_summary_nsbm(fit_cov, species, spatial, lcpo_val_cov, model = "covariate")
  if(cv.folds > 1) {
    cv_rows <- data.frame(
      Field = c("CV folds:", "AUC mean ± sd:"),
      Value = c(cv.folds, paste0(round(cv_res_cov$auc_mean,3), " ± ", round(cv_res_cov$auc_sd,3))),
      stringsAsFactors = FALSE
    )
    summary_cov <- rbind(summary_cov, cv_rows)
  }
  summary_df <- rbind(summary_glo, separator, summary_cov[-1,])

  #
  sabina <- list(
    Species.Name = species,
    args = list(output = output,
       family = family, 
       link = link,
       prior.range = prior.range,
       prior.sigma = prior.sigma,
       spatial = spatial, 
       rm.corr = rm.corr,
       corcut = corcut,
       proj.new.env = proj.new.env,
       cv.folds = cv.folds, 
       n.threads = n.threads, 
       seed = seed,
       save.independent = save.independent
    ),
    Selected.Variables.Global = nsbm_obj$Selected.Variables.Global,
    Selected.Variables.Regional = nsbm_obj$Selected.Variables.Regional,
    Selected.Variables.Covariate <- names(sp_covreg),
    current.projections = current.projections,
    new.projections = if(!is.null(proj_list)) rapply(proj_list, terra::wrap, how = "list") else NULL,
    Summary = summary_df
  )

  class(sabina) <- "nsbm.inlabru"
  return(sabina)
}

