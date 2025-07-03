#' @name NSBM.pure
#'
#' @title Nested species distribution modeling (pure bayes hierarchical)...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru.
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param mesh An INLA mesh object created externally with `create_mesh()`. Required only if \code{spatial = TRUE}. Ignored if \code{spatial = FALSE}.
#' @param output Character, `"intensity"` or `"probability"` (default).
#' @param spatial Logical; include spatial latent field (SPDE) in the model (default: TRUE). 
#' @param prior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param prior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param nested_intercept Logical; if TRUE = model regional intercept as deviation from global.
#' @param auto_select_terms Logical; if TRUE, select for each covariate automatically between “const”, “linear” or “rw2”.  #@@@JMB no funciona bien aun
#' @param family  Character; `"binomial"` (default), supported `"poisson"`, `"nbinomial"`, or `"cp"`.
#' @param link  Character; link function `"logit"` (default for `family = "binomial"`), otherwise `"log"`.
#' @param cv.folds Number of k-folds for cross-validation (default: 1 = no CV). If >1, returns mean ± sd AUC for global and regional models.
#' @param n.threads Number of threads for analysis
#' @param proj.new.env Logical. Whether to compute predictions under new scenarios (default: TRUE).
#' @param seed Optional integer to set the random seed for reproducibility.
#' @param save.output Logical. If TRUE, saves key model outputs (predictions, evaluation, summary ...).
#'
#' @return A named list of class `nsbm.inlabru` with the following elements:
#' \item{Species.Name}{Species name}
#' \item{args}{List of arguments used in the model fitting.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{current.projections}{List with: prediction (`pred`) and spatial field (`pred_sp`).}
#' \item{new.projections}{List of projections to new.env (if `proj.new.env = TRUE`).}
#' \item{Summary}{\code{data.frame} with key evaluation metrics and significant variables.} #@@@JMB revisar y refinar
#'
#' @export
NSBM.pure <- function(nsbm_obj, 
                         output = "probability",
                         family = "binomial",        # #@@@JMB arg util si vamos a aceptar counts
                         link = "logit",             #if(family=="binomial") "logit" else "log", 
                         prior.range = c(5, 0.01),
                         prior.sigma = c(1, 0.01),
                         nested_intercept = TRUE,    # TRUE: IRegional se modela como desviación de IGlobal
                         auto_select_terms = FALSE,
                         spatial = TRUE,
                         mesh = NULL, 
                         n.threads = 1, 
                         cv.folds = 1,
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
  #if(auto_select_terms && (is.null(spline.k))) {
  #  stop("When auto_select_terms = TRUE you must supply both spline.k (e.g. spline.k = 10).")
  #}
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

  sabina <- list()

  # Data preparation
  pres_reg <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Regional, coords = c("x", "y"))
  pres_glo <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  sp_covglo <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)

  crs <- sf::st_crs(sp_covglo)
  sf::st_crs(pres_glo) <- crs
  sf::st_crs(pres_reg) <- crs
  
  pres_glo <- sf::st_transform(pres_glo, crs)
  pres_reg <- sf::st_transform(pres_reg, crs)

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

  pres_reg$region <- 1L  
  pres_reg$presence <- 1L

  pres_glo$region <- 1L
  pres_glo$presence <- 1L  

  abs_reg <- sf::st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
  abs_reg$presence <- 0L
  abs_reg$region <- 1L
  sf::st_crs(abs_reg) <- crs
  abs_reg <- sf::st_transform(abs_reg, crs)

  abs_glo <- sf::st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
  sf::st_crs(abs_glo) <- crs
  abs_glo$presence <- 0L
  abs_glo$region <- 1L
  abs_glo <- sf::st_transform(abs_glo, crs) 

  pp_reg <- rbind(pres_reg, abs_reg)
  pp_glo <- rbind(pres_glo, abs_glo)

  # SPDE
  if(spatial) {
    matern <- INLA::inla.spde2.pcmatern(
      mesh,
      prior.range = prior.range,
      prior.sigma = prior.sigma
    )
  }

  # nested intercept
  if(nested_intercept) {
    # IRegional anidado como desviación de IGlobal
    intercept_terms <- c("IGlobal(1)", "IRegional(1, model='iid', group=region)")
  } else {
    # I independientes
    intercept_terms <- c("IGlobal(1)", "IRegional(1)")
  }
  base_intercepts <- paste(intercept_terms, collapse = " + ")

  # prefilter select terms
  if(auto_select_terms) {
    sm_g <- prefilter_mgcv(pp_glo, sp_covglo, k = 10)  #@@@JMB poner k como argumento??
    sm_r <- prefilter_mgcv(pp_reg, sp_covreg, k = 10)
  } else {
    sm_g <- setNames(rep("const", length(names(sp_covglo))), names(sp_covglo))
    sm_r <- setNames(rep("const", length(names(sp_covreg))), names(sp_covreg))
  }

  # Model components
  cmp_cov <- fcov(obj = nsbm_obj, 
                  spobjglo = "sp_covglo", 
                  spobjreg = "sp_covreg",
                  splinegl = sm_g,
                  splinereg = sm_r)
 
  # Formula
  cmp_formula <- if(spatial) {
    paste0("~ ", base_intercepts, " + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ ", base_intercepts, " + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  if(output == "probability") {
    f_spatial <- if(spatial) " + spatial" else ""
    # Likelihoods
    lik_glo <- inlabru::like(
      family = family,
      formula = as.formula(paste0("presence ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      #formula = as.formula(      #@@@JMB para un futuro uso con datos de abundancia por ejemplo, pero hay que adaptar pres_regional$region <- 1L, global y abs regional y global
      #  paste0(
      #    if(family=="binomial") "presence" else  
      #    if(family%in%c("poisson","nbinomial")) "count" else "geometry",
      #    " ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal
      #  )
      #),
      data = pp_glo,
      samplers = bdy_glo,
      domain = list(geometry = mesh),
      control.family = list(link = link)
    )
    lik_reg <- inlabru::like(
      family = family,
      formula = as.formula(paste0("presence ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      #formula = as.formula(
      #  paste0(
      #    if(family=="binomial") "presence" else  
      #    if(family%in%c("poisson","nbinomial")) "count" else "geometry",
      #    " ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional
      #  )
      #),
      data = pp_reg,
      samplers = bdy_reg,
      domain = list(geometry = mesh),
      control.family = list(link = link)
    )
    pred_formula <- as.formula(paste0("~ 1 / (1 + exp(-(IRegional", f_spatial, " + ", cmp_cov$like$fregional,")))"))
    # Pred formula para dif family/link
    #pred_formula <- switch(family,
    #  binomial = as.formula(paste0("~ 1/(1 + exp(-(IRegional", f_spatial, " + ", cmp_cov$like$fregional, ")))")),
    #  poisson = as.formula(paste0("~ exp(IRegional", f_spatial, " + ", cmp_cov$like$fregional, ")")),
    #  nbinomial = as.formula(paste0("~ exp(IRegional", f_spatial, " + ", cmp_cov$like$fregional, ")")),
    #  stop("Family unavailable")
    #)
  } else {
    lik_glo <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pres_glo,
      samplers = bdy_glo,
      domain = list(geometry = mesh)
    )
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      data = pres_reg,
      samplers = bdy_reg,
      domain = list(geometry = mesh)
    )
    pred_formula <- as.formula(paste0("~ exp(IRegional", f_spatial, " + ", cmp_cov$like$fregional,")"))
  }

  # k-fold CV
  if(cv.folds > 1) {
    # stratified k folds
    folds_g <- make_stratified_kfolds(pp_glo$presence, cv.folds)
    folds_r <- make_stratified_kfolds(pp_reg$presence, cv.folds)

    aucs <- numeric(cv.folds)

    for(k in seq_len(cv.folds)) {
      train_g <- pp_glo[folds_g != k, ]
      test_g <- pp_glo[folds_g == k, ]
      train_r <- pp_reg[folds_r != k, ]
      test_r <- pp_reg[folds_r == k, ]

      # fit fold k
      fit_k <- inlabru::bru(
        components = cmp,
        inlabru::like(
          family = family,
          formula = lik_glo$formula,
          data = train_g,
          samplers = lik_glo$samplers,
          domain = lik_glo$domain,
          control.family = list(link = link)
        ),
        inlabru::like(
          family = family,
          formula = lik_reg$formula,
          data = train_r,
          samplers = lik_reg$samplers,
          domain = lik_reg$domain,
          control.family = list(link = link)
        ),
        options = list(control.compute = list(cpo = FALSE))
      )

      # rm NAs
      coords_t <- sf::st_coordinates(test_g)
      cov_g <- terra::extract(sp_covglo, coords_t)[, names(sp_covglo), drop = FALSE]
      cov_r <- terra::extract(sp_covreg, coords_t)[, names(sp_covreg), drop = FALSE]
      keep <- complete.cases(cov_g, cov_r)

      # pred and AUC
      pk <- predict(fit_k, test_r, pred_formula)
      aucs[k] <- as.numeric(pROC::auc(test_r$presence, pk$mean))
    }

    cv_res <- list(
      cv.folds = cv.folds, 
      auc_mean = mean(aucs, na.rm=TRUE),
      auc_sd = sd(aucs, na.rm=TRUE)
    )
  } else {
    cv_res <- NULL
  }

  # si cv.folds == 1, final fit complete
  # Model fitting
  fit <- inlabru::bru(
    components = cmp,
    lik_glo,
    lik_reg,
    options = list(control.compute = list(cpo = TRUE, waic = TRUE, dic  = TRUE))
  )

  # log sum of conditional predictive ordinates
  lcpo_val <- if(!is.null(fit$cpo$cpo)) round(sum(log(fit$cpo$cpo)), 2) else "Not computed"

  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)
  pred.df$region <- 1L

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
  if(proj.new.env && !is.null(nsbm_obj$Scenarios)) {
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
    # Create directories
    values_path <- file.path("Results", "NSBM_pure", "Values")
    projections_path <- file.path("Results", "NSBM_pure", "Projections")
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

  # summary          #@@@JMB pendiente revisar/completar...
  summary_df <- generate_summary_nsbm(fit, species, spatial, lcpo_val, model = "pure") 
  if(!is.null(cv_res)) {
    cv_rows <- data.frame(
      Field = c("CV folds:", "AUC mean ± sd:"),
      Value = c(as.character(cv_res$cv.folds), sprintf("%.2f ± %.2f", cv_res$auc_mean, cv_res$auc_sd)),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, cv_rows)
  }

  #
  sabina <- list(
    Species.Name = nsbm_obj$SpeciesName,
    args = list(
      family = family,
      link = link,
      prior.range = prior.range,
      prior.sigma = prior.sigma,
      spatial = spatial,
      nested_intercept = nested_intercept,
      proj.new.env = proj.new.env,
      cv.folds = cv.folds,
      n.threads = n.threads,
      seed = seed
    ),
    Selected.Variables.Global = nsbm_obj$Selected.Variables.Global,
    Selected.Variables.Regional = nsbm_obj$Selected.Variables.Regional,
    current.projections = list(
      pred = terra::wrap(pred),
      pred_sp = if(!is.null(pred_sp)) terra::wrap(pred_sp) else NULL
    ),
    new.projections = rapply(proj_list, terra::wrap, how = "list"),
    Summary = summary_df
  )
 
  attr(sabina, "class") <- "nsbm.inlabru"
  return(sabina)

}


### Helps/Auxiliars

# prefilter mgcv term selection
               #@@@JMB no es multivariante, probar funcion gamsel
prefilter_mgcv <- function(pp,       # sf con columnas $presence (0/1) y %region (entero)
                           cov_rast, # SpatRaster
                           k) {      # número de bases spline para cada GAM

  # df without geom
  df <- sf::st_drop_geometry(pp)
  coords <- sf::st_coordinates(pp)
  cov_vals <- terra::extract(cov_rast, coords)[, names(cov_rast), drop = FALSE]
  df <- cbind(df, cov_vals)
  
  # Formula multivariante con splines??
  vars <- names(cov_rast)
  smooth <- paste0("s(", vars, ", k=", k, ", bs='tp')")
  fmla <- as.formula(paste("presence ~", paste(smooth, collapse = " + ")))
  
  # Fit GAM. if select=TRUE penaliza EDF y puede bajarlos a 0???
  gm <- mgcv::gam(
    formula = fmla,
    data = df,
    family = binomial(link = "logit"),
    select = TRUE,
    method = "REML"
  )
  
  # EDF and p-value)
  st <- summary(gm)$s.table
  
  # Decision by covariable
  sel <- setNames(rep("const", length(vars)), vars)
  for (v in vars) {
    rownm <- paste0("s(", v, ")")
    if(!rownm %in% rownames(st)) next
    edf <- st[rownm, "edf"]
    pval <- st[rownm, "p-value"]
    if(pval > 0.05) sel[v] <- "const"
    else if(edf < 1.5) sel[v] <- "linear"
    else sel[v] <- "rw2"
  }
  
  return(sel)
}

# formulas
# original (desuso)
fcov0 <- function(obj, 
                 spobjglo, 
                 spobjreg) {

  vars <- unique(c(obj$Selected.Variables.Global, obj$Selected.Variables.Regional))
  cmp1 <- paste(paste(vars, "(1)", sep = ""), collapse = " + ") # efecto constante

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

# adaptado a auto spline
fcov <- function(obj, 
                 spobjglo,
                 spobjreg,
                 splinegl,
                 splinereg) {

  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional
  cmp1 <- paste0(unique(c(vg, vr)), "(1)", collapse = " + ")
  
  cmpglobal <- paste(sapply(vg, function(X) {
    m <- splinegl[X]
    paste0(X, "GL(main=", spobjglo,
           ", main_layer='", X,
           "', model='", m, "')")
  }), collapse = " + ")
  
  cmpregional <- paste(sapply(vr, function(X) {
    m <- splinereg[X]
    paste0(X, "RE(main=", spobjreg,
           ", main_layer='", X,
           "', model='", m, "'",
           if(m == "rw2") ", group=region" else "",
           ")")
  }), collapse = " + ")
  
  fglobal <- paste(sprintf("%s * %sGL", vg, vg), collapse = " + ")
  fregional <- paste(sprintf("%s * %sRE", vr, vr), collapse = " + ")
  
  list(
    cmp = paste(c(cmp1, cmpglobal, cmpregional), collapse = " + "),
    like = list(fglobal = fglobal, fregional = fregional)
  )
}



###----------------###

# prepare summary
generate_summary_nsbm <- function(fit, species=species, spatial, lcpo_val, model) {

  # Fixed effects and hyperpar
  summary_fixed <- fit$summary.fixed
  hyper <- fit$summary.hyperpar

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
  dic_val <- if(!is.null(fit$dic$dic) && !is.na(fit$dic$dic)) round(fit$dic$dic, 2) else "Not computed"
  waic_val <- if(!is.null(fit$waic$waic) && !is.na(fit$waic$waic)) round(fit$waic$waic, 2) else "Not computed"
  mlik_val <- if(!is.null(fit$mlik) && !is.na(fit$mlik[1,1])) round(fit$mlik[1, 1], 2) else "Not computed"
  
  # 
  spatial_range <- if(!is.null(hyper) && "Range for spatial" %in% rownames(hyper)) {
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
        paste0("NSBM.", model," (with SPDE)")
      } else {
        paste0("NSBM.", model," (no SPDE)")
      },
      dic_val,
      waic_val,
      mlik_val,
      lcpo_val,
      spatial_range,
      if(length(var_labels) > 0) paste(var_labels, collapse = ", ") else "None"
    ),
    stringsAsFactors = FALSE
  )

  return(summary_df)
}


###----------------###

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

###----------------###


make_stratified_kfolds <- function(y, K) {
  # y: vector de 1/0
  # K: número de folds
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  f <- integer(length(y))
  f[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  f[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  f
}



