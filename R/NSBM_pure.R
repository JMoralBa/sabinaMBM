#' @name NSBM.pure
#'
#' @title Nested species distribution modeling (pure bayes hierarchical)...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru...
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param family  A standard R \code{family} object (e.g. \code{binomial(link="logit")}), or the character string \code{"cp"} to fit a Cox point process (intensity). 
#' @param spde.mesh An INLA mesh object created externally with `create_mesh()`. includes spatial latent field (SPDE) in the model (default: NULL). 
#' @param spde.pcprior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param spde.pcprior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param nested.intercept Logical; if TRUE = model regional intercept as deviation from global.
#' @param covariate.pcprior.smoothness PC-prior for RW2 smoothing of covariates:
#'   - `NULL` (default): no smoothing, all effects constant
#'   - `list(u, alpha)`: choose your smoothing (e.g., list(u=0.5, alpha=0.01)
#'   - `"auto"`: analyse a few (u, alpha) combinations and take the one with best fit (lowest WAIC)
#' @param proj.new.env Logical. Whether to compute predictions under new scenarios (default: TRUE).
#' @param cv.folds Number of k-folds for cross-validation (default: 1 = no CV). If >1, returns mean ± sd AUC for global and regional models.
#' @param n.threads Number of threads for analysis
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
#' \item{Summary}{\code{data.frame} with key evaluation metrics and significant variables.}  #@@@JMB revisar y refinar
#'
#' @details
#' - binomial(logit): Use for presence–absence (1/0) data. Output: probability.    #@@@JMB revisar todo lo relacionado con family-link (opciones, inverse link in pred, ...)
#' - binomial(cloglog): Use for presence–absence (1/0) data when 1s are very rare (lots of 0s). Output: probability.
#' - poisson(log): Use for non-overdispersed counts (e.g. abundance counts 0,1,2…). Output: expected count
#' - nbinomial(log): Use for overdispersed counts (variance > mean) (e.g. xx). Output: expected count
#' - gaussian(identity): Use for continuous responses normal distribution (e.g. xx). Output: predicted value on original scale.
#' - gaussian(log): Use for strictly positive, skewed continuous data (e.g. xx). Output: predicted value on original scale.
#' - beta(logit): Use for proportions in (0,1) (e.g. rescaled canopy cover). Output: predicted proportion
#' - tweedie(log): Use for semicontinuous data with many zeros and a continuous positive tail (e.g. xx). Output: expected value
#' - cp: (Cox process) Use for presence-only data or spatial point patterns. Output: intensity.
#'
#' @export
NSBM.pure <- function(nsbm_obj, 
                      family = binomial(link = "logit"), # family object binomial(), poisson(), etc., o "cp" para intensity (procesos puntuales)
                      nested.intercept = TRUE,    # TRUE: IRegional se modela como desviación de IGlobal, FALSE interceptos independientes
                      spde.mesh = NULL,           # NULL or mesh objetct para calcular efecto spatial
                      spde.pcprior.range = NULL,
                      spde.pcprior.sigma = NULL,
                      covariate.pcprior.smoothness = NULL,  #NULL, list(u=0.5, alpha=0.01) o "auto". #@@@JMB !!! está en las mismas unidades en las que está construida la mesh. Pensar.....
                      proj.new.env = TRUE,
                      cv.folds = 1,
                      n.threads = 1, 
                      seed = NULL,
                      save.output = FALSE) {

  if(!inherits(nsbm_obj, "nsdm.vinput")) {
    stop("The 'nsbm_obj' must be of class 'nsdm.vinput', Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  if(inherits(family, "family")) {
    fam <- family$family
    lnk <- family$link
  } else if(is.character(family) && length(family) == 1 && family == "cp") {
    fam <- "cp"
    lnk <- NULL
  } else {
    stop("`family` must be either\n",
         "  - a standard family() object (e.g. binomial(link = 'logit'), poisson(link = 'log'), etc.)\n",
         "  - the string 'cp' for a Cox process.")  #@@@JMB poner permitidos o enviar a ?NSBM.pure details?
  }
  valid_links <- list(binomial = c("logit", "cloglog"),
                      poisson = "log",
                      nbinomial = "log",
                      gaussian = c("identity", "log"),
                      beta = "logit",
                      tweedie = "log",
                      cp = NULL)
  if(!fam %in% names(valid_links)) {
    stop("Unsupported family ", fam, ".\n",
         "Please, see ?NSBM.pure details for more.")
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    stop(paste0("Link ", lnk, " is not allowed for family ", fam))
  }
  if(is.null(spde.mesh)) {
    warning("`spde.mesh` is NULL, so the spatial (SPDE) component will be omitted. \n",
            "To include it, create a mesh with `create_mesh()` and pass it to `spde.mesh`.")
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
  # Disable inla’s rw2 min.diff check when PC-prior smoothing is enabled
  if(!is.null(covariate.pcprior.smoothness)) {  #@@@JMB es tramposo quitar este chequeo?
    m <- get("inla.models", inla.get.inlaEnv())
    m$latent$rw2$min.diff <- NULL
    assign("inla.models", m, inla.get.inlaEnv())
  }

  sabina <- list()

  # Data preparation
  sp_covglo <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)

  crs <- sf::st_crs(sp_covglo)

  pp_glo <- rbind(
    cbind(nsbm_obj$SpeciesData.XY.Global,
          resp = if(!is.null(nsbm_obj$Response.Global)) nsbm_obj$Response.Global else 1L),  #@@@JMB nsbm_obj$Response.Global y Regional habría que generarlos en sabinaNSDM input y arrastrar si hay algo
    cbind(nsbm_obj$Background.XY.Global,                                                    # si lo hacemo así poner algún check con stop/warning para que datos y family sean coherentes
          resp = 0L)
  )
  pp_reg <- rbind(
    cbind(nsbm_obj$SpeciesData.XY.Regional,
          resp = if(!is.null(nsbm_obj$Response.Regional)) nsbm_obj$Response.Regional else 1L),  
    cbind(nsbm_obj$Background.XY.Regional,
          resp = 0L)
  )

  pp_glo <- sf::st_as_sf(pp_glo, coords = c("x","y"), crs = crs) %>%  
            sf::st_transform(crs)
  pp_reg <- sf::st_as_sf(pp_reg, coords = c("x","y"), crs = crs) %>% 
            sf::st_transform(crs)

  pp_glo$region <- 1L
  pp_reg$region <- 1L

  # Spatial domain definition
  pts_reg <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs)

  aux <- sf::st_sf(geometry = c(sf::st_geometry(pp_glo), sf::st_geometry(pts_reg)))
  sf::st_crs(aux) <- crs

  bdy_glo <- sf::st_convex_hull(sf::st_union(aux))
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(raster::rasterToPolygons(raster::raster(sp_covreg)))))

  sf::st_crs(bdy_glo) <- crs
  sf::st_crs(bdy_reg) <- crs

  # SPDE
  if(!is.null(spde.mesh)) {
    matern <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = spde.pcprior.range,
      prior.sigma = spde.pcprior.sigma
    )
  }

  # nested intercept
  if(nested.intercept) {
    # IRegional anidado como desviación de IGlobal
    intercept_terms <- c("IGlobal(1)", "IRegional(1, model='iid', group=region)")
  } else {
    # I independientes
    intercept_terms <- c("IGlobal(1)", "IRegional(1)")
  }
  base_intercepts <- paste(intercept_terms, collapse = " + ")

  # auto pc-prior spline rw2 selection
  if(identical(covariate.pcprior.smoothness, "auto")) {
    # grid pairs
    u_vals <- c(0.1, 0.5, 1.0, 2.0, 5.0)   #@@@JMB revisar valores, coherencia y que está en grados? controla curvatura (0.1≈lineal,1≈defecto,5≈muy curva)
    alpha_vals <- c(0.01, 0.05)            #@@@JMB es correcto? regula penalización (0.01 fuerte,0.05 laxa)
    grid <- expand.grid(u = u_vals, alpha = alpha_vals)

    waics <- numeric(nrow(grid))
    for (i in seq_len(nrow(grid))) {
      u_i <- grid$u[i]
      alpha_i <- grid$alpha[i]
      pc_i <- list(u = u_i, alpha = alpha_i)

      # reconstruimos la llamada a NSBM.pure() con pc_i
      mc <- match.call() # captura la llamada original
      mc$covariate.pcprior.smoothness <- pc_i
      mc$save.output <- FALSE  # no guarda archivos intermedios
      mc$cv.folds <- 1   # no hace cv intermedias
      
      waics[i] <- tryCatch({
        fit_i <- eval(mc, parent.frame())
        as.numeric(fit_i$Summary$Value[fit_i$Summary$Field == "WAIC:"])
      }, error = function(e) {
        warning(sprintf("auto PC-prior fail for u=%.2f, α=%.2f: %s", u_i, alpha_i, e$message))
        NA_real_
      })

      # progress
      cat(sprintf(
        "[%d/%d] u=%.2f  α=%.2f  →  WAIC=%s\n",
        i, nrow(grid), u_i, alpha_i,
        if (is.na(waics[i])) "ERROR" else sprintf("%.2f", waics[i])
      ))
    }

    # valid grid pairs
    valid <- !is.na(waics)
    if(!any(valid)) {
      stop("All (u,α) combinations in 'auto' has failed.")
    }

    # select the best u and alpha (min WAIC)
    best_idx <- which.min(waics)
    best_u <- grid$u[best_idx]
    best_alpha <- grid$alpha[best_idx]
    covariate.pcprior.smoothness <- list(u = best_u, alpha = best_alpha)

    message(sprintf(
      " - Best covariate PC-Prior u=%.2f, alpha=%.2f (WAIC=%.1f)",
      best_u, best_alpha, waics[best_idx]
    ))
  }

  # Model components
  cmp_cov <- fcov(obj = nsbm_obj, 
                  spobjglo = "sp_covglo", 
                  spobjreg = "sp_covreg",
                  covariate.pcprior.smoothness = covariate.pcprior.smoothness)
 
  # Formula
  cmp_formula <- if(!is.null(spde.mesh)) {
    paste0("~ ", base_intercepts, " + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ ", base_intercepts, " + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  f_spatial <- if(!is.null(spde.mesh)) " + spatial" else ""
  I_nested <- if(nested.intercept==TRUE) " + IGlobal" else ""

  if(fam == "cp") {
    # for intensity, Cox process (cp)
    pres_glo <- pp_glo[pp_glo$resp != 0L, ]
    pres_reg <- pp_reg[pp_reg$resp != 0L, ]
    # Likelihoods
    lik_glo <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pres_glo,
      samplers = bdy_glo,
      domain = list(geometry = spde.mesh)
    )
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", I_nested, f_spatial, " + ", cmp_cov$like$fregional)),
      data = pres_reg,
      samplers = bdy_reg,
      domain = list(geometry = spde.mesh)
    )
    eta <- paste0("IRegional", I_nested, f_spatial, " + ", cmp_cov$like$fregional)
    pred_formula <- as.formula(paste0("~ exp(",eta,")"))
  } else {
    # for any type of family-object (e.g, binomial(), etc.)
    # Likelihoods
    lik_glo <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pp_glo,
      samplers = bdy_glo,
      domain = list(geometry = spde.mesh),
      control.family = list(link = lnk)
    )
    lik_reg <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IRegional", I_nested, f_spatial, " + ", cmp_cov$like$fregional)), #@@@JMB si I_nested entra, entra el componente spatial de modelo global?
      data = pp_reg,
      samplers = bdy_reg,
      domain = list(geometry = spde.mesh),
      control.family = list(link = lnk)
    )
    eta <- paste0("IRegional", I_nested, f_spatial, " + ", cmp_cov$like$fregional)
    pred_formula <- switch(lnk,
      "logit" = as.formula(paste0("~ 1 / (1 + exp(-(", eta, ")))")),
      "cloglog" = as.formula(paste0("~ 1 - exp(-exp(", eta, "))")), 
      "log" = as.formula(paste0("~ exp(", eta, ")")),
      "identity" = as.formula(paste0("~ ", eta)))
  }
  # k-fold CV
  if(cv.folds > 1) {
    # stratified k folds
    folds_g <- make_stratified_kfolds(pp_glo$resp, cv.folds)
    folds_r <- make_stratified_kfolds(pp_reg$resp, cv.folds)

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
          family = fam,
          formula = lik_glo$formula,
          data = train_g,
          samplers = lik_glo$samplers,
          domain = lik_glo$domain,
          control.family = list(link = lnk)
        ),
        inlabru::like(
          family = fam,
          formula = lik_reg$formula,
          data = train_r,
          samplers = lik_reg$samplers,
          domain = lik_reg$domain,
          control.family = list(link = lnk)
        ),
        options = list(control.compute = list(cpo = FALSE, config = TRUE))
      )

      # rm NAs
      coords_t <- sf::st_coordinates(test_g)
      cov_g <- terra::extract(sp_covglo, coords_t)[, names(sp_covglo), drop = FALSE]
      cov_r <- terra::extract(sp_covreg, coords_t)[, names(sp_covreg), drop = FALSE]
      keep <- complete.cases(cov_g, cov_r)

      # pred and AUC
      pk <- predict(fit_k, test_r, pred_formula)
      aucs[k] <- as.numeric(pROC::auc(test_r$resp, pk$mean))
    }

    cv_res <- list(
      cv.folds = cv.folds, 
      auc_mean = mean(aucs, na.rm=TRUE),
      auc_sd = sd(aucs, na.rm=TRUE)
    )
  } else {
    cv_res <- NULL
  }

  # Model fitting complete
  fit <- inlabru::bru(
    components = cmp,
    lik_glo,
    lik_reg,
    options = list(control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE))
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

  if(!is.null(spde.mesh)) {
    pred_sp <- predict(fit, pred.df, ~ spatial)
    pred_sp <- pred_as_tif(pred_sp, sp_covreg)
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
  spatial <- if(!is.null(spde.mesh)) TRUE else FALSE
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
      family = fam,
      link = lnk,
      spde.mesh = if(!is.null(spde.mesh)) TRUE else FALSE,
      spde.pcprior.range = if(!is.null(spde.mesh)) spde.pcprior.range else NULL,
      spde.pcprior.sigma = if(!is.null(spde.mesh)) spde.pcprior.sigma else NULL,
      nested.intercept = nested.intercept,
      covariate.pcprior.smoothness = covariate.pcprior.smoothness,
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

# formulas
fcov <- function(obj, 
                 spobjglo, 
                 spobjreg,
                 covariate.pcprior.smoothness = covariate.pcprior.smoothness) {
  # vars
  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional

  cmp1 <- paste0(unique(c(vg, vr)), "(1)", collapse = " + ")

  if(is.null(covariate.pcprior.smoothness)) {
    # "const" (default)
    cmpglobal <- paste(sapply(vg, function(X) {
      paste0(X, "GL(main = ", spobjglo, ", main_layer = \"", X, "\", model = \"const\")")
    }), collapse = " + ")

    cmpregional <- paste(sapply(vr, function(X) {
      paste0(X, "RE(main = ", spobjreg, ", main_layer = \"", X, "\", model = \"const\")")
    }), collapse = " + ")
  } else {
    # custom list(u, alpha) or "auto" (best PC-prior RW2 selection for smooth of covariates) 
    u <- covariate.pcprior.smoothness$u     #@@@JMB cambiar de list a c(u, alpha)??
    alpha <- covariate.pcprior.smoothness$alpha

    cmpglobal <- paste(sapply(vg, function(X) {
      paste0(X, "GL(main = ", spobjglo,
        ", main_layer = '", X, "',",
        " model = 'rw2', scale.model = TRUE,",
        " hyper = list(prec = list(prior = 'pc.prec', param = c(", u, ", ", alpha, "))))"
      )
    }), collapse = " + ")

    cmpregional <- paste(sapply(vr, function(X) {
      paste0(
        X, "RE(main = ", spobjreg,
        ", main_layer = '", X, "',",
        " model = 'rw2', scale.model = TRUE,",
        " hyper = list(prec = list(prior = 'pc.prec', param = c(", u, ", ", alpha, "))),",
        " group = region)"
      )
    }), collapse = " + ")
  }

  # Formula likelihood
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
      if(!is.null(spde.mesh)) {
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



