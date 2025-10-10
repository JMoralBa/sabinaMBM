#' @name NSBM.pure
#'
#' @title Nested species distribution modeling (pure bayes hierarchical)...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru...
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param family A standard R \code{family} object (e.g. \code{binomial(link="logit")}), or the character string \code{"cp"} to fit a Cox point process (intensity). 
#' @param spde.mesh An INLA mesh object created externally with `create_mesh()`. includes spatial latent field (SPDE) in the model (default: NULL). 
#' @param spde.pcprior.range Numeric vector length 2. Pc-prior on spatial range (e.g., `c(5, 0.01)`).
#' @param spde.pcprior.sigma Numeric vector length 2. Pc-prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param latent.pcprior.range Numeric vector of length 2. PC-prior on the spatial range of the latent SPDE for the global covariate (e.g., `c(0.05, 0.05)` in degrees). If `NULL` (default), no latent SPDE is created; if only this is supplied, `latent.pcprior.sigma` defaults to `c(1, 0.01)`.
#' @param latent.pcprior.sigma Numeric vector of length 2. PC-prior on the marginal standard deviation of the latent SPDE for the global covariate (e.g., `c(1, 0.01)`). If `NULL` (default), no latent SPDE is created; if supplied without `latent.pcprior.range`, the range prior defaults to `5 × resolution` of the global raster.
#' @param nested.intercept Logical; if TRUE = model regional intercept as deviation from global.
#' @param covariate.effects Optional named list to control the global/regional covariate effects (see details). If `NULL` (default), no smoothing, all covariate effects remain constant (linear).
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
#' family/link:
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
#' Spatial effects:
#' - Residual spatial field `spatial(geometry, model = matern)`: captures local spatial autocorrelation in the response not explained by covariates. Requires `spde.mesh`.
#' - Global latent field (optional) : if `latent.pcprior.range` and/or `latent.pcprior.sigma` are provided, a single SPDE field `GLspde` represents large-scale global structure and enters the linear predictor via `beta_GL * GLspde`.
#'
#' covariate.effects: control per-covariate effects at global/regional scale.
#' - Option 1 `NULL` (default): no smoothing, All covariates enter linearly (`"const"`) in both scales.
#' - Option 2 named `list`: Provide a list with entries `global`, `regional`, and/or `default`.
#'   Structure:
#'     covariate.effects = list(
#'       global = list(var1 = "const" | "drop" | list(model="rw2", u=..., alpha=...),
#'                     var2 = "const", ...),
#'       regional = list(varA = "const" | "drop" | list(model="rw2", u=..., alpha=...),
#'                       varB = "const", ...),
#'       default = "const" | "drop" | list(model="rw2", u=..., alpha=...)
#'     )
#'     Allowed values per covariate:
#'       - `"const"`: linear effect (we use the term "const" for consistency with inlabru).
#'       - `"drop"`: exclude that covariate.
#'       - `list(model="rw2", u=..., alpha=...)`: rw2 smoothing with a PC-prior on precision. Use the same units as the rasters for `u` and `alpha`.
#'     If a covariate is not listed under `global`/`regional`, it inherits from `default` (recommended `"const"`).
#'     Example:
#'       covariate.effects = list(
#'         global = list(bio4 = "const",
#'                       bio1 = list(model="rw2", u=..., alpha=...),
#'                       bio12 = "drop" ),
#'         regional = list(bio12 = "const",
#'                         bio1 = list(model="rw2", u=..., alpha=...),
#'                         bio4 = "drop"),
#'         default = "const")    # default rule if a covariate is not listed
#'
#' @export
NSBM.pure <- function(nsbm_obj, 
                      family = binomial(link = "logit"), # family object binomial(), poisson(), etc., o "cp" para intensity (procesos puntuales)
                      nested.intercept = TRUE,    # TRUE: IRegional se modela como desviación de IGlobal, FALSE interceptos independientes
                      spde.mesh = NULL,           # NULL or mesh objetct para calcular efecto spatial
                      spde.pcprior.range = NULL,
                      spde.pcprior.sigma = NULL,
                      latent.pcprior.range = NULL,   # activate latent SPDE if not NULL
                      latent.pcprior.sigma = NULL,   # activate latent SPDE if not NULL (plug-in if both NULL)
                      covariate.effects = NULL,
                      proj.new.env = TRUE,
                      cv.folds = 1,
                      n.threads = 1, 
                      seed = NULL,
                      save.output = FALSE) {

  # checks
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
         "Supported families are: ", paste(names(valid_links), collapse = ", "), "\n",
         "Please, see ?NSBM.pure details for more.")
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    stop("Link `", lnk, "` is not allowed for family `", fam, "`.\n",
         "Allowed links for '", fam, "': ", paste(valid_links[[fam]], collapse = ", "), ".\n",
         "Please, see ?NSBM.pure details for more.")
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
  if(!is.null(covariate.effects)) {
    if(!is.list(covariate.effects)) {
      stop("`covariate.effects` must be a list or NULL.")
    }
    allowed_top <- c("global", "regional", "default")
    unknown_top <- setdiff(names(covariate.effects), allowed_top)
    if(length(unknown_top) > 0) {
      stop("`covariate.effects` has invalid top-level entries: ",
         paste(unknown_top, collapse = ", "),
         ". Allowed entries are: 'global', 'regional', 'default'.")
    }
    if(!is.null(covariate.effects$default)) {
      def <- covariate.effects$default
      if(is.list(def)) {
        if(is.null(def$model) || !def$model %in% c("const","rw2","drop")) {
          stop("`covariate.effects$default` must include `model = 'const'|'rw2'|'drop'`.")
        }
        if(identical(def$model, "rw2") && (is.null(def$u) || is.null(def$alpha))) {
          stop("`covariate.effects$default` with `model = 'rw2'` requires both `u` and `alpha`.")
        }
      } else if(is.character(def)) {
        if(!def %in% c("const","drop")) {
          stop("`covariate.effects$default` must be 'const' or 'drop'. To use 'rw2' by default, provide a list with `model = 'rw2'`, `u`, and `alpha`.")
        }
      } else {
        stop("`covariate.effects$default` must be a character or a list.")
      }
    }
  }

  sabina <- list()

  .opt_plus <- function(s) if(!is.null(s) && nzchar(s)) paste0(" + ", s) else ""

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

  pp_glo <- sf::st_as_sf(pp_glo, coords = c("x","y"), crs = crs)
  pp_glo <- sf::st_transform(pp_glo, crs)
  pp_reg <- sf::st_as_sf(pp_reg, coords = c("x","y"), crs = crs) 
  pp_reg <- sf::st_transform(pp_reg, crs)

  pp_glo$region <- 0L
  pp_reg$region <- 1L

  # Spatial domain definition
  pts_reg <- sf::st_as_sf(terra::as.points(sp_covreg, values = FALSE))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs)

  aux <- sf::st_sf(geometry = c(sf::st_geometry(pp_glo), sf::st_geometry(pts_reg)), crs = crs)

  bdy_glo <- sf::st_convex_hull(sf::st_union(aux))
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(terra::as.polygons(!is.na(sp_covreg[[1]]), dissolve = TRUE))))
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

  # latent SPDE for global covariate
  latent_global <- !is.null(latent.pcprior.range) || !is.null(latent.pcprior.sigma)
  if(latent_global) {
    # calc defaults based on global covariate resolution
    rast_gl <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)[[1]]
    res_xy <- terra::res(rast_gl)
    mean_res <- mean(res_xy)
    # default range = 5 × resolution
    if(is.null(latent.pcprior.range)) {
      latent.pcprior.range <- c(mean_res * 5, 0.05)
      message(sprintf("latent.pcprior.range by default: c(%.5f, 0.05)", mean_res*5))
    }
    # default sigma
    if(is.null(latent.pcprior.sigma)) {
      latent.pcprior.sigma <- c(1, 0.01)
      message("latent.pcprior.sigma by default: c(1, 0.01)")
    }
    spde_cov <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = latent.pcprior.range,
      prior.sigma = latent.pcprior.sigma)
  } else {
    spde_cov <- NULL
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

  # Model components
  cmp_cov <- fcov(obj = nsbm_obj, 
                  spobjglo = "sp_covglo", 
                  spobjreg = "sp_covreg",
                  covariate.effects = covariate.effects,
                  use_latent = latent_global,
                  spde_cov = spde_cov,
                  sp_covglo = sp_covglo,
                  sp_covreg = sp_covreg,
                  pp_glo_sf = pp_glo,
                  pp_reg_sf = pp_reg)
 
  # Formula
  cmp_formula <- if(!is.null(spde.mesh)) {
    paste0("~ ", base_intercepts, " + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ ", base_intercepts, " + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  f_spatial <- if(!is.null(spde.mesh)) " + spatial" else ""
  I_nested <- if(nested.intercept==TRUE) " + IGlobal" else ""
  extra_latent <- if(latent_global) " + beta_GL * GLspde" else ""

  if(fam == "cp") {
    # for intensity, Cox process (cp)
    pres_glo <- pp_glo[pp_glo$resp != 0L, ]
    pres_reg <- pp_reg[pp_reg$resp != 0L, ]
    # Likelihoods
    lik_glo <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, extra_latent, .opt_plus(cmp_cov$like$fglobal))),
      data = pres_glo,
      samplers = bdy_glo,
      domain = list(geometry = spde.mesh)
    )
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", I_nested, f_spatial, extra_latent, .opt_plus(cmp_cov$like$fregional))),
      data = pres_reg,
      samplers = bdy_reg,
      domain = list(geometry = spde.mesh)
    )
    eta <- paste0("IRegional", I_nested, f_spatial, extra_latent, .opt_plus(cmp_cov$like$fregional))
    pred_formula <- as.formula(paste0("~ exp(",eta,")"))
  } else {
    # for any type of family-object (e.g, binomial(), etc.)
    # Likelihoods
    lik_glo <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IGlobal", f_spatial, extra_latent, .opt_plus(cmp_cov$like$fglobal))),
      data = pp_glo,
      samplers = bdy_glo,
      domain = list(geometry = spde.mesh),
      control.family = list(link = lnk)
    )
    lik_reg <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IRegional", I_nested, f_spatial, extra_latent, .opt_plus(cmp_cov$like$fregional))), #@@@JMB si I_nested entra, entra el componente spatial de modelo global?
      data = pp_reg,
      samplers = bdy_reg,
      domain = list(geometry = spde.mesh),
      control.family = list(link = lnk)
    )
    eta <- paste0("IRegional", I_nested, f_spatial, extra_latent, .opt_plus(cmp_cov$like$fregional))
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
        options = list(control.compute = list(cpo = FALSE, config = TRUE),
                       control.inla = list(int.strategy = "eb"),
                       control.mode = list(restart = TRUE))
      )

      # rm NAs
      coords_t <- sf::st_coordinates(test_r)
      both <- c(sp_covglo, sp_covreg)
      cov_all <- terra::extract(both, coords_t)
      keep <- stats::complete.cases(cov_g, cov_r)

      # pred and AUC
      test_r2 <- test_r[keep, ]
      if(nrow(test_r2) > 0) {
        pk <- predict(fit_k, test_r2, pred_formula)
        aucs[k] <- as.numeric(pROC::auc(test_r2$resp, pk$mean))
      } else {
        aucs[k] <- NA_real_
      }
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
    options = list(control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
                   control.inla = list(int.strategy = "eb"),
                   control.mode = list(restart = TRUE))
  )

  lcpo_val <- if(!is.null(fit$cpo$cpo)) round(sum(log(fit$cpo$cpo)), 2) else NA_real_

  # Predictions
  pred.df <- sf::st_as_sf(terra::as.points(sp_covreg, values = FALSE))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)
  pred.df$region <- 1L

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
    for(sc in names(nsbm_obj$Scenarios)) {
      #sc <- 1
      scen_rast <- terra::unwrap(nsbm_obj$Scenarios[[sc]])
      scen_df <- sf::st_as_sf(terra::as.points(scen_rast, values = FALSE))
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
    # save CPO values (one per observation)   #@@@JMB useful for leave-one-out diagnostics or model comparison??
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
    Species.Name = species,
    args = list(
      family = fam,
      link = lnk,
      spde.mesh = !is.null(spde.mesh),
      spde.pcprior.range = if(!is.null(spde.mesh)) spde.pcprior.range else NULL,
      spde.pcprior.sigma = if(!is.null(spde.mesh)) spde.pcprior.sigma else NULL,
      latent.pcprior.range = latent.pcprior.range,
      latent.pcprior.sigma = latent.pcprior.sigma,
      nested.intercept = nested.intercept,
      covariate.effects = covariate.effects,
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
    new.projections = if(length(proj_list) > 0) lapply(proj_list, terra::wrap) else list(),
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
                 use_latent = FALSE,
                 spde_cov = NULL,
                 sp_covglo,
                 sp_covreg,
                 covariate.effects = NULL,
                 pp_glo_sf = NULL,
                 pp_reg_sf= NULL) {

  # selected vars
  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional

  cmp1 <- paste0(unique(c(vg, vr)), "(1)", collapse = " + ")

  # resolve per-covariate spec (covariate.effects)
  .spec_cache <- new.env(parent = emptyenv())
  resolve_spec <- function(varname, scale, ccform, default_model = "const") {
    key <- paste(scale, varname, sep = "||")
    hit <- .spec_cache[[key]]
    if(!is.null(hit)) return(hit)
    if(is.null(ccform) || !is.list(ccform)) {
      res <- list(model = default_model, u = NA, alpha = NA)
     .spec_cache[[key]] <- res
      return(res)
    }
    path <- paste0("covariate.effects$", scale, "$", varname)
    spec <- if(!is.null(ccform[[scale]]) && !is.null(ccform[[scale]][[varname]])) {
      ccform[[scale]][[varname]]
    } else if(!is.null(ccform$default)) {
      ccform$default
    } else default_model

    if(is.character(spec)) {
      if(identical(spec, "const")) { res <- list(model = "const", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(spec, "drop")) { res <- list(model = "drop", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(spec, "rw2")) {
        stop(paste0("`", path, "` uses model='rw2' but is missing `u` and/or `alpha`",
             "Use `list(model='rw2', u=..., alpha=...)`."))
      }
      stop(paste0("`", path, "` has invalid value '", spec, "'. ",
                  "Use 'const'|'drop' or a list with `model='rw2'`, `u`, and `alpha`."))
    } else if(is.list(spec)) {
      m <- spec$model
      if(is.null(m)) stop(path, " is missing `model`.")
      if(identical(m, "const")) { res <- list(model = "const", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(m, "drop")) { res <- list(model = "drop", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(m, "rw2")) {
        if(is.null(spec$u) || is.null(spec$alpha)) {
          stop(paste0("`", path, "` with `model = 'rw2'` requires both `u` and `alpha`."))
        }
        res <- list(model = "rw2", u = spec$u, alpha = spec$alpha); .spec_cache[[key]] <- res; return(res)
      }
      stop(paste0("`", path, "` has invalid `model = '", m, "'`. Use 'const'|'rw2'|'drop'."))
    } else {
      stop(paste0("`", path, "` has an unsupported specification type."))
    }
  }

  # control "drop"     #@@@JMB subir a checks de NSBM_pure()
  if(is.list(covariate.effects)) {
    all_dropped_gl <- length(obj$Selected.Variables.Global) > 0 && all(vapply(
      obj$Selected.Variables.Global,
      function(X) resolve_spec(X, "global", covariate.effects, "const")$model == "drop",
      logical(1)))
    all_dropped_re <- length(obj$Selected.Variables.Regional) > 0 && all(vapply(
      obj$Selected.Variables.Regional,
      function(X) resolve_spec(X, "regional", covariate.effects, "const")$model == "drop",
      logical(1)))
    if(all_dropped_gl && all_dropped_re && is.null(spde_cov) && !use_latent) {
      stop("covariate.effects: all covariates dropped in both scales and no spatial/latent field; model would be intercept-only.")
    }
    if(all_dropped_gl && !use_latent) {
      warning("covariate.effects: all global covariates are dropped; only regional (and spatial/latent if present) will contribute.")
    }
    if(all_dropped_re) {
      warning("covariate.effects: all regional covariates are dropped; only global (and spatial/latent if present) will contribute.")
    }
  }

  # rw2: thin knots to enforce min relative spacing (INLA check 1e-3) 
  thin_knots <- function(x, min_ratio = 1e-3) {
    x <- sort(unique(as.numeric(x)))
    if(length(x) <= 2) return(x)
    r <- diff(range(x))
    if(!is.finite(r) || r == 0) return(unique(x))
    out <- x[1]
    for(xi in x[-1]) {
      if((xi - out[length(out)]) / r >= min_ratio) out <- c(out, xi)
    }
    if(length(out) < 3L) {
      out <- seq(min(x), max(x), length.out = min(max(5L, length(x)), 50L))
    }
    unique(out)
  }

  #rw2 defaults (inla recommends quantile grouping and K = 150-300 (balance stability and flexibility)
  rw2_K <- 300L              #@@@JMB pensar si dejamos esto por defecto
  rw2_method <- "quantile" 

  # build rw2 support via inla.group + thinning
  build_rw2_values <- function(rast_layer, coords, K = rw2_K, method = rw2_method) {
    vals <- suppressWarnings(as.numeric(terra::extract(rast_layer, coords)[, 1]))
    vals <- vals[is.finite(vals)] 
    rng <- range(vals)
    if(diff(rng) == 0) {
      stop("rw2: covariate has zero range; rw2 requires variability.")
    }    
    g <- INLA::inla.group(vals, n = K, method = method)
    v <- sort(unique(as.numeric(levels(g))))
    v <- thin_knots(v, min_ratio = 1e-3)
    if(length(v) < 3L) {
      stop("RW2: <3 knots after grouping/thinning.") #@@@JMB aquí también se podría ajustar K, pero demasiados args en mi opinión
    }
    v
  }

  fmt_vec <- function(v) paste0("c(", paste(format(v, digits = 7), collapse = ","), ")")

  # detect if rw2 is request
  rw2_effects <- is.list(covariate.effects) && !is.null(covariate.effects)
  need_rw2_global <- rw2_effects && length(vg) > 0 &&
    any(vapply(vg, function(X) resolve_spec(X, "global",  covariate.effects, "const")$model == "rw2", logical(1)))
  need_rw2_regional <- rw2_effects && length(vr) > 0 &&
    any(vapply(vr, function(X) resolve_spec(X, "regional", covariate.effects, "const")$model == "rw2", logical(1)))
 
  coords_all <- NULL
  coords_r <- NULL
  if(need_rw2_global) {
    coords_all <- rbind(sf::st_coordinates(pp_glo_sf),
                        sf::st_coordinates(pp_reg_sf))
  }
  if(need_rw2_regional) {
    coords_r <- sf::st_coordinates(pp_reg_sf)
  }

  # no covariate.effects (all const)
  if(is.null(covariate.effects)) {
    cmpglobal <- if(length(vg) > 0) paste0(
      paste0(vg, "GL(main = ", spobjglo, ", main_layer = '", vg, "', model = 'const')"),
      collapse = " + ") else ""
    fglobal <- if(length(vg) > 0) paste0(vg, " * ", vg, "GL", collapse = " + ") else ""
    cmpregional <- if(length(vr) > 0) paste0(
      paste0(vr, "RE(main = ", spobjreg, ", main_layer = '", vr, "', model = 'const')"),
      collapse = " + ") else ""
    fregional <- if (length(vr) > 0) paste0(vr, " * ", vr, "RE", collapse = " + ") else ""
    return(list(
      cmp = paste(Filter(nzchar, c(cmp1, cmpglobal, cmpregional)), collapse = " + "),
      like = list(fglobal = if(nzchar(fglobal)) fglobal else "",
                  fregional = if(nzchar(fregional)) fregional else "")
    ))
  }

  # global
  cmpglobal <- ""
  fglobal <- ""

  if(length(vg) > 0) {
    if(use_latent && !is.null(spde_cov)) {
      # With latent SPDE (shared global spde)
      # if any global cov requests rw2, add rw2 over global spde. Else GLspde + beta_GL
      rw2_pairs <- list()
      if(rw2_effects) {
        for(X in vg) {
          specX <- resolve_spec(X, "global", covariate.effects, "const")
          if(identical(specX$model, "rw2")) {
            rw2_pairs[[length(rw2_pairs) + 1]] <- c(specX$u, specX$alpha)
          }
        }
      }
      covblk <- "GLspde(main = geometry, model = spde_cov)"
      if(length(rw2_pairs) == 0) {
        cmpglobal <- paste(covblk, "beta_GL(1)", sep = " + ")
        fglobal <- ""
      } else {
        ua <- do.call(rbind, rw2_pairs)
        ua_uniq <- unique(ua)
        if(nrow(ua_uniq) > 1) {
          stop("Whit latent global SPDE, all global RW2 must sharee (u, alpha). Provide the same pair for all global covariates (e.g., via `covariate.effects$global` or `covariate.effects$default`).")
        }
        u_g <- ua_uniq[1, 1]
        a_g <- ua_uniq[1, 2]
        rw2blk <- paste0(
          "s(GLspde, model='rw2', scale.model=TRUE, ",
          "hyper=list(prec=list(prior='pc.prec',param=c(", u_g, ",", a_g, "))))"
        )
        cmpglobal <- paste(covblk, rw2blk, sep = " + ")
        fglobal <- ""  # with latent no X * XGL terms
      }
    } else {
      # Without latent (per-covariate GL components)
      gl_terms <- character(0)
      fgl_terms <- character(0)

      for(X in vg) {
        specX <- resolve_spec(X, "global", covariate.effects, "const")
        if(identical(specX$model, "drop")) next

        if(identical(specX$model, "const")) {
          gl_terms <- c(gl_terms, paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'const')"))
        } else if(identical(specX$model, "rw2")) {
          vals <- if(!is.null(coords_all)) build_rw2_values(sp_covglo[[X]], coords_all) else numeric(0)
          core <- paste0(
            X, "GL(main = ", spobjglo, ", main_layer = '", X, "', ",
            "model = 'rw2', scale.model = TRUE, ",
            "hyper=list(prec=list(prior='pc.prec',param=c(", specX$u, ",", specX$alpha, ")))",
            if(length(vals) >= 3L) paste0(", values = ", fmt_vec(vals)) else "",
            ")"
          )
          gl_terms <- c(gl_terms, core)
        }
        fgl_terms <- c(fgl_terms, paste0(X, " * ", X, "GL"))
      }

      if(length(gl_terms) > 0) cmpglobal <- paste(gl_terms, collapse = " + ")
      if(length(fgl_terms) > 0) fglobal <- paste(fgl_terms, collapse = " + ")
    }
  }

  # regional
  cmpregional <- ""
  fregional <- ""

  if(length(vr) > 0) {
    re_terms <- character(0)
    fre_terms <- character(0)

    for(X in vr) {
      specX <- resolve_spec(X, "regional", covariate.effects, "const")
      if(identical(specX$model, "drop")) next

      if(identical(specX$model, "const")) {
        re_terms <- c(re_terms, paste0(X, "RE(main = ", spobjreg, ", main_layer = '", X, "', model = 'const')"))
      } else if(identical(specX$model, "rw2")) {
        vals <- if(!is.null(coords_r)) build_rw2_values(sp_covreg[[X]], coords_r) else numeric(0)
        core <- paste0(
          X, "RE(main = ", spobjreg, ", main_layer = '", X, "', ",
          "model = 'rw2', scale.model = TRUE, ",
          "hyper = list(prec = list(prior = 'pc.prec', param = c(", specX$u, ",", specX$alpha, "))), ",
          "group = region",
          if(length(vals) >= 3L) paste0(", values = ", fmt_vec(vals)) else "",
          ")"
        )
        re_terms <- c(re_terms, core)
      }
      fre_terms <- c(fre_terms, paste0(X, " * ", X, "RE"))
    }
    if(length(re_terms) > 0) cmpregional <- paste(re_terms, collapse = " + ")
    if(length(fre_terms) > 0) fregional <- paste(fre_terms, collapse = " + ")
  }

  # output fcov
  list(
    cmp = paste(Filter(nzchar, c(cmp1, cmpglobal, cmpregional)), collapse = " + "),
    like = list(fglobal = if(nzchar(fglobal)) fglobal else "",
                fregional = if(nzchar(fregional)) fregional else "")
  ) 

}


###----------------###




