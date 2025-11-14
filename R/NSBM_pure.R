#' @name NSBM.pure
#'
#' @title Nested species distribution modeling (pure bayes hierarchical)...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru...
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param family A standard R \code{family} object (e.g. \code{binomial(link="logit")}), or the character string \code{"cp"} to fit a Cox point process (intensity). 
#' @param spde.mesh An INLA mesh object created externally with `create_mesh()`. Requiered if spatial and/or latent SPDE components are used. If `NULL` (default), the model runs without spatial structure. 
#' @param spde.pcprior.range Numeric vector length 2. Pc-prior on spatial range (e.g., `c(5, 0.01)`).
#' @param spde.pcprior.sigma Numeric vector length 2. Pc-prior on marginal standard deviation (σ) (e.g., `c(1, 0.01)`).
#' @param latent.pcprior.range Numeric vector of length 2. PC-prior on the range of the latent global field SPDE for the global covariate (e.g., `c(0.05, 0.05)` in degrees). If `NULL` (default), no latent SPDE is created.
#' @param latent.pcprior.sigma Numeric vector of length 2. PC-prior on the marginal standard deviation (σ) of the latent global field SPDE for the global covariate (e.g., `c(1, 0.01)`). If `NULL` (default), no latent SPDE is created.
#' @param nested.intercept Logical; if TRUE = model regional intercept as deviation from global.
#' @param covariate.effects Optional named list to control the global/regional covariate effects (see details). If `NULL` (default), no smoothing, all covariate effects remain constant (linear).
#' @param proj.new.env Logical. Whether to compute predictions under new scenarios (default: TRUE).
#' @param cv.folds Number of k-folds for cross-validation (default: 1 = no CV). If >1, returns mean ± sd AUC for global and regional models.
#' @param n.threads Number of threads used by INLA (default = 1).
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
#' \item{Summary}{\code{data.frame} summarizing model fit, diagnostics, and significant variables.}  #@@@JMB revisar y refinar
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
#' SPDE structure and priors:
#' - The mesh (`spde.mesh`) defines the domain for spatial and latent fields.  
#' - If neither spatial nor latent priors are defined, the model runs without SPDE components.
#' - Residual spatial field `spatial(geometry, model=matern)`: captures local spatial autocorrelation in the response not explained by covariates. Requires `spde.mesh`.
#' - Latent global field:` captures broad-scale structure  across regions and enters the linear predictor via `beta_GL * GLspde`.
#' - If both spatial and latent fields are used, the latent range should be ≥ 3× the spatial range to avoid overlap (Bakka et al., 2018).
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

  spatial_local <- !is.null(spde.mesh) && (!is.null(spde.pcprior.range) || !is.null(spde.pcprior.sigma))
  latent_global <- !is.null(spde.mesh) && (!is.null(latent.pcprior.range) || !is.null(latent.pcprior.sigma))

  # checks
  if(!inherits(nsbm_obj, "nsdm.vinput")) {
    stop("❌ The 'nsbm_obj' must be of class 'nsdm.vinput'. Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  if(inherits(family, "family")) {
    fam <- family$family
    lnk <- family$link
  } else if(is.character(family) && length(family) == 1 && family == "cp") {
    fam <- "cp"
    lnk <- NULL
  } else {
    stop("❌ `family` must be either:\n",
         "  - a standard family() object (e.g. binomial(link = 'logit'), poisson(link = 'log'), etc.)\n",
         "  - the string 'cp' for a Cox process.\n")  #@@@JMB poner permitidos o enviar a ?NSBM.pure details?
  }
  valid_links <- list(binomial = c("logit", "cloglog"),
                      poisson = "log",
                      nbinomial = "log",
                      gaussian = c("identity", "log"),
                      beta = "logit",
                      tweedie = "log",
                      cp = NULL)
  if(!fam %in% names(valid_links)) {
    stop("❌ Unsupported family ", fam, ".\n",
         "  Supported families are: ", paste(names(valid_links), collapse = ", "), "\n",
         "  Please, see ?NSBM.pure details for more.\n")
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    stop("❌ Link `", lnk, "` is not allowed for family `", fam, "`.\n",
         "  Allowed links for '", fam, "': ", paste(valid_links[[fam]], collapse = ", "), ".\n",
         "  Please, see ?NSBM.pure details for more.\n")
  }
  if(is.null(spde.mesh)) {
    warning("⚠️ No `spde.mesh` provided, so the spatial and latent SPDE components will be omitted. \n",
            "  To include them, create a mesh with `create_mesh()` and pass it to `spde.mesh`.\n")
  } else if(!inherits(spde.mesh, "inla.mesh")) {
    stop("❌ `spde.mesh` must be a valid INLA mesh object (class 'inla.mesh'). Create a mesh with `create_mesh()`.")
  }
  if(spatial_local) {
    if(is.null(spde.pcprior.range) || is.null(spde.pcprior.sigma)) {
      stop("❌ Missing priors for spatial field. Define both `spde.pcprior.range` and `spde.pcprior.sigma`.")
    }
  }
  if(latent_global) {
    if(is.null(latent.pcprior.range) || is.null(latent.pcprior.sigma)) {
      stop("❌ Missing priors for latent field. Define both `latent.pcprior.range` and `latent.pcprior.sigma`.")
    }
  }
  if(spatial_local && latent_global) {
    if(latent.pcprior.range[1] < spde.pcprior.range[1] * 3) {
      warning("⚠️ `latent.pcprior.range` < 3× `spde.pcprior.range`: fields may overlap, causing double-counting of spatial variance.")
    }
  }
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1) {
      stop("❌ 'seed' must be a single numeric value.\n")
    }
    set.seed(seed)
  }
  available_cores <- parallel::detectCores(logical = TRUE)
  if(!is.null(n.threads) && n.threads > available_cores) {
    stop(paste0("❌ Requested `n.threads` = ", n.threads, " exceeds available cores (", available_cores,").\n"))
  }
  INLA::inla.setOption(num.threads = n.threads)
  if(!is.null(covariate.effects)) {
    if(!is.list(covariate.effects)) {
      stop("❌ `covariate.effects` must be a list or 'NULL'.\n")
    }
    allowed_top <- c("global", "regional", "default")
    unknown_top <- setdiff(names(covariate.effects), allowed_top)
    if(length(unknown_top) > 0) {
      stop("❌ Invalid entries in `covariate.effects`: ",
         paste(unknown_top, collapse = ", "),
         ". Allowed: 'global', 'regional', 'default'.\n")
    }
    if(!is.null(covariate.effects$default)) {
      def <- covariate.effects$default
      if(is.list(def)) {
        if(is.null(def$model) || !def$model %in% c("const","rw2","drop")) {
          stop("❌ `covariate.effects$default` must include `model = 'const'|'rw2'|'drop'`.\n")
        }
        if(identical(def$model, "rw2") && (is.null(def$u) || is.null(def$alpha))) {
          stop("❌ 'rw2' model on `covariate.effects$default` requires both `u` and `alpha` parameters.\n")
        }
      } else if(is.character(def)) {
        if(!def %in% c("const","drop")) {
          stop("❌ `covariate.effects$default` must be 'const' or 'drop'.\n")
        }
      } else {
        stop("❌ `covariate.effects$default` must be a character or list.\n")
      }
    }
    finding_drops <- function(x, scale) {
      if(!is.null(covariate.effects[[scale]][[x]])) {
        covariate.effects[[scale]][[x]]
      } else if(!is.null(covariate.effects$default)) {
        covariate.effects$default
      } else "const"
    }
    all_dropped_gl <- length(nsbm_obj$Selected.Variables.Global) > 0 &&
      all(vapply(nsbm_obj$Selected.Variables.Global,
                 function(x) finding_drops(x, "global") == "drop",
                 logical(1)))
    all_dropped_re <- length(nsbm_obj$Selected.Variables.Regional) > 0 &&
      all(vapply(nsbm_obj$Selected.Variables.Regional,
                 function(x) finding_drops(x, "regional") == "drop",
                 logical(1)))
    if(all_dropped_gl && all_dropped_re) {
      message("ℹ️ All covariates dropped at both global and regional scales.\n")
    } else if(all_dropped_gl) {
      message("ℹ️ All global covariates dropped.\n")
    } else if(all_dropped_re) {
      message("ℹ️ All regional covariates dropped.\n")
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


  # SPDE domain definition
  pts_reg <- sf::st_as_sf(terra::as.points(sp_covreg, values = FALSE))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs)

  aux <- sf::st_sf(geometry = c(sf::st_geometry(pp_glo), sf::st_geometry(pts_reg)), crs = crs)

  bdy_glo <- sf::st_convex_hull(sf::st_union(aux))
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(terra::as.polygons(!is.na(sp_covreg[[1]]), dissolve = TRUE))))
  sf::st_crs(bdy_glo) <- crs
  sf::st_crs(bdy_reg) <- crs


  # spatial local
  if(spatial_local)  {
    matern <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = spde.pcprior.range,
      prior.sigma = spde.pcprior.sigma
    )
  }


  # Latent global 
  if(latent_global) {
    spde_cov <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = latent.pcprior.range,
      prior.sigma = latent.pcprior.sigma)
  } else {
    spde_cov <- NULL
  }


  # Nested intercept
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

  if(latent_global && any(sapply(covariate.effects$global, function(x) x$model == "rw2"))) {
    warning("⚠️ Both a latent field and global rw2 smoothers are active. This may double-count spatial smoothing and dilute uncertainty. Consider disabling rw2 for global covariates.")
  }
 

  # Formula
  cmp_formula <- if(spatial_local) {
    paste0("~ ", base_intercepts, " + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ ", base_intercepts, " + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  dom <- if(spatial_local | latent_global) list(geometry = spde.mesh) else NULL
  f_spatial <- if(spatial_local) " + spatial" else ""
  I_nested <- if(nested.intercept==TRUE) " + IGlobal" else ""
  f_latent <- if(latent_global) " + beta_GL * GLspde" else ""

  if(fam == "cp") {
    # for intensity, Cox process (cp)
    pres_glo <- pp_glo[pp_glo$resp != 0L, ]
    pres_reg <- pp_reg[pp_reg$resp != 0L, ]
    # Likelihoods
    lik_glo <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, f_latent, .opt_plus(cmp_cov$like$fglobal))),
      data = pres_glo,
      samplers = bdy_glo,
      domain = dom
    )
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", I_nested, f_spatial, f_latent, .opt_plus(cmp_cov$like$fregional))),
      data = pres_reg,
      samplers = bdy_reg,
      domain = dom
    )
    eta <- paste0("IRegional", I_nested, f_spatial, f_latent, .opt_plus(cmp_cov$like$fregional))
    pred_formula <- as.formula(paste0("~ exp(",eta,")"))
  } else {
    # for any type of family-object (e.g, binomial(), etc.)
    # Likelihoods
    lik_glo <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IGlobal", f_spatial, f_latent, .opt_plus(cmp_cov$like$fglobal))),
      data = pp_glo,
      samplers = bdy_glo,
      domain = dom,
      control.family = list(link = lnk)
    )
    lik_reg <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ IRegional", I_nested, f_spatial, f_latent, .opt_plus(cmp_cov$like$fregional))), #@@@JMB si I_nested entra, entra el componente spatial de modelo global?
      data = pp_reg,
      samplers = bdy_reg,
      domain = dom,
      control.family = list(link = lnk)
    )
    eta <- paste0("IRegional", I_nested, f_spatial, f_latent, .opt_plus(cmp_cov$like$fregional))
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
      keep <- stats::complete.cases(cov_all)

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


  # latent SPDE for global covariate diagnostics/warnings
  if(latent_global) {
    hyp <- fit$summary.hyperpar
    # hyperpars
    i_lat_range <- which(grepl("Range.*GLspde", rownames(hyp), ignore.case = TRUE))
    i_lat_sigma <- which(grepl("Stdev.*GLspde|Sigma.*GLspde", rownames(hyp), ignore.case = TRUE))
    i_sp_range <- grep("Range.*(spatial|matern)", rownames(hyp), ignore.case = TRUE)
    i_sp_sigma <- grep("(Stdev|Sigma).*(spatial|matern)", rownames(hyp), ignore.case = TRUE)

    # posterior mean
    lat_range <- if(length(i_lat_range)) hyp$mean[i_lat_range] else NA_real_
    lat_sigma <- if(length(i_lat_sigma)) hyp$mean[i_lat_sigma] else NA_real_
    sp_range <- if(length(i_sp_range)) hyp$mean[i_sp_range] else NA_real_
    sp_sigma <- if(length(i_sp_sigma)) hyp$mean[i_sp_sigma] else NA_real_

    # spatial and latent overlap?
    if(is.finite(lat_range) && is.finite(sp_range)) {
      ratio <- lat_range / sp_range
      if(ratio < 1) {
        warning(paste0(
          "⚠️ Potential overlap between latent GLspde and residual SPDE (range ratio = ", round(ratio, 2), ").\n",
          "   Latent and residual fields may be capturing similar spatial scales. This can cause identifiability issues or double-smoothing.\n",
          "   Consider increasing latent.pcprior.range to at least 3–5× the residual range prior, or reduce spatial.pcprior.range to ensure clear scale separation."
        ), call. = FALSE)
      }
    }
  } 


  # Predictions
  pred.df <- sf::st_as_sf(terra::as.points(sp_covreg, values = FALSE))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)
  pred.df$region <- 1L

  pred <- predict(fit, pred.df, pred_formula)
  pred <- pred_as_tif(pred, sp_covreg) # from sf to tif

  if(spatial_local) {
    pred_sp <- predict(fit, pred.df, ~ spatial)
    pred_sp <- pred_as_tif(pred_sp, sp_covreg)
  } else {
    pred_sp <- NULL
  }

  if(latent_global) {
    pred_lat <- predict(fit, pred.df, ~ beta_GL * GLspde)
    pred_lat  <- pred_as_tif(pred_lat, sp_covreg)
  } else {
    pred_lat <- NULL
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
      scen_df$region <- 1L
      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[sc]] <- pred_as_tif(proj_pred, template = scen_rast)  # from sf to tif
    }
  }


  # diagnostics
  coords_reg <- sf::st_coordinates(pp_reg)
  data_used  <- data.frame(x = coords_reg[,1], y = coords_reg[,2], resp = pp_reg$resp)
  diag_block <- nsbm_diagnostics(
                  fit, 
                  data_used, 
                  priors = list(spde.pcprior.range = spde.pcprior.range,
                                spde.pcprior.sigma = spde.pcprior.sigma,
                                latent.pcprior.range = latent.pcprior.range,
                                latent.pcprior.sigma = latent.pcprior.sigma),
                  pred_sp = if(spatial_local) pred_sp else NULL,
                  pred_lat = if(latent_global) pred_lat else NULL)
  if(length(diag_block$warnings)) {
    message(paste(unique(diag_block$warnings), collapse = " | "))
  }


  # save outputs
  species <- nsbm_obj$Species.Name
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
    # save latent contribution
    if(!is.null(pred_lat)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_latent.tif"))
      terra::writeRaster(terra::unwrap(pred_lat), file_path, overwrite = TRUE)
    }
    # save new scenarios
    if(length(proj_list) > 0 && !is.null(nsbm_obj$Scenarios)) {
      for(i in seq_along(proj_list)) {
        sc_name <- names(proj_list)[i]
        file_path <- file.path(projections_path, paste0(species, "_", sc_name, ".tif"))
        terra::writeRaster(terra::unwrap(proj_list[[i]]), file_path, overwrite = TRUE)
      }
    }
    # save diagnostics if spdf field
    if(spatial_local || latent_global) {
      file_path <- file.path(values_path, paste0(species, "_diagnostics.png"))
      ggplot2::ggsave(file_path, diag_block$composite, width = 8, height = 12, dpi = 300)
    }

    message("Results saved in the following local folder(s):")
    message(paste(
    "  - Current projection (pred), spatial field (pred_sp), and new scenarios: ", projections_path, "\n",
    " - Fixed/random spatial effects, hyperparameters, evaluation, diagnostic: ", values_path, "\n",
    " - Full model object: ", file.path(values_path, paste0(species, "_model_fit.rds")), "\n"
    ))
  }

  # summary          #@@@JMB pendiente revisar/completar...
  summary_df <- generate_summary_nsbm(fit, species, lcpo_val) 
  if(!is.null(cv_res)) {
    cv_rows <- data.frame(
      Field = c("","--------- Cross Validation --------", "CV folds:", "AUC mean ± sd:"),
      Value = c("","", cv_res$cv.folds, sprintf("%.2f ± %.2f", cv_res$auc_mean, cv_res$auc_sd)),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, cv_rows)
  }
  if(!is.null(diag_block)) {
    summary_df <- rbind(summary_df, diag_block$metrics)
  }

  # return
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
      pred_sp = if(!is.null(pred_sp)) terra::wrap(pred_sp) else NULL,
      pred_latent = if(!is.null(pred_lat)) terra::wrap(pred_lat) else NULL
    ),
    new.projections = if(length(proj_list) > 0) lapply(proj_list, terra::wrap) else list(),
    Summary = summary_df
  )
 
  attr(sabina, "class") <- "nsbm.inlabru"
  return(sabina)

}


### Helps/Auxiliars
# -----------------------------

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
        stop(paste0("❌ `", path, "` uses model='rw2' but is missing `u` and/or `alpha`",
             "  Use `list(model='rw2', u=..., alpha=...)`."))
      }
      stop(paste0("❌ `", path, "` has invalid value '", spec, "'. ",
                  "Use 'const'|'drop' or a list with `model='rw2'`, `u`, and `alpha`."))
    } else if(is.list(spec)) {
      m <- spec$model
      if(is.null(m)) stop(path, " is missing `model`.")
      if(identical(m, "const")) { res <- list(model = "const", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(m, "drop")) { res <- list(model = "drop", u = NA, alpha = NA); .spec_cache[[key]] <- res; return(res) }
      if(identical(m, "rw2")) {
        if(is.null(spec$u) || is.null(spec$alpha)) {
          stop(paste0("❌ `", path, "` with `model = 'rw2'` requires both `u` and `alpha`."))
        }
        res <- list(model = "rw2", u = spec$u, alpha = spec$alpha); .spec_cache[[key]] <- res; return(res)
      }
      stop(paste0("❌ `", path, "` has invalid `model = '", m, "'`. Use 'const'|'rw2'|'drop'."))
    } else {
      stop(paste0("❌ `", path, "` has an unsupported specification type."))
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
      stop("❌  rw2: covariate has zero range; rw2 requires variability.")
    }    
    g <- INLA::inla.group(vals, n = K, method = method)
    v <- sort(unique(as.numeric(levels(g))))
    v <- thin_knots(v, min_ratio = 1e-3)
    if(length(v) < 3L) {
      stop("❌  RW2: <3 knots after grouping/thinning.") #@@@JMB aquí también se podría ajustar K, pero demasiados args en mi opinión
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
    latent_block <- ""
    if(isTRUE(use_latent) && !is.null(spde_cov)) {
      latent_block <- "GLspde(main = geometry, model = spde_cov) + beta_GL(1)"
    }
    return(list(
      cmp = paste(Filter(nzchar, c(cmp1, latent_block, cmpglobal, cmpregional)), collapse = " + "),
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


# -----------------------------

# disgnstics
nsbm_diagnostics <- function(fit,
                             data_used,
                             sp_covreg,
                             priors = NULL,
                             pred_sp = NULL,
                             pred_lat = NULL) {

  spatial_local <- "spatial" %in% names(fit$summary.random)
  latent_global <- "GLspde" %in% names(fit$summary.random)
  hyp <- fit$summary.hyperpar

  val <- function(param) if (param %in% rownames(hyp)) hyp[param, "0.5quant"] else NA_real_
  ci_ratio <- function(param) {
    if (!param %in% rownames(hyp)) return(NA_real_)
    ciw <- hyp[param, "0.975quant"] - hyp[param, "0.025quant"]
    med <- hyp[param, "0.5quant"]
    if (!is.finite(ciw) || !is.finite(med) || med == 0) return(NA_real_)
    ciw / abs(med)
  }

  # extract hyperparameters
  range_res <- if(spatial_local) val("Range for spatial") else NA_real_
  sigma_res <- if(spatial_local) val("Stdev for spatial") else NA_real_
  range_lat <- if(latent_global) val("Range for GLspde") else NA_real_
  sigma_lat <- if(latent_global) val("Stdev for GLspde") else NA_real_

  # CI/median ratios
  ci_ratios <- c(range_residual = ci_ratio("Range for spatial"),
                 range_latent = ci_ratio("Range for GLspde"),
                 sigma_residual = ci_ratio("Stdev for spatial"),
                 sigma_latent = ci_ratio("Stdev for GLspde"))

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
    range_residual = if(!is.null(priors$spde.pcprior.range))
      close_rel(range_res, priors$spde.pcprior.range) else FALSE,
    sigma_residual = if(!is.null(priors$spde.pcprior.sigma))
      close_rel(sigma_res, priors$spde.pcprior.sigma) else FALSE,
    range_latent = if(!is.null(priors$latent.pcprior.range))
      close_rel(range_lat, priors$latent.pcprior.range) else FALSE,
    sigma_latent = if(!is.null(priors$latent.pcprior.sigma))
      close_rel(sigma_lat, priors$latent.pcprior.sigma) else FALSE
  )

  # variance and range ratios
  # sigma latent / sigma residual > 1.5 ==> latent field dominates (Blangiardo & Cameletti 2015)
  sigma_ratio <- if(spatial_local && latent_global && is.finite(sigma_lat) && is.finite(sigma_res) && sigma_res > 0)
    sigma_lat / sigma_res else NA_real_
  # 0.5 < range_latent / range_residual < 2 ==> poor scale separation (Bakka et al. 2018)
  range_ratio <- if(spatial_local && latent_global && is.finite(range_lat) && is.finite(range_res) && range_res > 0)
    range_lat / range_res else NA_real_

  # correlation spatial–latent fields
  # r > 0.7 00> high correlation
  field_correlation <- NA_real_
  if(spatial_local && latent_global) {
    f_sp <- fit$summary.random$spatial$mean
    f_lat <- fit$summary.random$GLspde$mean
    n <- min(length(f_sp), length(f_lat))
    field_correlation <- stats::cor(f_sp[seq_len(n)], f_lat[seq_len(n)], use = "pairwise.complete.obs")
  } 
    
  # Moran’s I residual autocorrelation
  moran_I <- NA_real_
  if(all(c("x", "y", "resp") %in% names(data_used))) {
    pred_mu <- fit$summary.fitted.values$mean
    if(length(pred_mu) >= nrow(data_used)) {
      rs <- data_used$resp - pred_mu[seq_len(nrow(data_used))]
      xy <- as.matrix(data_used[, c("x", "y")])
      #maxdist <- if(spatial_local) range_res else if(latent_global) range_lat else NA_real_
      maxdist <- if (spatial_local) {
        range_res
      } else if (latent_global) {
        range_lat
      } else {               #@@@JMB sin spatial ni latent usa 1/4 de la diagonal???
        bb <- apply(xy, 2, range, na.rm = TRUE)
        sqrt(sum((bb[2,] - bb[1,])^2)) / 4
      }   
      if(is.finite(maxdist) && maxdist > 0) {
        nb <- spdep::dnearneigh(xy, 0, maxdist, longlat = FALSE)
        lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
        mi <- spdep::moran(rs, lw, n = length(rs), S0 = spdep::Szero(lw))
        moran_I <- as.numeric(mi$I)
      }
    }
  }

  # bayesian pseudo-R2, correlation, RMSE
  y_obs <- data_used$resp
  y_pred <- fit$summary.fitted.values$mean[seq_len(nrow(data_used))]
  pseudo_R2 <- 1 - (stats::var(y_obs - y_pred, na.rm = TRUE) / stats::var(y_obs, na.rm = TRUE))
  pred_cor <- stats::cor(y_obs, y_pred, use = "complete.obs")
  rmse <- sqrt(mean((y_obs - y_pred)^2, na.rm = TRUE))

  # warnings
  warns <- character()
  if(is.finite(max_CIratio) && max_CIratio > 25)   #@@@JMB rev CI/median thresholds 10–25–50 Bakka et al. 2018) 
    warns <- c(warns, "⚠️ Weak identifiability: max CI/median > 25. Strengthen priors or simplify model.")
  if(spatial_local && latent_global && is.finite(sigma_ratio) && sigma_ratio > 1.5)
    warns <- c(warns, "⚠️ Latent field dominates variance (σ_latent > 1.5 × σ_residual). Reduce σ_latent or increase its range.")
  if(spatial_local && latent_global && is.finite(range_ratio) && range_ratio > 0.5 && range_ratio < 2)
    warns <- c(warns, "⚠️ Poor range separation (0.5 < latent/residual < 2). Enforce latent range ≥ 3–5× residual range.")
  if (!is.null(field_correlation) && is.finite(field_correlation) && field_correlation > 0.7)
    warns <- c(warns, "⚠️ High correlation between latent and residual fields (r > 0.7). Possible redundancy.")
  if(is.finite(moran_I) && moran_I > 0.10)
    warns <- c(warns, paste0("⚠️ Residual spatial autocorrelation (Moran’s I ≈ ", round(moran_I, 2), "). Model may miss local dependence."))
  if(any(unlist(prior_close)))
    warns <- c(warns, "⚠️ Posterior ≈ prior thresholds detected: weak data information for some hyperparameters.")

  #
  diag_metrics <- data.frame(
    Field = c(
      "",
      "-------- Model diagnostics --------",
      "  Parameter identifiability:",
      "      CI/median (Range: residual | latent)",
      "      CI/median (σ: residual | latent)",
      "      Max CI/median ratio (overall)",
      "",
      "  Spatial structure:",
      "      Field correlation (latent-residual, r)",
      "      Range ratio (latent / residual)",
      "      Variance ratio (σ_latent / σ_residual)",
      "",
      "  Residual spatial dependence:",
      "      Residual Moran’s I",
      "",
      "  Predictive performance:",
      "      Bayesian pseudo-R²",
      "      Observed-predicted correlation (r)",
      "      RMSE"
    ),
    Value = c(
      "", "", "",
      paste0(round(ci_ratios["range_residual"], 1), " | ",
             round(ci_ratios["range_latent"], 1)),
      paste0(round(ci_ratios["sigma_residual"], 1), " | ",
             round(ci_ratios["sigma_latent"], 1)),
      round(max_CIratio,2),
      "",
      "",
      round(field_correlation,2),
      round(range_ratio,2),
      round(sigma_ratio,2),
      "",
      "",
      round(moran_I,3),
      "",
      "",
      round(pseudo_R2,3),
      round(pred_cor,3),
      round(rmse,3)
    ),
    stringsAsFactors = FALSE
  )


  # plots
  # pA Hyperparameters: posterior marginals + PC prior thresholds
  post_df <- NULL
  if(!is.null(fit$marginals.hyperpar)) {
    for(nm in names(fit$marginals.hyperpar)) {
      sm <- INLA::inla.smarginal(fit$marginals.hyperpar[[nm]])
      post_df <- rbind(post_df, data.frame(x = sm$x, y = sm$y, par = nm))
    }
  }

  # Create prior reference ticks for red dashed lines
  prior_ticks <- do.call(rbind, Filter(Negate(is.null), list(
    if(!is.null(priors$spde.pcprior.range)) data.frame(x = priors$spde.pcprior.range[1], par = "Range for spatial"),
    if(!is.null(priors$spde.pcprior.sigma)) data.frame(x = priors$spde.pcprior.sigma[1], par = "Stdev for spatial"),
    if(!is.null(priors$latent.pcprior.range)) data.frame(x = priors$latent.pcprior.range[1], par = "Range for GLspde"),
    if(!is.null(priors$latent.pcprior.sigma)) data.frame(x = priors$latent.pcprior.sigma[1], par = "Stdev for GLspde")
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
                       ggplot2::aes(x = x, y = y, label = "PC prior (u)"),
                       vjust = -0.8, hjust = 0.9, size = 2.8,
                       color = "#c0392b", angle = 90) +
    ggplot2::labs(
      title = "A) Hyperparameters: posterior marginals",
      subtitle = "Red dashed lines = PC prior thresholds (u)",
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

  # pB spatial residual vs latent fields
  ras_to_df <- function(r, nm) {
    rr <- terra::unwrap(r)[["mean"]]
    df <- terra::as.data.frame(rr, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", "mean")
    df$which <- nm
    df
  }
  maps_df <- data.frame()
  if(!is.null(pred_sp)) maps_df <- rbind(maps_df, ras_to_df(pred_sp, "Residual spatial field"))
  if(!is.null(pred_lat)) maps_df <- rbind(maps_df, ras_to_df(pred_lat, "Latent global field"))

  hierarchical <- any(grepl("copy", rownames(fit$summary.hyperpar), ignore.case = TRUE)) ||
                   "IGlobal" %in% names(fit$summary.random)
  # Define midpoint dynamically
  if(hierarchical) {
    midpoint_val <- mean(maps_df$mean, na.rm = TRUE)  # recentra en la media estimada
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
        title = "B) Spatial fields (posterior mean)",
        subtitle = if(hierarchical)
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
        title = "B) Spatial fields (posterior mean)",
        subtitle = "No spatial or latent fields present in the model"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, color = "#2c3e50"),
        plot.subtitle = ggplot2::element_text(size = 9, color = "#5d6d7e")
      )
  }

  # pC Residual correlogram
  cor_df <- data.frame()
  if(is.finite(moran_I)) {
    coords <- as.matrix(data_used[, c("x", "y")])
    rs <- data_used$resp - fit$summary.fitted.values$mean[seq_len(nrow(data_used))]
    dmat <- as.matrix(stats::dist(coords))
    dvec <- dmat[upper.tri(dmat)]
    nbins <- 15L
    brks <- seq(0, max(dvec, na.rm = TRUE), length.out = nbins + 1)
    mid <- 0.5 * (brks[-1] + brks[-length(brks)])
    rho <- rep(NA_real_, nbins)
    ut_r <- row(dmat)[upper.tri(dmat)]
    ut_c <- col(dmat)[upper.tri(dmat)]
    for (i in seq_len(nbins)) {
      sel <- dvec >= brks[i] & dvec < brks[i + 1]
      if (sum(sel) > 30) {
        idx <- which(sel)
        r1 <- rs[ut_r[idx]]; r2 <- rs[ut_c[idx]]
        rho[i] <- stats::cor(r1, r2, use = "pairwise.complete.obs")
      }
    }
    cor_df <- data.frame(dist_mid = mid, rho = rho)
  }
  range_eff <- if(is.finite(range_res)) range_res else if(is.finite(range_lat)) range_lat else NA_real_
  pC <- ggplot2::ggplot(cor_df, ggplot2::aes(x = dist_mid, y = rho)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(na.rm = TRUE, size = 1.2, color = "#1a5276") +
    ggplot2::geom_line(na.rm = TRUE, color = "#1a5276", linewidth = 0.6) +
    ggplot2::labs(
      title = "C) Residual correlogram",
      subtitle = if(is.finite(range_eff))
        paste0("Model range ≈ ", round(range_eff, 3), " (map units)\n(distance where correlation vanishes)")
      else
        "No spatial/latent field: full extent shown",
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
  rs <- data_used$resp - fit$summary.fitted.values$mean[seq_len(nrow(data_used))]
  res_mean <- mean(rs, na.rm = TRUE)
  df_r <- data.frame(resid = rs)
  center_label <- if (hierarchical) {
    paste0("Mean residual (≈ ", round(res_mean, 3), ")")
  } else {
    "Model mean = 0"
  }
  line_x <- if(hierarchical) res_mean else 0

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
      title = "D) Residual histogram + Q-Q plot",
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

  composite <- patchwork::wrap_plots(pA, pB, (pC + pD1 + pD2), ncol = 1)+
    patchwork::plot_layout(heights = c(1.4, 0.8, 1)
  )

  # output
  out <- list(
    metrics = diag_metrics,
    ci_ratios = ci_ratios,
    warnings = warns,
    composite = composite
  )

  return(out)
}



###----------------###




