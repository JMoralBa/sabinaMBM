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
#' @param coupling.intercept Logical; if TRUE = model regional intercept as deviation from global.
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
                      coupling.intercept = "additive",    # NULL regional-only, "additive" interceptos independientes, ""hierarchical" IRegional se modela como desviación de IGlobal
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
    stop("❌  The 'nsbm_obj' must be of class 'nsdm.vinput'. Please see sabinaNSDM::NSDM.SelectCovariates().")
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
         "  Supported families are: ", paste(names(valid_links), collapse = ", "), "\n")
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    stop("❌ Inalid link `", lnk, "` for this family `", fam, "`.\n",
         "  Allowed links for '", fam, "': ", paste(valid_links[[fam]], collapse = ", "), ".\n")
  }
  if(is.null(coupling.intercept)) {
    if(length(nsbm_obj$Selected.Variables.Regional) == 0) {
      stop("❌ 'coupling.intercept = NULL' (regional-only) but no regional covariates are present in 'nsbm_obj'.\n\n")
    }
    if(latent_global) {
      stop("❌ Latent global SPDE not allowed in regional-only model (no global scale exists).\n\n")
    }
    message("ℹ️ 'coupling.intercept = NULL'; running a regional-only model (no global component will be used).\n\n")
  } else {
    if(!is.character(coupling.intercept) || !coupling.intercept %in% c("additive", "hierarchical")) {
      stop("❌  The 'coupling.intercept' must be NULL, 'additive', or 'hierarchical'.\n\n")
    }
    if(is.null(spde.mesh)) {
      stop("❌  You provided SPDE priors but no 'spde.mesh'.\n",
           "  To include them, create a mesh with `create_mesh()`.\n\n")
    }
    if(length(nsbm_obj$Selected.Variables.Global) == 0) {
      stop("❌ 'coupling.intercept = ", coupling.intercept, "' requires a global component,\n",
           "   but no global covariates are present in 'nsbm_obj'.\n",
           "   Use 'coupling.intercept = NULL' for a regional-only model.\n\n")
    }
  }
  if (is.null(spde.mesh) && (spatial_local || latent_global)) {
    stop("❌ SPDE priors were provided but 'spde.mesh' is NULL. Create a mesh with `create_mesh()` and pass it to `spde.mesh`.\n\n")
  }
  if(!is.null(spde.mesh) && !inherits(spde.mesh, "inla.mesh")) {
    stop("❌ `spde.mesh` must be a valid INLA mesh object (class 'inla.mesh'). Use `create_mesh()` to build it.\n\n")
  }
  if(spatial_local && (is.null(spde.pcprior.range) || is.null(spde.pcprior.sigma))) {
    stop("❌ Missing priors for spatial field. Must provide both `spde.pcprior.range` and `spde.pcprior.sigma`.\n\n")
  }
  if(latent_global && (is.null(latent.pcprior.range) || is.null(latent.pcprior.sigma))) {
      stop("❌ Missing priors for latent field. Mus provide both `latent.pcprior.range` and `latent.pcprior.sigma`.\n\n")
  }
  if(spatial_local && latent_global) {
    if(latent.pcprior.range[1] < spde.pcprior.range[1] * 3) {
      warning("⚠️ `latent.pcprior.range` < 3× `spde.pcprior.range`: fields may overlap, causing double-counting of spatial variance.\n\n")
    }
  }
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1) {
      stop("❌ 'seed' must be a single numeric value.\n\n")
    }
    set.seed(seed)
  }
  available_cores <- parallel::detectCores(logical = TRUE)
  if(!is.null(n.threads) && n.threads > available_cores) {
    stop(paste0("❌ Requested `n.threads` = ", n.threads, " exceeds available cores (", available_cores,").\n\n"))
  }
  INLA::inla.setOption(num.threads = n.threads)
#
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
  if(latent_global && !is.null(covariate.effects)) {
    vg_all <- nsbm_obj$Selected.Variables.Global
    if(length(vg_all) > 0) {
      has_rw2_global <- any(vapply(
        vg_all,
        function(v) {
          spec <- resolve_spec(v, "global", covariate.effects)
          identical(spec$model, "rw2")
        },
        logical(1)
      ))
      if(has_rw2_global) {
        stop("❌ RW2 global effects are not allowed when the latent SPDE is active.\n",
             "   Both mechanisms smooth broad-scale structure and become redundant.\n",
             "   Disable RW2 for global covariates (in `covariate.effects`) or disable the latent SPDE\n",
             "   (`latent.pcprior.range` / `latent.pcprior.sigma` = NULL).\n\n")
      }
    }
  }

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


  # SPDE components
  if(spatial_local)  {
    matern <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = spde.pcprior.range,
      prior.sigma = spde.pcprior.sigma
    )
  }

  if(latent_global) {
    spde_cov <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = latent.pcprior.range,
      prior.sigma = latent.pcprior.sigma)
  } else {
    spde_cov <- NULL
  }


  # Nested intercept
  IID_PC_PRIOR <- "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"
  if(is.null(coupling.intercept)) {
    #intercept_terms <- c("IRegional(1, model='iid')")
    intercept_terms <- c(paste0("IRegional(1, model='iid', ", IID_PC_PRIOR, ")"))
  } else if(coupling.intercept == "additive") {
    intercept_terms <- c(
      #"IGlobal(1, model='iid')",
      #"IRegional(1, model='iid')"
      paste0("IGlobal(1, model='iid', ", IID_PC_PRIOR, ")"),
      paste0("IRegional(1, model='iid', ", IID_PC_PRIOR, ")")
    )
  } else if(coupling.intercept == 'hierarchical') {
    intercept_terms <- c(
      #"IGlobal(1, model='iid')",
      paste0("IGlobal(1, model='iid', ", IID_PC_PRIOR, ")"),
      paste0("IRegional(1, copy='IGlobal', fixed=FALSE, ",
             "hyper=list(beta=list(prior='normal', param=c(1,0.001))))"))
  }
  base_intercepts <- paste(intercept_terms, collapse = " + ")


  # Model components
  cmp_cov <- fcov(obj = nsbm_obj, 
                  spobjglo = "sp_covglo", 
                  spobjreg = "sp_covreg",
                  sp_covglo = sp_covglo,
                  sp_covreg = sp_covreg,
                  covariate.effects = covariate.effects,
                  pp_glo_sf = pp_glo,
                  pp_reg_sf = pp_reg)


  # component formula (intercept + spataial + latent + fcovs
  cmp <- c(
    base_intercepts,
    if(spatial_local) "spatial(geometry, model = matern)" else NULL, # residual field local (mesh)
    if(latent_global) "GLspde(main = geometry, model = spde_cov) + beta_GL(1)" else NULL, # global field + coef
    cmp_cov$cmp   # bio1GL(), bio1RE(),...
  )
  cmp <- paste(cmp[!is.na(cmp) & nzchar(cmp)], collapse = " + ")
  cmp <- as.formula(paste("~", cmp))

  dom <- if (spatial_local || latent_global) list(geometry = spde.mesh) else NULL
  f_spatial <- if(spatial_local) " + spatial" else ""
  f_latent <- if(latent_global) " + beta_GL * GLspde" else ""

  .opt_plus <- function(s) if(!is.null(s) && nzchar(s)) paste0(" + ", s) else ""
 
  rhs_glo <- if(is.null(coupling.intercept)) {
    NULL  # no for regional-only
  } else {
    paste0(
      "IGlobal",
      f_spatial,
      f_latent,
      .opt_plus(cmp_cov$like$fglobal)
    )
  }
  rhs_reg <- if(is.null(coupling.intercept)) {
    paste0(
      "IRegional",
      f_spatial,
      f_latent,
      .opt_plus(cmp_cov$like$fregional)
    )
  } else {
    paste0(
      "IGlobal + IRegional",
      f_spatial,
      f_latent,
      .opt_plus(cmp_cov$like$fregional)
    )
  }

  # likelihoods
  if(fam == "cp") {
    # for intensity, Cox process (cp)
    pres_glo <- pp_glo[pp_glo$resp != 0L, ]
    pres_reg <- pp_reg[pp_reg$resp != 0L, ]
    if(!is.null(coupling.intercept)) {
      lik_glo <- inlabru::like(
        family = "cp",
        formula = as.formula(paste0("geometry ~ ", rhs_glo)),
        data = pres_glo,
        samplers = bdy_glo,
        domain = dom)
    } else {
      lik_glo <- NULL
    }
    lik_reg <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ ", rhs_reg)),
      data = pres_reg,
      samplers = bdy_reg,
      domain = dom)
  } else {
    # for any type of family-object (e.g, binomial(), etc.)
    if(!is.null(coupling.intercept)) {
      lik_glo <- inlabru::like(
        family = fam,
        formula = as.formula(paste0("resp ~ ", rhs_glo)),
        data = pp_glo,
        samplers = bdy_glo,
        domain = dom,
        control.family = list(link = lnk))
    } else {
      lik_glo <- NULL
    }
    lik_reg <- inlabru::like(
      family = fam,
      formula = as.formula(paste0("resp ~ ", rhs_reg)),
      data = pp_reg,
      samplers = bdy_reg,
      domain = dom,
      control.family = list(link = lnk))
  }

  eta_terms <- c(
    if(!is.null(coupling.intercept)) "IGlobal" else NULL,
    "IRegional",
    if(spatial_local) "spatial" else NULL,
    if(latent_global) "beta_GL * GLspde" else NULL,
    cmp_cov$like$fregional
  )
  eta <- paste(stats::na.omit(eta_terms), collapse = " + ")

  if(fam == "cp") {
    pred_formula <- as.formula(paste0("~ exp(",eta,")"))
  } else {
    pred_formula <- switch(
      lnk,
      "logit" = as.formula(paste0("~ 1 / (1 + exp(-(", eta, ")))")),
      "cloglog" = as.formula(paste0("~ 1 - exp(-exp(", eta, "))")), 
      "log" = as.formula(paste0("~ exp(", eta, ")")),
      "identity" = as.formula(paste0("~ ", eta)))
  } 


  # Model fitting complete
  if(is.null(coupling.intercept)) {
    fit <- inlabru::bru(
      components = cmp,
      lik_reg,
      options = list(
        control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
        control.inla = list(int.strategy = "eb"),
        control.mode = list(restart = TRUE)))
  } else {
    fit <- inlabru::bru(
      components = cmp,
      lik_glo,
      lik_reg,
      options = list(
        control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
        control.inla = list(int.strategy = "eb"),
        control.mode = list(restart = TRUE)))
  }


  # k-fold CV (only fam no cp)
  if(cv.folds > 1) {
    if (fam == "cp") {
      warning("⚠️ Cross-validation (cv.folds > 1) is not implemented for family = 'cp'. CV results will be NULL.")
      cv_res <- NULL
    } else {
      # stratified k folds
      if(!is.null(coupling.intercept)) {
        folds_g <- make_stratified_kfolds(pp_glo$resp, cv.folds)
      } else {
        folds_g <- NULL
      }
      folds_r <- make_stratified_kfolds(pp_reg$resp, cv.folds)

      aucs <- numeric(cv.folds)

      for(k in seq_len(cv.folds)) {
        if(!is.null(coupling.intercept)) {
          train_g <- pp_glo[folds_g != k, ]
          test_g <- pp_glo[folds_g == k, ]
        } else {
          train_g <- NULL
          test_g  <- NULL
        }    
        train_r <- pp_reg[folds_r != k, ]
        test_r <- pp_reg[folds_r == k, ]

        # likelihood k
        if(!is.null(coupling.intercept)) {
          lik_g_k <- inlabru::like(
            family = fam,
            formula = lik_glo$formula,
            data = train_g,
            samplers = lik_glo$samplers,
            domain = lik_glo$domain,
            control.family = list(link = lnk)
          )
        } else {
          lik_g_k <- NULL
        }
        lik_r_k <- inlabru::like(
          family  = fam,
          formula = lik_reg$formula,
          data = train_r,
          samplers = lik_reg$samplers,
          domain = lik_reg$domain,
          control.family = list(link = lnk)
        )

        # fit fold k
        if(is.null(coupling.intercept)) {
          fit_k <- inlabru::bru(
            components = cmp,
            lik_r_k,
            options = list(control.compute = list(cpo = FALSE, config = TRUE),
                           control.inla = list(int.strategy = "eb"),
                           control.mode = list(restart = TRUE)))
        } else {
          fit_k <- inlabru::bru(
            components = cmp,
            lik_g_k,
            lik_r_k,
            options = list(control.compute = list(cpo = FALSE, config = TRUE),
                           control.inla = list(int.strategy = "eb"),
                           control.mode = list(restart = TRUE)))
        }

        # rm NAs
        coords_t <- sf::st_coordinates(test_r)
        both_cov <- c(sp_covglo, sp_covreg)
        cov_all <- terra::extract(both_cov, coords_t)
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
    }
  } else {
    cv_res <- NULL
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
                  pred_lat = if(latent_global) pred_lat else NULL,
                  coupling.intercept = coupling.intercept)

  if(length(diag_block$warnings)) {
    message(paste(unique(diag_block$warnings), collapse = "\n\n"))
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
      Value = c(fit$waic$waic, fit$dic$dic, fit$mlik[1], diag_block$bayes_fit$lcpo_val)
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
    # save diagnostic plots
   if(!is.null(diag_block$plots$hyperparams)) {
     ggplot2::ggsave(
       file.path(values_path, paste0(species, "_hyperparams.png")),
       plot = diag_block$plots$hyperparams,
       width = 6, height = 5, dpi = 300
      )
    }
    if(!is.null(diag_block$plots$correlogram)) {
      ggplot2::ggsave(
        file.path(values_path, paste0(species, "_correlogram.png")),
        plot = diag_block$plots$correlogram,
        width = 6, height = 5, dpi = 300
      )
    }
    if(!is.null(diag_block$plots$hist)) {
      ggplot2::ggsave(
        file.path(values_path, paste0(species, "_residualHistogram.png")),
        plot = diag_block$plots$hist,
        width = 6, height = 5, dpi = 300
      )
    }
    if(!is.null(diag_block$plots$qq)) {
      ggplot2::ggsave(
        filename = file.path(values_path, paste0(species, "_QQplot.png")),
        plot = diag_block$plots$qq,
        width = 6, height = 4, dpi = 300, bg = "white"
      )
    }
    if(!is.null(diag_block$plots$spatialfields)) { 
      ggplot2::ggsave(
        file.path(values_path, paste0(species, "_SPDEfields.png")),
        plot = diag_block$plots$spatialfields,
        width = 6, height = 5, dpi = 300
      )
    }
    if(!is.null(diag_block$plots$semivariogram)) {
      ggplot2::ggsave(
        filename = file.path(values_path, paste0(species, "_semivariogram.png")),
        plot = diag_block$plots$semivariogram,
        width = 6, height = 4, dpi = 300, bg = "white"
      )
    }

    message("ℹ️ Results saved in the following local folder(s):")
    message(paste(
    "  - Projections (current pred, spatial field pred_sp, latent field pred_lat and scenarios: ", projections_path, "\n",
    " - Model values (fixed/random effects, hyperparameters, evaluation and diagnostics: ", values_path, "\n",
    " - Full model object (.rds): ", file.path(values_path, paste0(species, "_model_fit.rds")), "\n"
    ))
  }


  # summary          #@@@JMB pendiente revisar/completar...
  summary_df <- nsbm_generate_summary(fit, species, fam, lnk, 
                                      diag_block, cv_res = cv_res)


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
      coupling.intercept = coupling.intercept,
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

# Components GL/RE for covariates
fcov <- function(obj, 
                 spobjglo, 
                 spobjreg,
                 sp_covglo,
                 sp_covreg,
                 covariate.effects = NULL,
                 pp_glo_sf = NULL,
                 pp_reg_sf= NULL) {

  # selected vars
  vg <- obj$Selected.Variables.Global
  vr <- obj$Selected.Variables.Regional

  cmp1 <- paste0(unique(c(vg, vr)), "(1)", collapse = " + ")

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

  build_rw2_values <- function(rast_layer, coords, K = rw2_K, method = rw2_method) {
    vals <- suppressWarnings(as.numeric(terra::extract(rast_layer, coords)[, 1]))
    vals <- vals[is.finite(vals)] 
    rng <- range(vals)
    if(diff(rng) == 0) {
      stop("❌  RW2 requires variability (covariate range = 0).\n",
           "   Use 'const' instead for this covariate or check the raster values.\n\n")
    }  
    g <- INLA::inla.group(vals, n = K, method = method)
    v <- sort(unique(as.numeric(levels(g))))
    v <- thin_knots(v)
    if(length(v) < 3L) {
      stop("❌  RW2 requires at least 3 distinct support values for the covariate.\n", #@@@JMB aquí también se podría ajustar K, pero demasiados args en mi opinión
           "   Consider lowering smoothing or checking covariate variability.\n\n")
    }
    v
  }

  # detect if rw2 is request
  need_rw2_global <- FALSE
  need_rw2_regional <- FALSE
  if(is.list(covariate.effects)) {
    if(length(vg) > 0) {
      need_rw2_global <- any(vapply(vg, function(x) 
        resolve_spec(x, "global",  covariate.effects)$model == "rw2", logical(1)))
    }
    if(length(vr) > 0) {
      need_rw2_regional <- any(vapply(vr, function(x) 
        resolve_spec(x, "regional", covariate.effects)$model == "rw2", logical(1)))
    }
  }

  coords_all <- if(need_rw2_global) {
    rbind(sf::st_coordinates(pp_glo_sf),
          sf::st_coordinates(pp_reg_sf))
  } else NULL
  coords_r <- if(need_rw2_regional) {
    sf::st_coordinates(pp_reg_sf)
  } else NULL

  # components global
  cmpglobal <- character(0)
  fglobal <- character(0)

  fmt_vals <- function(v) {
    paste0("c(", paste(format(v, digits = 7), collapse = ","), ")")
  }

  for(X in vg) {
    specX <- resolve_spec(X, "global", covariate.effects)
    if(specX$model == "drop") next
    if(specX$model == "const") {
      cmpglobal <- c(cmpglobal,
        paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'const')"))
    } else if(specX$model == "rw2") {
      vals <- build_rw2_values(sp_covglo[[X]], coords_all)
      cmpglobal <- c(cmpglobal,
        paste0(X, "GL(main = ", spobjglo, ", main_layer = '", X, "', model = 'rw2', scale.model = TRUE, ",
               "hyper = list(prec = list(prior='pc.prec', param=c(", specX$u, ",", specX$alpha, "))), ",
               "values = ", fmt_vals(vals),")"))
    }
    fglobal <- c(fglobal, paste0(X, " * ", X, "GL"))
  }

  # components regional
  cmpregional <- character(0)
  fregional <- character(0)

  for(X in vr) {
    specX <- resolve_spec(X, "regional", covariate.effects)
    if(specX$model == "drop") next
    if(specX$model == "const") {
      cmpregional <- c(cmpregional,
        paste0(X, "RE(main = ", spobjreg, ", main_layer = '", X, "', model = 'const')"))
    } else if(specX$model == "rw2") {
      vals <- build_rw2_values(sp_covreg[[X]], coords_r)
      cmpregional <- c(cmpregional,
        paste0(X, "RE(main = ", spobjreg,
          ", main_layer = '", X, "', model='rw2', scale.model=TRUE, ",
          "hyper = list(prec=list(prior='pc.prec', param=c(",
          specX$u, ",", specX$alpha, "))), ",
          "group = region, ",
          "values = ", fmt_vals(vals), ")"))
    }
    fregional <- c(fregional, paste0(X, " * ", X, "RE"))
  }

  # out fcov
  list(
    #cmp = paste(c(cmpglobal, cmpregional), collapse = " + "),
    cmp = paste(c(cmp1, cmpglobal, cmpregional), collapse = " + "),
    like = list(fglobal = if(length(fglobal) > 0) paste(fglobal, collapse = " + ") else "",
                fregional = if(length(fregional) > 0) paste(fregional, collapse = " + ") else "")) 
}


# -----------------------------


# interprete covariate_effects
resolve_spec <- function(varname, scale, covariate.effects, default_model = "const") {

  path_label <- paste0("covariate.effects$", scale, "$", varname)

  # If no covariate.effects, all const by default
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
    return(list(model = "const", u = NA, alpha = NA))
  }

  # if spec is a string -> "const"/"drop"
  if(is.character(spec)) {
    if(spec == "const") { 
      res <- list(model = "const", u = NA, alpha = NA)
    }
    if(spec == "drop") { 
      res <- list(model = "drop", u = NA, alpha = NA)
    } 
    if(spec == "rw2") {
      stop("❌  Invalid RW2 specification in `", path_label, "`.\n",
           "   RW2 must be expressed as a list: list(model='rw2', u=..., alpha=...).\n",
           "   For linear effects use 'const'; to exclude use 'drop'.\n\n")
    } 
    stop("❌ Invalid keyword '", spec, "' in `", path_label, "`.\n",
         "   Valid options: 'const', 'drop', or list(model='rw2', u=..., alpha=...) (smoothed RW2 effect).\n\n")
  }

  # if spec is a list, only rw2 admited 
  if(is.list(spec)) {
    if(is.null(spec$model)) {
      stop("❌  Invalid list specification in `", path_label, "`.\n",
           "   Lists are only allowed for RW2 effects. Use: list(model='rw2', u=..., alpha=...).\n",
           "   For linear effects use 'const'; to exclude use 'drop'.\n\n")
    }
    if(!identical(spec$model, "rw2")) {
      # RW2 reequieres u alpha
      if(is.null(spec$u) || is.null(spec$alpha)) {
        stop("❌  Incomplete RW2 specification in `", path_label, "`.\n",
             "   Provide both `u` and `alpha`..\n\n")
      }
      return(list(model = "rw2", u = spec$u, alpha = spec$alpha))
    }
    stop("❌  Invalid `model` in `", path_label, "`.\n",
         "   When using a list, only model='rw2' is permitted.\n",
         "   For linear effects use 'const'; to exclude use 'drop'.\n\n")
  }

  #unsupported type
  stop("❌  Unsupported type in `", path_label, "`.\n",
       "   Must be string ('const'/'drop') or list(model='rw2', u=..., alpha=...).\n\n")
}


# -----------------------------


# diagnstics
nsbm_diagnostics <- function(fit,
                             data_used,
                             priors = NULL,
                             pred_sp = NULL,
                             pred_lat = NULL,
                             coupling.intercept) {

  spatial_local <- "spatial" %in% names(fit$summary.random)
  latent_global <- "GLspde" %in% names(fit$summary.random)

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
  range_res_mean <- if(spatial_local) hyper["Range for spatial", "mean"] else "—"
  sigma_res_mean <- if(spatial_local) hyper["Stdev for spatial", "mean"] else "—"
  range_lat_mean <- if(latent_global) hyper["Range for GLspde", "mean"] else "—"
  sigma_lat_mean <- if(latent_global) hyper["Stdev for GLspde", "mean"] else "—"
  range_res_sd <- if(spatial_local) hyper["Range for spatial", "sd"] else "—"
  sigma_res_sd <- if(spatial_local) hyper["Stdev for spatial", "sd"] else "—"
  range_lat_sd <- if(latent_global) hyper["Range for GLspde", "sd"] else "—"
  sigma_lat_sd <- if(latent_global) hyper["Stdev for GLspde", "sd"] else "—"

  range_res <- if(spatial_local) paste0(round(range_res_mean, 2), " ± ", round(range_res_sd, 2)) else "—"
  sigma_res <- if(spatial_local) paste0(round(sigma_res_mean, 2), " ± ", round(sigma_res_sd, 2)) else "—"
  range_lat <- if(latent_global) paste0(round(range_lat_mean, 2), " ± ", round(range_lat_sd, 2)) else "—"
  sigma_lat <- if(latent_global) paste0(round(sigma_lat_mean, 2), " ± ", round(sigma_lat_sd, 2)) else "—"
  format_hyper <- function(param) {
   if(!param %in% rownames(hyper)) return("—")
      paste0(round(hyper[param, "mean"], 2), " ± ", round(hyper[param, "sd"], 2))
  }
  prec_IGlobal <- if("Precision for IGlobal" %in% rownames(hyper)) format_hyper("Precision for IGlobal") else "—"
  beta_IRegional <- if("Beta for IRegional" %in% rownames(hyper)) format_hyper("Beta for IRegional") else "—"

  # spatial and latent overlap?
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
  sig_vars <- signif_vars(fit)

  y_obs <- data_used$resp
  y_pred <- fit$summary.fitted.values$mean[seq_len(nrow(data_used))]

  # auc_full (if binary)
  if(all(y_obs %in% c(0,1))) {
    auc_full <- suppressMessages(as.numeric(pROC::auc(y_obs, y_pred)))
  } else {
    auc_full <- "—"
  }

  # brier
  brier <- mean((y_pred - y_obs)^2, na.rm = TRUE)

  # calibration & slopes
  if(all(y_pred > 0 & y_pred < 1)) {
    df_cal <- data.frame(
      logit_p = qlogis(y_pred),
      y = y_obs
    )
    cal_mod <- lm(logit_p ~ y, data = df_cal)
    cal_intercept <- coef(cal_mod)[1]
    cal_slope <- coef(cal_mod)[2]
  } else {
    cal_intercept <- "—"
    cal_slope <- "—"
  }

  # pit
  pit_vals <- fit$cpo$pit
  if(!is.null(pit_vals)) {
    ks_pit <- suppressWarnings(ks.test(pit_vals, "punif")$p.value)
  } else {
    ks_pit <- "—"
  }

  # cpo_failures
  if(!is.null(fit$cpo$cpo)) {
    cpo_failures <- sum(fit$cpo$failure, na.rm = TRUE)
  } else {
    cpo_failures <- "—"
  }

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
  # sigma latent mean / sigma residual mean > 1.5 ==> latent field dominates (Blangiardo & Cameletti 2015)
  sigma_ratio <- if(spatial_local && latent_global && is.finite(sigma_lat_mean) && is.finite(sigma_res_mean) && sigma_res_mean > 0)
    sigma_lat_mean / sigma_res_mean else "—"
  # 0.5 < range_latent / range_residual < 2 ==> poor scale separation (Bakka et al. 2018)
  range_ratio <- if(spatial_local && latent_global && is.finite(range_lat_mean) && is.finite(range_res_mean) && range_res_mean > 0)
    range_lat_mean / range_res_mean else "—"

  # correlation spatial–latent fields
  # r > 0.7 00> high correlation
  field_correlation <- "—"
  if(spatial_local && latent_global) {
    f_sp <- fit$summary.random$spatial$mean
    f_lat <- fit$summary.random$GLspde$mean
    n <- min(length(f_sp), length(f_lat))
    field_correlation <- stats::cor(f_sp[seq_len(n)], f_lat[seq_len(n)], use = "pairwise.complete.obs")
  } 
    
  # Moran’s I residual autocorrelation
  moran_I <- NA_real_
  rs <- y_obs - y_pred
  xy <- as.matrix(data_used[, c("x", "y")])
  maxdist <- if(spatial_local) {
    range_res_mean
  } else if(latent_global) {
    range_lat_mean
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

  # pred correlation, RMSE
  pred_cor <- stats::cor(y_obs, y_pred, use = "complete.obs")
  rmse <- sqrt(mean((y_obs - y_pred)^2, na.rm = TRUE))

  # coverage 95%
  if(all(c("0.025quant", "0.975quant") %in% colnames(fit$summary.fitted.values))) {
    low95 <- fit$summary.fitted.values[seq_len(nrow(data_used)), "0.025quant"]
    up95 <- fit$summary.fitted.values[seq_len(nrow(data_used)), "0.975quant"]
    cov95 <- mean(y_obs >= low95 & y_obs <= up95, na.rm = TRUE)
  } else {
    cov95 <- "—"
  }

  # coverage 50%
  if(all(c("0.25quant", "0.75quant") %in% colnames(fit$summary.fitted.values))) {
    low50 <- fit$summary.fitted.values[seq_len(nrow(data_used)), "0.25quant"]
    up50 <- fit$summary.fitted.values[seq_len(nrow(data_used)), "0.75quant"]
    cov50 <- mean(y_obs >= low50 & y_obs <= up50, na.rm = TRUE)
  } else {
    cov50 <- "—"
  }

  # warnings
  warns <- character()
  # overlap latent vs residual
  if(spatial_local && latent_global && is.finite(range_ratio) && range_ratio > 0.5 && range_ratio < 3) {
    warns <- c(warns,
      paste("⚠️  Insufficient spatial scale separation: latent/spatial range ratio = ",round(range_ratio, 2), ").\n",
             "   When 0.5 < ratio < 3, both fields may capture the same spatial structure, which causes identifiability issues and double-smoothing.\n",
             "   Recomended: latent.pcprior.range ≥ 3–5× residual range, or tighten spatial.pcprior.range.\n")
    )
  }  
  # weak identifiability (CI/median ratio)
  if(is.finite(max_CIratio) && max_CIratio > 25) {
    warns <- c(warns,
      "⚠️ Weak hyperparameter identifiability detected (CI/median > 25).\n",
      "   Recommended: use stronger PC priors.\n"
    )
  }
  # latent field dominating variance
  if(spatial_local && latent_global && is.finite(sigma_ratio) && sigma_ratio > 1.5) {
    warns <- c(warns,
      "⚠️ Variance imbalance between latent and residual/spatial SPDE. Variance ratio (sigma latent / sigma residual) = ", round(sigma_ratio, 2), ".\n",
      "   The latent field dominates the spatial variability, making the residual SPDE redundant.\n",
      "   Recommended: reduce latent.pcprior.sigma or increase spatial.pcprior.sigma.\n"
    )
  }
  # high correlation between latent and residual
  if(spatial_local && latent_global && is.finite(field_correlation) && field_correlation > 0.7) {
    warns <- c(warns,
      paste0("⚠️ High correlation between latent and residual spatial fields (r = ",
             round(field_correlation, 2), ").\n",
             "   Possible redundancy: both SPDE components capture the same spatial pattern.\n",
             "   Recommended: strengthen priors to separate scales, or remove one SPDE fields.\n")
    )
  }
  # residual autocorrelation
  if(is.finite(moran_I) && moran_I > 0.10) {
    warns <- c(warns,
      paste0("⚠️ Residual spatial autocorrelation detected (Moran’s I ≈ ", round(moran_I, 2), ").\n",
             "   Model missing local spatial structure.\n",
             "   Recommended: add a local SPDE component, refine mesh resolucion (smaller max.edges), or include missing covariates.\n")  #@@@JMB no estoy segura
    )
  }
  # posterior ≈ prior (weak data information)
  if(any(unlist(prior_close))) {
    warns <- c(warns,
      "⚠️ Posterior close to PC-prior mode: weak data information relative to prior strength.\n",
      "   Recommended: relax priors or increase data resolution.\n")    #@@@JMB rev recommendation??
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
          "   Recommended: check for outliers or review model specification/priors.\n\n"))
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
  res_mean <- mean(rs, na.rm = TRUE)
  df_r <- data.frame(resid = rs)
  center_label <- if(hierarchical) {
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
  if(spatial_local) {
    sp_vals <- fit$summary.random$spatial$mean
    sp_vals <- sp_vals[seq_len(nrow(data_used))]
  } else {
    sp_vals <- 0
  }
  if(latent_global) {
    lat_vals <- fit$summary.random$GLspde$mean
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
  if(spatial_local) {
    pE <- pE + ggplot2::geom_vline(xintercept = range_res_mean, linetype = "dashed", color = "#c0392b", linewidth = 0.5) +
      ggplot2::geom_text(
        data = data.frame(x = range_res_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Spatial range"),
        aes(x = x, y = y, label = label), 
        color = "#c0392b", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }
  if(latent_global) {
    pE <- pE + ggplot2::geom_vline(xintercept = range_lat_mean, linetype = "dashed", color = "#2980b9", linewidth = 0.5) +
      ggplot2::geom_text(
      data = data.frame(x = range_lat_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Latent range"),
      aes(x = x, y = y, label = label),
      color = "#2980b9", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }
  

  # covariates importance
  #pF <- nsbm_vars_importance(fit)


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
                       copy_beta = beta_IRegional), #copy_effect
    fixed_covariates = sig_vars,
    predictive = list(auc_full = auc_full,
                       brier = brier,
                       rmse = rmse, 
                       corr_obs_pred = pred_cor),
    calibration = list(slope = cal_slope,
                       ks_pit = ks_pit),
    coverage = list(cov50 = cov50,
                    cov95 = cov95),
    diagnostics = list(moran_I = moran_I,
                       max_CIratio = max_CIratio, # relative uncertainty
                       range_ratio = range_ratio, # scale separation ratio
                       sigma_ratio = sigma_ratio, # variance ratio
                       field_correlation = field_correlation), # prior influence
    warnings = warns,
    plots = list(hyperparams = pA,
                 spatialfields = pB,
                 correlogram = pC,
                 hist = pD1,
                 qq = pD2,
                 #vars_importance = pF,
                 semivariogram = pE)
  )

  return(out)
}






###----------------###




