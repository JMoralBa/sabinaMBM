#' @name NSBM.pure
#'
#' @title Nested species distribution modeling (Hierarchical Bayesian approach)
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling with spatial structure via \code{INLA} and \code{inlabru}. It allows the integration of global and regional scales with rigorous mathematical coupling structures.
#'
#' @param nsbm_obj An object of class \code{nsdm.vinput}, resulting from \code{sabinaNSDM::NSDM.SelectCovariates()}.
#' @param family A standard R \code{family} object (e.g. \code{binomial(link="logit")}), or the character string \code{"cp"} to fit a Cox point process (intensity). 
#' @param spde.mesh An \code{inla.mesh} object created externally with \code{create_mesh()}. Requiered if spatial and/or latent SPDE components are used. If `NULL` (default), the model runs without spatial structure. 
#' @param max.iter Numeric. Maximum number of Newton-Raphson iterations for the \code{inlabru} optimization (default: 10). Increasing this value helps convergence in complex models.
#' @param spde.pcprior.range Numeric vector length 2. PC-prior on the spatial range (e.g., \code{c(5, 0.01)}).
#' @param spde.pcprior.sigma Numeric vector length 2. PC-prior on the marginal standard deviation (e.g., \code{c(1, 0.01)}).
#' @param latent.pcprior.range Numeric vector length 2. PC-prior on the range of the latent global field SPDE. If \code{NULL}, no latent SPDE is created.
#' @param latent.pcprior.sigma Numeric vector length 2. PC-prior on the marginal standard deviation of the latent global field.
#' @param covariate.effects Optional named list to control the functional form of covariates (e.g., \code{"linear"} for linear, or \code{"rw2"} for non-linear splines)(see details). If \code{NULL} (default), all covariate effects remain constant (linear).
#' @param coupling.intercept Character. Controls how the regional intercept inherits information from the global intercept. Options:
#'   \itemize{
#'     \item \code{"unpooled"} (default): Independent intercepts (no borrowing of strength).
#'     \item \code{"ordered_hierarchical"}: Regional intercept borrows strength from the global (Zhou & Bradley, 2023).
#'     \item \code{"bayesian_feedback"}: Sequential updating using global posterior as regional prior via moment matching (mean and precision). A warning is issued if the global posterior is strongly skewed (|skewness| > 1), as moment matching may be unreliable in that case. See Figueira et al. (2024).
#'     \item \code{NULL}: Regional-only model (no global component).
#'   }
#' @param coupling.predictors Character or list. Controls how global and regional covariate effects are related. Options:
#'   \itemize{
#'     \item \code{"unpooled"} (default): Independent estimation of global and regional coefficients.
#'     \item \code{"ordered_hierarchical"}: Regional coefficient is modelled as a scaled copy of the global coefficient (\code{beta_RE | beta_GL ~ N(beta_GL, 0.5^2)}), implementing true hierarchical borrowing of strength via INLA \code{copy} (Zhou & Bradley, 2024). Requires the variable to be present at both scales.
#'     \item \code{"scale_decomposed"}: Shared covariates are decomposed into a large-scale trend (\code{X_ls = X_glo resampled}) and a small-scale anomaly (\code{X_ss = X_reg - X_glo}), following Cressie & Wikle (2011) and Zhou & Bradley (2024). In the output, coefficients are named \code{XGL_ls} and \code{XRE_ss}.
#'     \item \code{"bayesian_feedback"}: Sequential updating using global posteriors as regional priors.
#'     \item \code{NULL}: No coupling; effects are estimated independently.
#'   }
#' @param background.weights Controls how background points are weighted in the likelihood. Options:
#'   \itemize{
#'     \item \code{NULL} (default): no weights applied. Fast but produces a biased intercept when using presence-background data (Warton & Shepherd, 2010).
#'     \item \code{"area_weighted"}: Weights are computed automatically as \code{w = 1} for presences and \code{w = A / n_bg} for background points, where \code{A} is the total domain area (km2) and \code{n_bg} is the number of background points. This rigorously approximates a Poisson point process integration and corrects intercept bias (Renner et al., 2015).
#'     \item Numeric vector: user-supplied weights of length equal to the number of observations (presences + background) per scale.
#'   }
#'   Note: This argument is automatically ignored if \code{family = "cp"} since the Cox process intrinsically handles spatial integration via the continuous SPDE mesh.
#' @param proj.new.env Logical. Whether to compute predictions under new environmental scenarios (default: \code{TRUE}).
#' @param cv.folds Integer. Number of k-folds for spatial cross-validation (default: 1 = no CV).
#' @param n.threads Integer. Number of parallel threads to be used by INLA/inlabru (default: 1).
#' @param inla.int.strategy Character. INLA hyperparameter integration strategy. One of \code{"eb"} (Empirical Bayes, default), \code{"ccd"} (Central Composite Design), or \code{"grid"} (full grid integration). \code{"eb"} is fast but underestimates uncertainty by fixing hyperparameters at their posterior mode. \code{"ccd"} integrates over hyperparameters and is recommended for final/publication results. See Rue et al. (2009) and Simpson et al. (2017).
#' @param seed Optional integer to set the random seed for reproducibility.
#' @param save.output Logical. If \code{TRUE}, saves key model outputs to disk.
#'
#' @return A named list of class `nsbm.inlabru` with the following elements:
#' \item{Species.Name}{Species name}
#' \item{args}{List of arguments used in the model fitting.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{current.projections}{List with: current prediction (\code{pred}), residual spatial field (\code{pred_sp}), and latent global field (\code{pred_latent}).}
#' \item{marginals}{List with posterior marginals: \code{hyperpar} (SPDE hyperparameters), \code{random} (IGlobal, IRegional intercepts), \code{fixed} (covariate coefficients). Use with \code{plot(x, which = "hyperparams")} or \code{plot(x, which = "intercepts")}.}
#' \item{scale_params}{Named list with \code{mean} and \code{sd} used to standardize each covariate internally. Used for back-transforming marginals to original scale in plots.}
#' \item{new.projections}{List of projections to new.env (if `proj.new.env = TRUE`).}
#' \item{Summary}{Named list of \code{data.frame}s summarizing model metadata, model fit, hyperparameters, intercepts, fixed effects, predictive performance, calibration, and diagnostics.}  #@@@JMB revisar y refinar
#'
#' @details
#' family/link:
#' - \code{binomial(logit)}: Use for presence-absence (1/0) data. Output: occurrence probability.
#' - binomial(cloglog): Use for presence–absence (1/0) data when 1s are very rare (lots of 0s). Output: probability.
#' - poisson(log): Use for non-overdispersed counts (e.g. abundance counts 0,1,2…). Output: expected count
#' - \code{nbinomial(log)}: Use for overdispersed counts (variance > mean). Output: expected count.
#' - \code{gaussian(identity)}: Use for continuous normally distributed responses (e.g. log-transformed abundance). Output: predicted value on original scale.
#' - \code{gaussian(log)}: Use for strictly positive, right-skewed continuous data. Output: predicted value on original scale.
#' - beta(logit): Use for proportions in (0,1) (e.g. rescaled canopy cover). Output: predicted proportion
#' - \code{tweedie(log)}: Use for semicontinuous data with many zeros and a continuous positive tail (e.g. biomass). Output: expected value.
#' - cp: (Cox process) Use for presence-only data or spatial point patterns. Output: intensity.
#'
#' SPDE structure and priors:
#' - The mesh (`spde.mesh`) defines the domain for spatial and latent fields.  
#' - If neither spatial nor latent priors are defined, the model runs without SPDE components.
#' - Residual spatial field `spatial(geometry, model=matern)`: captures local spatial autocorrelation in the response not explained by covariates. Requires `spde.mesh`.
#' - Latent global field: captures broad-scale structure across regions and enters the linear predictor via \code{beta_GL * GLspde}.
#' - If both spatial and latent fields are used, the latent range should be ≥ 3× the spatial range to avoid overlap (Bakka et al., 2018).
#'
#' coupling.intercept = "ordered_hierarchical"
#' This option implements a true hierarchical Bayesian nested structure in which the regional intercept borrows strength from the global one. 
#' The regional intercept is expressed as a deviation from the global one using an INLA `copy` structure (IRegional = β * IGlobal + ε_regional), with a prior β ~ Normal(1, sd = 0.5), meaning the regional level is encouraged (but does not forced) to follow the global signal.
#' Ecologically, this formulation ensures that the global model captures broad-scale suitability (the species fundamental niche), while the regional intercept refines this large-scale signal/pattern to reflect local microclimatic conditions such and high-resolution predictors/conditions.
#' This hierarchical structure enables principled information sharing between scales while still allowing regional deviations.
#' 
#' covariate.effects: control per-covariate effects at global/regional scale.
#' - Option 1 `NULL` (default): no smoothing, All covariates enter linearly (`"linear"`) in both scales.
#' - Option 2 named `list`: provide a list with entries `global`, `regional`, and/or `default`.
#'   Structure:
#'     covariate.effects = list(
#'       global = list(var1 = "linear" | "drop" | list(model="rw2", u=..., alpha=...),
#'                     var2 = "linear", ...),
#'       regional = list(varA = "linear" | "drop" | list(model="rw2", u=..., alpha=...),
#'                       varB = "linear", ...),
#'       default = "linear" | "drop" | list(model="rw2", u=..., alpha=...)
#'     )
#'     Allowed values per covariate:
#'       - `"linear"`: linear effect (we use the term "linear" for consistency with inlabru).
#'       - `"drop"`: exclude that covariate.
#'       - `list(model="rw2", u=..., alpha=...)`: rw2 smoothing with a PC-prior on precision. Use the same units as the rasters for `u` and `alpha`.
#'     If a covariate is not listed under `global`/`regional`, it inherits from `default` (recommended `"linear"`).
#'     Example:
#'       covariate.effects = list(
#'         global = list(bio4 = "linear",
#'                       bio1 = list(model="rw2", u=..., alpha=...),
#'                       bio12 = "drop" ),
#'         regional = list(bio12 = "linear",
#'                         bio1 = list(model="rw2", u=..., alpha=...),
#'                         bio4 = "drop"),
#'         default = "linear")    # default rule if a covariate is not listed
#'
#' coupling.predictors: control of cross-scale covariate effects  
#' This argument defines how each predictor behaves across global and regional scales.  
#' It is only applied to variables that appear in both global and regional scales.
#' - Option 1 `NULL` (default): no coupling. Global and regional effects are estimated independently unless one is dropped via `covariate.effects`.
#' - Option 2 \code{"unpooled"}: global and regional effects enter the linear predictor as two independent components (η = β_GL * X_GL + β_RE * X_RE + …). No information is shared between scales.
#' - Option 3 \code{"ordered_hierarchical"}: implements a Bayesian soft constraint where the regional coefficient has a prior centred on 1 (\code{beta_RE ~ N(1, 0.5^2)}), encouraging — but not forcing — similarity with the global effect. With standardized covariates (applied internally), a coefficient of 1 implies the regional effect mirrors the global scale. True hierarchical copy over linear components is not reliably supported in INLA; this formulation is a practical approximation. Following Zhou & Bradley (2024).
#' - Option 4 named `list`: provide a list with entries `variables` and `default`.        #@@@JMB PARA PENSAR EN v2???? considerar permitir meter un SpatRaster de predicciones hechas fuera de aqui como spatial offset en los regional predictor.
#'   Structure:
#'     coupling.predictors = list(
#'       default   = "unpooled" | "NULL" | "ordered_hierarchical",
#'       variables = list(var1 = "...", var2 = "...", ...)
#'     )
#' Validity rules (automatically checked):
#' - Hierarchical coupling.predictors requires the predictor to exist at both scales.
#' - If `covariate.effects` drops the global effect of a variable, hierarchical coupling is not allowed (no parent effect available).
#' - If the regional effect is dropped, hierarchical coupling is also invalid.
#' - Global RW2 smoothing cannot be combined with hierarchical coupling for the same variable (avoids identifiability and double-smoothing issues).
#'
#' @seealso \code{\link{create_mesh}}, \code{\link{plot.nsbm.inlabru}}, \code{\link{summary.nsbm.inlabru}}
#'
#' @references
#' Cressie, N. & Wikle, C.K. (2011). \emph{Statistics for Spatio-Temporal Data}. Wiley.
#'
#' Zhou, S. & Bradley, J.R. (2024). Bayesian hierarchical modeling for bivariate multiscale
#' spatial data with application to blood test monitoring.
#' \emph{Spatial and Spatio-temporal Epidemiology}, 50, 100661.
#'
#' Bakka, H. et al. (2018). Spatial modelling with R-INLA: A review.
#' \emph{WIREs Computational Statistics}, 10, e1443.
#'
#' Warton, D.I. & Shepherd, L.C. (2010). Poisson point process models solve the
#' pseudo-absence problem for presence-only data in ecology.
#' \emph{The Annals of Applied Statistics}, 4(3), 1383--1402.
#'
#' @export
NSBM.pure <- function(nsbm_obj, 
                      family = binomial(link = "logit"), # family object binomial(), poisson(), etc., o "cp" para intensity (procesos puntuales)
                      covariate.effects = NULL,
                      coupling.intercept = "unpooled",
                      coupling.predictors = "unpooled",
                      spde.mesh = NULL,
                      spde.pcprior.range = NULL,
                      spde.pcprior.sigma = NULL,
                      latent.pcprior.range = NULL,
                      latent.pcprior.sigma = NULL,
                      background.weights = NULL,
                      proj.new.env = TRUE,
                      cv.folds = 1,
                      n.threads = 1,
                      max.iter = 10,
                      inla.int.strategy = "eb", 
                      seed = NULL,
                      save.output = FALSE) {

  vg <- nsbm_obj$Selected.Variables.Global
  vr <- nsbm_obj$Selected.Variables.Regional

  spatial_local <- !is.null(spde.mesh) && !is.null(spde.pcprior.range) && !is.null(spde.pcprior.sigma)
  latent_global <- !is.null(spde.mesh) && !is.null(latent.pcprior.range) && !is.null(latent.pcprior.sigma)

  ## Checks
  if(!inherits(nsbm_obj, "nsdm.vinput")) {
    stop("❌  The 'nsbm_obj' must be of class 'nsdm.vinput'. Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  #
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
                      poisson = "log", nbinomial = "log",
                      gaussian = c("identity", "log"), beta = "logit",
                      tweedie = "log", cp = NULL)
  if(!fam %in% names(valid_links)) {
    stop("❌ Unsupported family ", fam, ".\n",
         "  Supported families are: ", paste(names(valid_links), collapse = ", "), "\n")
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    stop("❌ Inalid link `", lnk, "` for this family `", fam, "`.\n",
         "  Allowed links for '", fam, "': ", paste(valid_links[[fam]], collapse = ", "), ".\n")
  }
  #
  valid_ic <- c("unpooled", "ordered_hierarchical", "bayesian_feedback")
  if(is.null(coupling.intercept)) {
    if(length(vr) == 0) stop("❌ 'coupling.intercept = NULL' (regional-only) but no regional covariates are present.\n\n")
    if(latent_global) stop("❌ Latent global SPDE not allowed in regional-only model (no global scale exists).\n\n")
    message("ℹ️ 'coupling.intercept = NULL'; running a regional-only model (no global component will be used).\n\n")
  } else {
    if(!is.character(coupling.intercept) || length(coupling.intercept) != 1 || !coupling.intercept %in% valid_ic) {
      stop("❌ The 'coupling.intercept' must be NULL, 'unpooled', 'ordered_hierarchical' or 'bayesian_feedback'.\n\n")
    }
    if(length(vg) == 0) {
      stop("❌ 'coupling.intercept = ", coupling.intercept, "' requires a global component,but none are present.\n",
           "   Use 'coupling.intercept = NULL' for a regional-only model.\n\n")
    }
  }
  if(is.null(spde.mesh) && (!is.null(spde.pcprior.range) || !is.null(latent.pcprior.range))) {
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
  #
  if(!is.null(covariate.effects)) {
    if(!is.list(covariate.effects)) stop("❌ `covariate.effects` must be a list or 'NULL'.\n")
    if(any(!names(covariate.effects) %in% c("global", "regional", "default"))) stop("❌ Invalid entries in `covariate.effects`. Allowed: 'global', 'regional', 'default'.\n\n")

    def <- covariate.effects$default
    if(!is.null(def)) {
      if(is.list(def)) {
        if(is.null(def$model) || !def$model %in% c("linear","rw2","drop")) stop("❌ `covariate.effects$default` must include `model = 'linear'|'rw2'|'drop'`.\n")
      } else if(is.character(def)) {
        if(!def %in% c("linear","drop")) stop("❌ `covariate.effects$default` must be 'linear' or 'drop'.\n")
      } else stop("❌ `covariate.effects$default` must be a character or list.\n")
    }
    # dropped vars
    all_dropped_gl <- !is.null(covariate.effects) && length(vg) > 0 && all(vapply(vg, function(x) .resolve_covariate_effects(x, "global", covariate.effects)$model == "drop", logical(1)))
    all_dropped_re <- !is.null(covariate.effects) && length(vr) > 0 && all(vapply(vr, function(x) .resolve_covariate_effects(x, "regional", covariate.effects)$model == "drop", logical(1)))

    if(all_dropped_gl && all_dropped_re) message("ℹ️ All covariates dropped at both global and regional scales.\n")
    else if(all_dropped_gl) message("ℹ️ All global covariates dropped.\n")
    else if(all_dropped_re) message("ℹ️ All regional covariates dropped.\n")

  }
  if(latent_global && length(vg) > 0) {
    if(any(vapply(vg, function(x) .resolve_covariate_effects(x, "global", covariate.effects)$model == "rw2", logical(1)))) {
      stop("❌ RW2 global effects are not allowed when the latent SPDE is active.\n",
           "   Both mechanisms smooth broad-scale structure and become redundant.\n",
           "   Disable RW2 for global covariates (in `covariate.effects`) or disable the latent SPDE\n",
           "   (`latent.pcprior.range` / `latent.pcprior.sigma` = NULL).\n\n")
    }
  }
  #
  valid_cp <- c("NULL", "unpooled", "ordered_hierarchical", "scale_decomposed", "bayesian_feedback")
  if(!is.null(coupling.predictors)) {
    if(is.character(coupling.predictors)) {
      if(!coupling.predictors %in% valid_cp) {
        stop("❌ `coupling.predictors` must be one of: ", paste(valid_cp, collapse=", "), ",\n",
             " or a named list specifying defaults and variable-specific modes.\n\n")
      }
    }
    else if(is.list(coupling.predictors)) {
      if(any(!names(coupling.predictors) %in% c("default", "variables"))) {
        stop("❌ Invalid entries in `coupling.predictors`. Allowed components: 'default', 'variables'.\n\n")
      }
      if(!is.null(coupling.predictors$default)) {
        if(!coupling.predictors$default %in% valid_cp) {
          stop("❌ `coupling.predictors$default` must be one of: ", paste(valid_cp, collapse=", "), ".\n\n")
        }
      }
      if(!is.null(coupling.predictors$variables)) {
        if(!is.list(coupling.predictors$variables)) {
          stop("❌ `coupling.predictors$variables` must be a named list.\n\n")
        }
        for(v in names(coupling.predictors$variables)) {
          mode_v <- coupling.predictors$variables[[v]]
          if(!mode_v %in% valid_cp) {
            stop("❌ Invalid coupling mode for variable '", v, "'. Allowed modes: ", paste(valid_cp, collapse=", "), ".\n\n")
          }
        }
      }
    }
    else {
      stop("❌ `coupling.predictors` must be NULL, a character mode, or a list(default=..., variables=list(...)).\n\n")
    }
  }
  # compatibility coupling.predictors x covariate.effects
  has_bf <- (!is.null(coupling.intercept) && coupling.intercept == "bayesian_feedback")
  has_joint <- (!is.null(coupling.intercept) && coupling.intercept %in% c("ordered_hierarchical", "scale_decomposed")) ||
    (length(vr) > 0 && any(vapply(vr, function(v) .resolve_coupling_predictor(v, coupling.predictors, vg) %in% c("ordered_hierarchical", "scale_decomposed"), logical(1))))
  vg_check <- if(is.null(coupling.intercept)) character(0) else vg
  if(length(vr) > 0) {
    for(v in vr) {
      cp_mode <- .resolve_coupling_predictor(v, coupling.predictors, vg_check)
      spec_re  <- .resolve_covariate_effects(v, "regional", covariate.effects)
      spec_gl  <- if(v %in% vg_check) .resolve_covariate_effects(v, "global", covariate.effects) else NULL
      # bayesian feedback
      if(cp_mode == "bayesian_feedback") {
        has_bf <- TRUE
        if(!(v %in% vg_check)) {
          stop("❌  `coupling.predictors`: variable '", v, "' explicitly requested 'bayesian_feedback', but it is missing in the global scale (or intercept is NULL).\n")
        }
        if(identical(spec_re$model, "rw2")) {
          stop("❌  `coupling.predictors`: variable '", v, "' cannot use 'bayesian_feedback' with 'rw2'. Supported for 'linear' only.\n")
        }
      # hierarchical joint methods
      } else if(cp_mode %in% c("ordered_hierarchical", "scale_decomposed")) {
        has_joint <- TRUE
        if(!(v %in% vg_check)) {
          stop("❌ `coupling.predictors`: variable '", v, "' cannot use '", cp_mode, "'. Missing in global scale.\n")
        }
        if(!is.null(spec_gl) && spec_gl$model == "drop") {
          stop("❌ `coupling.predictors`: variable '", v, "' cannot use '", cp_mode, "'. Global effect is dropped.\n")
        }
        if(spec_re$model == "drop") {
          stop("❌ `coupling.predictors`: variable '", v, "' cannot use '", cp_mode, "'. Regional effect is dropped.\n")
        }
        if(cp_mode == "ordered_hierarchical" && !is.null(spec_gl) && spec_gl$model == "rw2") {
          stop("❌ `coupling.predictors`: variable '", v, "' cannot use 'ordered_hierarchical'. Global effect uses RW2.\n")
        }
      }
    }
  }
  if(has_bf && has_joint) {
    stop("❌  Cannot mix 'bayesian_feedback' with 'ordered_hierarchical' or 'scale_decomposed'.\n",
         "    Feedback requires sequential fitting, while hierarchical couplings require joint fitting.\n\n")
  }
  #  
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1) stop("❌ 'seed' must be a single numeric value.\n\n")
    set.seed(seed)
  }
  #
  valid_int_strategy <- c("eb", "ccd", "grid")
  if(!inla.int.strategy %in% valid_int_strategy) {
    stop("❌ `inla.int.strategy` must be one of: ", paste(valid_int_strategy, collapse = ", "), ".\n\n")
  }
  if(inla.int.strategy == "eb") {
    warning("⚠️ `inla.int.strategy = 'eb'` (Empirical Bayes) is fast but underestimates uncertainty\n",
            "   by fixing hyperparameters at their posterior mode instead of integrating over them.\n",
            "   For final/publication results use `inla.int.strategy = 'ccd'`.\n",
            "   Reference: Rue et al. (2009) JRSS-B, Simpson et al. (2017) Statistical Science.\n\n")
  }
  #
  available_cores <- parallel::detectCores(logical = TRUE)
  if(!is.null(n.threads) && n.threads > available_cores) {
    stop(paste0("❌ Requested `n.threads` = ", n.threads, " exceeds available cores (", available_cores,").\n\n"))
  }
  INLA::inla.setOption(num.threads = n.threads)
  #
  if(!is.null(background.weights) && !is.numeric(background.weights)) {
    if(!is.character(background.weights) || length(background.weights) != 1 || background.weights != "area_weighted") {
      stop("❌ background.weights must be NULL, 'area_weighted', or a numeric vector.\n\n")
    }
  }
  if(!is.null(background.weights) && fam == "cp") {
    message("ℹ️ background.weights is ignored for family = 'cp': area integration is handled via continuous samplers.\n")
    background.weights <- NULL
  }


  ## Data preparation
  sp_covglo <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)
  crs <- sf::st_crs(sp_covglo)

    #@@@JMB ESTANDARIZACIÓN INTERNA??? var estandarizadas (media=0, sd=1) sobre el raster completo antes de entrar al modelo. 
    # Necesario para estabilidad INLA con copy en ordered_hierarchical porque daba problemas
    # Virgilio, estandarizar sobre raster completo o solo puntos de datos????
  # var standarization
  all_vars <- unique(c(names(sp_covglo), names(sp_covreg)))
  scale_params <- list()
  for(v in names(sp_covglo)) {
    vals <- terra::values(sp_covglo[[v]], na.rm = TRUE)
    scale_params[[v]] <- list(mean = mean(vals), sd = sd(vals))
    sp_covglo[[v]] <- (sp_covglo[[v]] - scale_params[[v]]$mean) / scale_params[[v]]$sd
  }
  for(v in names(sp_covreg)) {
    if(!v %in% names(scale_params)) {
      vals <- terra::values(sp_covreg[[v]], na.rm = TRUE)
      scale_params[[v]] <- list(mean = mean(vals), sd = sd(vals))
    }
    sp_covreg[[v]] <- (sp_covreg[[v]] - scale_params[[v]]$mean) / scale_params[[v]]$sd
  }

  # scale_decomposed pre-processing spatial anomalies 
  shared_vars <- intersect(vg, vr)
  
  # large-scale trend (_ls) and small-scale anomaly (_ss) — Zhou & Bradley (2024), Cressie & Wikle (2011)
  sd_vars <- shared_vars[vapply(shared_vars, function(X)
    .resolve_coupling_predictor(X, coupling.predictors, vg) == "scale_decomposed", logical(1))]

  if(length(sd_vars) > 0) {
    glo_resampled_stack <- terra::resample(sp_covglo[[sd_vars]], sp_covreg[[sd_vars[1]]])
    for(X in sd_vars) {
      sp_covreg[[paste0(X, "_ss")]] <- sp_covreg[[X]] - glo_resampled_stack[[X]]
      sp_covreg[[paste0(X, "_ls")]] <- glo_resampled_stack[[X]]
      for(sfx in c("_ls", "_ss")) {
        vname <- paste0(X, sfx)
        vals  <- terra::values(sp_covreg[[vname]], na.rm = TRUE)
        scale_params[[vname]] <- list(mean = mean(vals), sd = sd(vals))
        sp_covreg[[vname]] <- (sp_covreg[[vname]] - scale_params[[vname]]$mean) / scale_params[[vname]]$sd
      }
    }
  }
  
  pp_glo <- rbind(
    cbind(nsbm_obj$SpeciesData.XY.Global, resp = if(!is.null(nsbm_obj$Response.Global)) nsbm_obj$Response.Global else 1L),  #@@@JMB nsbm_obj$Response.Global y Regional habría que generarlos en sabinaNSDM input y arrastrar si hay algo
    cbind(nsbm_obj$Background.XY.Global, resp = 0L)     # si lo hacemo así poner algún check con stop/warning para que datos y family sean coherentes
  )
  pp_reg <- rbind(
    cbind(nsbm_obj$SpeciesData.XY.Regional, resp = if(!is.null(nsbm_obj$Response.Regional)) nsbm_obj$Response.Regional else 1L),  
    cbind(nsbm_obj$Background.XY.Regional, resp = 0L)          
  )

  pp_glo <- sf::st_as_sf(pp_glo, coords = c("x","y"), crs = crs)
  pp_glo <- sf::st_transform(pp_glo, crs)
  pp_reg <- sf::st_as_sf(pp_reg, coords = c("x","y"), crs = crs) 
  pp_reg <- sf::st_transform(pp_reg, crs)

  pp_glo$region <- 0L
  pp_reg$region <- 1L

  # rm NAs
  if(length(vg) > 0) {
    keep_glo <- stats::complete.cases(terra::extract(sp_covglo, pp_glo, ID = FALSE))
    pp_glo <- pp_glo[keep_glo, ]
  }
  if(length(vr) > 0) {
    ext_reg <- terra::extract(sp_covreg, pp_reg, ID = FALSE)
    if(has_joint) {
      ext_reg <- cbind(ext_reg, terra::extract(sp_covglo, pp_reg, ID = FALSE))
    }
    pp_reg <- pp_reg[stats::complete.cases(ext_reg), ]
  }

  ## Background weights
  w_glo <- NULL
  w_reg <- NULL
  if(is.character(background.weights) && background.weights == "area_weighted" && fam != "cp") {
    # global weights
    if(!is.null(coupling.intercept)) {
      n_pres_glo <- sum(pp_glo$resp != 0L)
      n_bg_glo <- sum(pp_glo$resp == 0L)
      A_glo <- as.numeric(terra::expanse(sp_covglo[[1]], unit = "km"))
      w_glo <- ifelse(pp_glo$resp != 0L, 1, A_glo / n_bg_glo)
    }
    # regional weights
    n_pres_reg <- sum(pp_reg$resp != 0L)
    n_bg_reg <- sum(pp_reg$resp == 0L)
    A_reg <- as.numeric(terra::expanse(sp_covreg[[1]], unit = "km"))
    w_reg <- ifelse(pp_reg$resp != 0L, 1, A_reg / n_bg_reg)

  } else if(is.numeric(background.weights)) {
    # user-supplied: split by scale
    n_glo_obs <- if(!is.null(coupling.intercept)) nrow(pp_glo) else 0L
    n_reg_obs <- nrow(pp_reg)
    if(length(background.weights) == n_glo_obs + n_reg_obs) {
      w_glo <- if(!is.null(coupling.intercept)) background.weights[seq_len(n_glo_obs)] else NULL
      w_reg <- background.weights[seq(n_glo_obs + 1L, n_glo_obs + n_reg_obs)]
    } else if(length(background.weights) == n_reg_obs) {
      w_reg <- background.weights
    } else {
      stop("❌ Length of background.weights (", length(background.weights), ") does not match ",
           "n_global + n_regional (", n_glo_obs + n_reg_obs, ") or n_regional (", n_reg_obs, ").\n\n")
    }
  }

  ## SPDE domain definition
  pred_sf <- sf::st_as_sf(terra::as.points(sp_covreg, values = FALSE))
  #sf::st_crs(pred_sf) <- crs  # !!!
  #pred_sf <- sf::st_set_crs(pred_sf, crs)
  pred_sf <- sf::st_transform(pred_sf, crs)

  #rm NAs
  ext_pred <- terra::extract(sp_covreg, pred_sf, ID = FALSE)
  if(has_joint) {
    ext_pred <- cbind(ext_pred, terra::extract(sp_covglo, pred_sf, ID = FALSE))
  }
  pred_sf <- pred_sf[stats::complete.cases(ext_pred), ]

  if(!is.null(coupling.intercept)) {
    # Use global raster bbox as domain — more conservative than convex hull over points,
    # avoids excluding raster cells that lack nearby presence/background records.
    bdy_glo <- sf::st_as_sfc(sf::st_bbox(sp_covglo))
    sf::st_crs(bdy_glo) <- crs
  } else {
    bdy_glo <- NULL
  }
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(terra::as.polygons(!is.na(sp_covreg[[1]]), dissolve = TRUE))))
  pred_sf$region <- 1L

  ## SPDE components
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


  ## Nested intercept
  # prior for intercepts iid. param = c(1, 0.01) -> P(σ > 1) = 0.01
  IID_PRIOR <- "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"
  # prior for beta_copy (regional deviation). beta_regional ~ Normal(mean = 1, sd = 0.5). INLA uses precisión, so sd = 0.5 -> tau = 1/(0.5^2) = 4
  COPY_BETA_PRIOR <- "hyper = list(beta = list(prior = 'normal', param = c(1, 4)))"

  if(is.null(coupling.intercept)) {
    intercept_terms <- c(paste0("IRegional(1, model='iid', ", IID_PRIOR, ")"))
  } else if(coupling.intercept == "unpooled") {
    intercept_terms <- c(
      paste0("IGlobal(1, model='iid', ", IID_PRIOR, ")"),
      paste0("IRegional(1, model='iid', ", IID_PRIOR, ")")
    )
  } else if(coupling.intercept == 'ordered_hierarchical') {
    intercept_terms <- c(
      paste0("IGlobal(1, model='iid', ", IID_PRIOR, ")"),
      paste0("IRegional(1, copy='IGlobal', fixed=FALSE, ", COPY_BETA_PRIOR, ")"))
  } else if(coupling.intercept == 'bayesian_feedback') {
    intercept_terms <- c(
      paste0("IGlobal(1, model='iid', ", IID_PRIOR, ")"),
      paste0("IRegional(main = 1, model='linear', mean.linear = bf_mean_int, prec.linear = bf_prec_int)")
    )
  }
  base_intercepts <- paste(intercept_terms, collapse = " + ")


  ## Model components
  nsbm_obj_fcov <- nsbm_obj
  if(is.null(coupling.intercept)) {
    nsbm_obj_fcov$Selected.Variables.Global <- character(0)
  }
  cmp_cov <- .fcov(obj = nsbm_obj_fcov, 
                  spobjglo = "sp_covglo", 
                  spobjreg = "sp_covreg",
                  sp_covglo = sp_covglo,
                  sp_covreg = sp_covreg,
                  covariate.effects = covariate.effects,
                  coupling.predictors = coupling.predictors,
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
    paste0("IGlobal", f_spatial, f_latent, .opt_plus(cmp_cov$like$fglobal))
  }
  rhs_reg <- paste0("IRegional", f_spatial, f_latent, .opt_plus(cmp_cov$like$fregional))

  # likelihoods
  liks <- .build_likelihoods(fam, lnk, rhs_glo, rhs_reg,
                              pp_glo, pp_reg, bdy_glo, bdy_reg, dom,
                              coupling.intercept,
                              w_glo = w_glo, w_reg = w_reg)
  lik_glo <- liks$lik_glo
  lik_reg <- liks$lik_reg

  eta_terms <- c(
    "IRegional",
    #if(!is.null(coupling.intercept) && coupling.intercept == "unpooled") "IGlobal" else NULL,
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


  ## Assemble likelihoods and fit model
  lik_list <- list()
  if(!is.null(coupling.intercept)) lik_list <- c(lik_list, list(lik_glo))
  lik_list <- c(lik_list, list(lik_reg))

  # if coupling.intercept = NULL, desactivate bayesian feedback
  needs_feedback <- !is.null(coupling.intercept) && (
    coupling.intercept == "bayesian_feedback" ||
    any(vapply(vr, function(v) .resolve_coupling_predictor(v, coupling.predictors, vg) == "bayesian_feedback", logical(1)))
  )
  
  fit <- .fit_nsbm(
    cmp = cmp, 
    lik_list = lik_list, 
    coupling.intercept = coupling.intercept, 
    coupling.predictors = coupling.predictors,
    needs_feedback = needs_feedback,
    vr = vr,
    n.threads = n.threads,
    seed = seed,
    int.strategy = inla.int.strategy,
    max.iter = max.iter
  )

  ## K-fold CV (only fam no cp)
  if(cv.folds > 1) {      #@@@JMB la cv actual usa k-folds random. Deberíamos cosiderar cV espacial bloqueado con blockCV::cv_spatial()??? 
    if(fam == "cp") {
      warning("⚠️ Cross-validation (cv.folds > 1) is not implemented for family = 'cp'. CV results will be NULL.")
      cv_res <- NULL
    } else {
      # stratified k folds
      folds_g <- if(!is.null(coupling.intercept)) .make_stratified_kfolds(pp_glo$resp, cv.folds) else NULL
      folds_r <- .make_stratified_kfolds(pp_reg$resp, cv.folds)
      cv_metrics <- numeric(cv.folds)

      both_cov <- c(sp_covglo, sp_covreg)

      for(k in seq_len(cv.folds)) {
        # build likelihoods for fold k
        train_g <- if(!is.null(coupling.intercept)) pp_glo[folds_g != k, ] else NULL
        train_r <- pp_reg[folds_r != k, ]
        test_r  <- pp_reg[folds_r == k, ]

        # recalculate weights for this k fold's training data
        w_glo_k <- NULL
        w_reg_k <- NULL
        if(is.character(background.weights) && background.weights == "area_weighted" && fam != "cp") {
          if(!is.null(coupling.intercept) && !is.null(train_g)) {
            n_bg_glo_k <- sum(train_g$resp == 0L)
            A_glo_k    <- as.numeric(terra::expanse(sp_covglo[[1]], unit = "km"))
            w_glo_k    <- ifelse(train_g$resp != 0L, 1, A_glo_k / n_bg_glo_k)
          }
          n_bg_reg_k <- sum(train_r$resp == 0L)
          A_reg_k    <- as.numeric(terra::expanse(sp_covreg[[1]], unit = "km"))
          w_reg_k    <- ifelse(train_r$resp != 0L, 1, A_reg_k / n_bg_reg_k)
        } else if(is.numeric(background.weights)) {
          # subset user weights to training indices
          w_glo_k <- if(!is.null(coupling.intercept) && !is.null(w_glo)) w_glo[folds_g != k] else NULL
          w_reg_k <- if(!is.null(w_reg)) w_reg[folds_r != k] else NULL
        }
        liks_k <- .build_likelihoods(fam, lnk, rhs_glo, rhs_reg,
                                      train_g, train_r,
                                      bdy_glo, bdy_reg, dom,
                                      coupling.intercept,
                                      w_glo = w_glo_k, w_reg = w_reg_k)
        lik_list_k <- Filter(Negate(is.null), list(liks_k$lik_glo, liks_k$lik_reg))

        # fit fold k
        fit_k <- .fit_nsbm(
          cmp = cmp,
          lik_list = lik_list_k,
          coupling.intercept = coupling.intercept,
          coupling.predictors = coupling.predictors,
          needs_feedback = needs_feedback,
          vr = vr,
          n.threads = n.threads,
          seed = seed,
          int.strategy = inla.int.strategy
        )

        # rm NAs from test set
        coords_t <- sf::st_coordinates(test_r)
        keep <- stats::complete.cases(terra::extract(both_cov, coords_t))

        # pred and metric
        test_r2 <- test_r[keep, ]
        if(nrow(test_r2) > 0) {
          pk <- predict(fit_k, test_r2, pred_formula)
          if(fam %in% c("binomial", "beta")) {
            if(length(unique(test_r2$resp)) > 1) {
              cv_metrics[k] <- suppressMessages(as.numeric(pROC::auc(test_r2$resp, pk$mean)))
            } else {
              cv_metrics[k] <- NA_real_
            }
          } else {
            cv_metrics[k] <- sqrt(mean((test_r2$resp - pk$mean)^2, na.rm = TRUE))
          }
        } else {
          cv_metrics[k] <- NA_real_
        }
      }

      metric_name <- if(fam %in% c("binomial", "beta")) "AUC" else "RMSE"
      cv_res <- list(
        cv.folds    = cv.folds,
        metric_name = metric_name,
        metric_mean = mean(cv_metrics, na.rm = TRUE),
        metric_sd   = sd(cv_metrics, na.rm = TRUE)
      )
    }
  } else {
    cv_res <- NULL
  }


  ## Predictions
  pred_combined <- predict(fit, pred_sf, pred_formula)
  pred <- .pred_as_tif(pred_combined, sp_covreg)
  pred_sp <- if(spatial_local) {
    .pred_as_tif(predict(fit, pred_sf, ~ spatial), sp_covreg)
  } else NULL
  pred_lat <- if(latent_global) {
    .pred_as_tif(predict(fit, pred_sf, ~ beta_GL * GLspde), sp_covreg)
  } else NULL


  ## New scenarios
  proj_list <- list()
  if(proj.new.env && !is.null(nsbm_obj$Scenarios)) {
    for(sc in names(nsbm_obj$Scenarios)) {
      #sc <- 1
      scen_rast <- terra::unwrap(nsbm_obj$Scenarios[[sc]])
      scen_df <- sf::st_as_sf(terra::as.points(scen_rast, values = FALSE))
      scen_df <- sf::st_transform(scen_df, crs)
      scen_df$region <- 1L
      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[sc]] <- .pred_as_tif(proj_pred, template = scen_rast)  # from sf to tif
    }
  }


  ## Diagnostics
  coords_reg <- sf::st_coordinates(pp_reg)
  data_used  <- data.frame(x = coords_reg[,1], y = coords_reg[,2], resp = pp_reg$resp)
  diag_block <- .nsbm_diagnostics(
                  fit,
                  fam,  
                  data_used,
                  n_glo = if(!is.null(coupling.intercept)) nrow(pp_glo) else 0L,
                  priors = list(spde.pcprior.range = spde.pcprior.range,
                                spde.pcprior.sigma = spde.pcprior.sigma,
                                latent.pcprior.range = latent.pcprior.range,
                                latent.pcprior.sigma = latent.pcprior.sigma),
                  pred_sp = if(spatial_local) pred_sp else NULL,
                  pred_lat = if(latent_global) pred_lat else NULL,
                  coupling.intercept = coupling.intercept,
                  scale_params = scale_params)

  if(length(diag_block$warnings)) {
    message(paste(unique(diag_block$warnings), collapse = "\n\n"))
  }


  ## Save outputs
  species <- nsbm_obj$Species.Name
  if(save.output) {
    # directories
    values_path <- file.path("Results", "NSBM_pure", "Values")
    projections_path <- file.path("Results", "NSBM_pure", "Projections")
    fs::dir_create(values_path, recurse = TRUE)
    fs::dir_create(projections_path, recurse = TRUE)
    # fixed effects
    if(!is.null(fit$summary.fixed)) {
      write.csv(fit$summary.fixed, file = file.path(values_path, paste0(species, "_fixed_effects.csv")), row.names = TRUE)
    }
    # spatial random effects
    if(!is.null(fit$summary.random)) {
      for(ran in names(fit$summary.random)) {
        write.csv(fit$summary.random[[ran]], file = file.path(values_path, paste0(species, "_random_", ran, ".csv")), row.names = TRUE)
      }
    }
    # hypermarams
    if(!is.null(fit$summary.hyperpar)) {
      write.csv(fit$summary.hyperpar, file = file.path(values_path, paste0(species, "_hyperparameters.csv")), row.names = TRUE)
    }
    # evaluation metrics
    eval_metrics <- data.frame(
      Metric = c("WAIC", "DIC", "MLik", "LCPO (sum log-CPO)"),
      Value = c(fit$waic$waic, fit$dic$dic, fit$mlik[1], diag_block$bayes_fit$lcpo_val))
    write.csv(eval_metrics, file = file.path(values_path, paste0(species, "_evaluation.csv")), row.names = FALSE)
    # CPO values (one per observation)   #@@@JMB useful for leave-one-out diagnostics or model comparison??
    write.csv(data.frame(CPO = fit$cpo$cpo), file = file.path(values_path, paste0(species, "_pointwise_CPO.csv")), row.names = FALSE)
    # full model object (fit)
    saveRDS(fit, file = file.path(values_path, paste0(species, "_model_fit.rds")))
    # pred current
    if(!is.null(pred)) {
      file_path <- file.path(projections_path, paste0(species, "_Current.tif"))
      terra::writeRaster(terra::unwrap(pred), file_path, overwrite = TRUE)   
    }
    # pred_sf
    if(!is.null(pred_sp)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_Spatial.tif"))
      terra::writeRaster(terra::unwrap(pred_sp), file_path, overwrite = TRUE)
    }
    # latent contribution
    if(!is.null(pred_lat)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_latent.tif"))
      terra::writeRaster(terra::unwrap(pred_lat), file_path, overwrite = TRUE)
    }
    # new scenarios
    if(length(proj_list) > 0 && !is.null(nsbm_obj$Scenarios)) {
      for(i in seq_along(proj_list)) {
        sc_name <- names(proj_list)[i]
        file_path <- file.path(projections_path, paste0(species, "_", sc_name, ".tif"))
        terra::writeRaster(terra::unwrap(proj_list[[i]]), file_path, overwrite = TRUE)
      }
    }
    # diagnostic plots
    plot_specs <- list(
      list(key = "hyperparams", file = "_hyperparams.png", w = 6, h = 5),
      list(key = "correlogram", file = "_correlogram.png", w = 6, h = 5),
      list(key = "hist", file = "_residualHistogram.png", w = 6, h = 5),
      list(key = "qq", file = "_QQplot.png", w = 6, h = 4),
      list(key = "spatialfields", file = "_SPDEfields.png", w = 6, h = 5),
      list(key = "semivariogram", file = "_semivariogram.png", w = 6, h = 4)
    )
    for(ps in plot_specs) {
      p <- diag_block$plots[[ps$key]]
      if(!is.null(p)) {
        ggplot2::ggsave(
          filename = file.path(values_path, paste0(species, ps$file)),
          plot = p, width = ps$w, height = ps$h, dpi = 300, bg = "white")
      }
    }
    # predictive metrics
    pred_metrics <- data.frame(
      Metric = c("AUC", "Tjur_R2", "Brier", "RMSE", "Cor_obs_pred",
                 "Calibration_slope", "Coverage_95", "MLPD", "PIT_KS_pvalue"),
      Value = c(diag_block$predictive$auc_full,
                diag_block$predictive$tjur_r2,
                diag_block$predictive$brier,
                diag_block$predictive$rmse,
                diag_block$predictive$corr_obs_pred,
                diag_block$calibration$slope,
                diag_block$coverage$cov95,
                diag_block$bayes_fit$mlpd_val,
                diag_block$calibration$ks_pit),
      stringsAsFactors = FALSE)
    write.csv(pred_metrics,
              file = file.path(values_path, paste0(species, "_Predictive_metrics.csv")),
              row.names = FALSE)


     message("ℹ️ Results saved in the following local folder(s):")
     message(paste(
     "  - Projections (current pred, spatial field pred_sp, latent field pred_lat and scenarios: ", projections_path, "\n",
     " - Model values (fixed/random effects, hyperparameters, evaluation and diagnostics: ", values_path, "\n",
     " - Full model object (.rds): ", file.path(values_path, paste0(species, "_model_fit.rds")), "\n"
     ))
  }


  # summary          #@@@JMB pendiente revisar/completar...
  summary_df <- .nsbm_generate_summary(
    fit = fit, 
    species = species, 
    fam = fam, 
    lnk = lnk, 
    coupling.intercept = coupling.intercept, 
    coupling.predictors = coupling.predictors, 
    diag_block = diag_block, 
    cv_res = cv_res,
    vg = vg,
    vr = vr
  )


  # return
  sabina <- list(
    Species.Name = species,
    args = list(
      family = fam,
      link = lnk,
      #spde.mesh = !is.null(spde.mesh),
      spde.pcprior.range = spde.pcprior.range,
      spde.pcprior.sigma = spde.pcprior.sigma,
      latent.pcprior.range = latent.pcprior.range,
      latent.pcprior.sigma = latent.pcprior.sigma,
      coupling.intercept = coupling.intercept,
      coupling.predictors = coupling.predictors,
      covariate.effects = covariate.effects,
      background.weights = background.weights,
      proj.new.env = proj.new.env,
      cv.folds = cv.folds,
      n.threads = n.threads,
      inla.int.strategy = inla.int.strategy,     
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
    marginals = list(
      hyperpar = fit$marginals.hyperpar,
      random = fit$marginals.random,
      fixed = fit$marginals.fixed
    ),
    scale_params = scale_params,
    pit_values = diag_block$calibration$pit_values,
    diagnostic_plots = list(
      correlogram  = diag_block$plots$correlogram,
      semivariogram = diag_block$plots$semivariogram
    ),
    Summary = summary_df
  )
 
  attr(sabina, "class") <- "nsbm.inlabru"
  return(sabina)

}


###----------------###




