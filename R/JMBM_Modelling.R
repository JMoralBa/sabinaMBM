#' @name JMBM.Modelling
#'
#' @title Nested species distribution modeling (Hierarchical Bayesian approach)
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling with spatial structure via \code{INLA} and \code{inlabru}. It allows the integration of global and regional scales with rigorous mathematical coupling structures.
#'
#' @param jmbm_obj An object of class \code{nsdm.vinput}, resulting from \code{sabinaNSDM::NSDM.SelectCovariates()}.
#' @param family A standard R \code{family} object (e.g. \code{binomial(link="logit")}), or the character string \code{"cp"} to fit a Cox point process (intensity). 
#' @param spde.mesh An \code{inla.mesh} object created externally with \code{create_mesh()}. Required if any spatial field (S_shared or S_re) is used. If \code{NULL} (default), the model runs without spatial structure.#' @param regional.pcprior.range Numeric vector length 2. PC-prior on the range of the local residual field S_re (e.g., \code{c(2, 0.01)}). If \code{NULL} and \code{shared.pcprior.range} is also \code{NULL}, no spatial fields are used.
#' @param shared.pcprior.range Numeric vector of length 2. PC-prior for the range of the shared spatial field S_shared, e.g. \code{c(10, 0.01)} means P(range < 10) = 0.01. Must be substantially larger than \code{regional.pcprior.range} to ensure scale separation (Bakka et al., 2018). If \code{NULL} (and \code{regional.pcprior.range} is also \code{NULL}), no spatial fields are used.
#' @param shared.pcprior.sigma Numeric vector of length 2. PC-prior for the marginal standard deviation of S_shared, e.g. \code{c(1, 0.01)} means P(sigma > 1) = 0.01. For multi-scale separation, set \code{shared.pcprior.sigma[1]} substantially higher than \code{regional.pcprior.sigma[1]} (e.g., c(1, 0.01) vs c(0.5, 0.01)).
#' @param regional.pcprior.range Numeric vector of length 2. PC-prior for the range of the local residual field S_re (regional predictor only). If \code{NULL} (and \code{shared.pcprior.range} is also \code{NULL}), no spatial fields are used. Providing \code{shared.pcprior.range} without \code{regional.pcprior.range} raises an error.
#' @param regional.pcprior.sigma Numeric vector of length 2. PC-prior for the marginal standard deviation of S_re.
#' @param covariate.effects Optional named list controlling the functional form of covariate effects. If \code{NULL} (default), all effects are linear. Accepts entries \code{"global"}, \code{"regional"}, and/or \code{"default"}, each set to \code{"linear"}, \code{"drop"}, \code{"rw2"} (non-linear RW2 with default PC prior \code{u = 5}, \code{alpha = 0.01}), or \code{list(model = "rw2", u = ..., alpha = ...)} for custom PC prior parameters.
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
#'     \item \code{"scale_decomposed"}: Shared covariates are decomposed into a large-scale trend (\code{X_glo_res = X_glo resampled}) and a small-scale anomaly (\code{X_reg_anom = X_reg - X_glo}), following Cressie & Wikle (2011) and Zhou & Bradley (2024). In the output, coefficients are named \code{XGL_glo_res} and \code{XRE_reg_anom}.
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
#' @return A named list of class `jmbm.inlabru` with the following elements:
#' \item{Species.Name}{Species name}
#' \item{args}{List of arguments used in the model fitting.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{current.projections}{List with: current prediction (\code{pred}), local residual spatial field (\code{pred_Sre}, if S_re was specified), and shared broad-scale spatial field (\code{pred_Sshared}, if S_shared was specified).}
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
#' - The mesh (\code{spde.mesh}) discretizes the spatial domain for the SPDE approximation (Lindgren et al., 2011).
#' - Spatial architecture is determined by which priors are supplied:
#'   - Both \code{NULL}: no spatial fields; purely covariate-driven model.
#'   - Only \code{regional.pcprior.*}: S_re only, in the regional predictor.
#'   - Both specified: S_shared (enters both global and regional predictors with weight 1) + S_re (regional only).
#'   - Only \code{shared.pcprior.*}: not allowed (raises an error).
#' - S_shared is a broad-scale Matern GRF estimated jointly from both likelihoods, capturing shared biogeographic autocorrelation not explained by covariates. It enters both \eqn{\eta_{GL}} and \eqn{\eta_{RE}} with fixed weight 1 (no scaling coefficient).
#' - S_re is a fine-scale Matern GRF in the regional predictor only, capturing residual local autocorrelation.
#' - To avoid double-smoothing, \code{shared.pcprior.range[1]} should be at least 3x \code{regional.pcprior.range[1]} (Bakka et al., 2018).
#' - Both fields use Penalised Complexity (PC) priors (Simpson et al., 2017; Fuglstad et al., 2019), which shrink toward a structureless base model unless data provide evidence of spatial structure.
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
#' bayesian_feedback for new scenarios:
#' When using \code{coupling.predictors = "bayesian_feedback"} with \code{proj.new.env = TRUE},
#' note that the regional priors are fixed based on the present-day global posteriors.
#' They are NOT updated for future scenarios. This is a current limitation.
#' For dynamic re-estimation of priors under different climate conditions, use
#' \code{coupling.predictors = "ordered_hierarchical"} (soft constraint) or
#' \code{coupling.predictors = "scale_decomposed"} (macro/micro decomposition) instead.
#'
#' @seealso \code{\link{create_mesh}}, \code{\link{plot.jmbm.inlabru}}, \code{\link{summary.jmbm.inlabru}}
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
JMBM.Modelling <- function(jmbm_obj, 
                      family = binomial(link = "logit"), # family object binomial(), poisson(), etc., o "cp" para intensity (procesos puntuales)
                      covariate.effects = NULL,
                      coupling.intercept = "unpooled",
                      coupling.predictors = "unpooled",
                      spde.mesh = NULL,
                      regional.pcprior.range = NULL,
                      regional.pcprior.sigma = NULL,
                      shared.pcprior.range = NULL,
                      shared.pcprior.sigma = NULL,
                      background.weights = NULL,
                      proj.new.env = TRUE,
                      cv.folds = 1,
                      n.threads = 1,
                      inla.int.strategy = "eb", 
                      seed = NULL,
                      save.output = FALSE,
                      verbose = TRUE) {

if (!is.null(jmbm_obj$Selected.Variables.Global) && length(jmbm_obj$Selected.Variables.Global) > 0) {
  vg <- sort(jmbm_obj$Selected.Variables.Global)
  shared_vr <- intersect(vg, jmbm_obj$Selected.Variables.Regional)
  exclusive_vr <- setdiff(jmbm_obj$Selected.Variables.Regional, vg)
  vr <- c(shared_vr, sort(exclusive_vr))
} else {
  vg <- character(0)
  shared_vr <- character(0)
  vr <- sort(jmbm_obj$Selected.Variables.Regional)
}

  has_Sre <- !is.null(spde.mesh) && !is.null(regional.pcprior.range)  && !is.null(regional.pcprior.sigma)
  has_Sshared <- !is.null(spde.mesh) && !is.null(shared.pcprior.range) && !is.null(shared.pcprior.sigma)


  .info("sabinaJMBM: Joint Multiscale Species Distribution Model")

  ## Checks
  if(!inherits(jmbm_obj, "nsdm.vinput")) {
    .stop("The 'jmbm_obj' must be of class 'nsdm.vinput'. Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  #
  if(inherits(family, "family")) {
    fam <- family$family
    lnk <- family$link
  } else if(is.character(family) && length(family) == 1 && family == "cp") {
    fam <- "cp"
    lnk <- NULL
  } else {
    .stop("`family` must be either a standard family() object or the string 'cp'.")  #@@@JMB poner permitidos o enviar a ?JMBM.Modelling details?
  }
  valid_links <- list(binomial = c("logit", "cloglog"),
                      poisson = "log", nbinomial = "log",
                      gaussian = c("identity", "log"), beta = "logit",
                      tweedie = "log", cp = NULL)
  if(!fam %in% names(valid_links)) {
    .stop(paste0("Unsupported family: ", fam, ". Supported: ", paste(names(valid_links), collapse = ", ")))
  }
  if(!is.null(lnk) && !lnk %in% valid_links[[fam]]) {
    .stop(paste0("Invalid link '", lnk, "' for family '", fam, "'. Allowed: ", paste(valid_links[[fam]], collapse = ", ")))
  }
  #
  valid_ic <- c("unpooled", "ordered_hierarchical", "bayesian_feedback")
  if(is.null(coupling.intercept)) {
    if(length(vr) == 0) .stop("'coupling.intercept = NULL' (regional-only) but no regional covariates are present.")
    if(has_Sshared) .stop("S_shared (shared spatial field) not allowed in regional-only model.")
  } else {
    if(!is.character(coupling.intercept) || length(coupling.intercept) != 1 || !coupling.intercept %in% valid_ic) {
      .stop("The 'coupling.intercept' must be NULL, 'unpooled', 'ordered_hierarchical' or 'bayesian_feedback'.")
    }
    if(length(vg) == 0) {
      .stop(paste0("'coupling.intercept = ", coupling.intercept, "' requires a global component. Use NULL for regional-only."))
    }
  }
  if(is.null(spde.mesh) && (!is.null(regional.pcprior.range) || !is.null(shared.pcprior.range))) {
    .stop("SPDE priors were provided but 'spde.mesh' is NULL. Create a mesh with create_mesh().")
  }
  if(!is.null(spde.mesh) && !inherits(spde.mesh, "inla.mesh")) {
    .stop("`spde.mesh` must be a valid INLA mesh object. Use create_mesh() to build it.")
  }
  if(has_Sre && (is.null(regional.pcprior.range) || is.null(regional.pcprior.sigma))) {
    .stop("Missing priors for Sre field. Provide both regional.pcprior.range and sigma.")
  }
  if(has_Sshared && (is.null(shared.pcprior.range) || is.null(shared.pcprior.sigma))) {
    .stop("Missing priors for Sshared field. Provide both shared.pcprior.range and sigma.")
  }
  if(has_Sre && has_Sshared) {
    if(shared.pcprior.range[1] < regional.pcprior.range[1] * 3) {
      .warn("`shared.pcprior.range` < 3x `regional.pcprior.range`: fields may overlap, causing variance inflation.")
    }
  }
  #
  if(!is.null(covariate.effects)) {
    if(!is.list(covariate.effects)) .stop("`covariate.effects` must be a list or 'NULL'.")
    if (any(!names(covariate.effects) %in% c("global", "regional", "default"))) {
      .stop("Invalid entries in `covariate.effects`. Allowed: 'global', 'regional', 'default'.")
    }

    def <- covariate.effects$default
    if(!is.null(def)) {
      if(is.list(def)) {
        if (is.null(def$model) || !def$model %in% c("linear", "rw2", "drop")) {
          .stop("`covariate.effects$default` must include `model = 'linear'|'rw2'|'drop'`.")
        }
      } else if(is.character(def)) {
        if (!def %in% c("linear", "drop")) .stop("`covariate.effects$default` must be 'linear' or 'drop'.")
      } else {
        .stop("`covariate.effects$default` must be a character or list.")
      }
    }
    # dropped vars
    all_dropped_gl <- !is.null(covariate.effects) && length(vg) > 0 && all(vapply(vg, function(x) .resolve_covariate_effects(x, "global", covariate.effects)$model == "drop", logical(1)))
    all_dropped_re <- !is.null(covariate.effects) && length(vr) > 0 && all(vapply(vr, function(x) .resolve_covariate_effects(x, "regional", covariate.effects)$model == "drop", logical(1)))
    if (all_dropped_gl && all_dropped_re) {
      .info("All covariates dropped at both global and regional scales.")
    } else if (all_dropped_gl) {
      .info("All global covariates dropped.")
    } else if (all_dropped_re) {
      .info("All regional covariates dropped.")
    }
  }
  if(has_Sshared && length(vg) > 0) {
    if(any(vapply(vg, function(x) .resolve_covariate_effects(x, "global", covariate.effects)$model == "rw2", logical(1)))) {
      .stop("RW2 global effects are not allowed when the Sshared SPDE is active. Both mechanisms smooth broad-scale structure. Disable RW2 or deactivate Sshared.")
    }
  }
  #
  valid_cp <- c("NULL", "unpooled", "ordered_hierarchical", "scale_decomposed", "bayesian_feedback")
  if(!is.null(coupling.predictors)) {
    if(is.character(coupling.predictors)) {
      if(!coupling.predictors %in% valid_cp) {
        .stop(paste0("`coupling.predictors` must be one of: ", paste(valid_cp, collapse = ", ")))
      }
    }
    else if(is.list(coupling.predictors)) {
      if(any(!names(coupling.predictors) %in% c("default", "variables"))) {
       .stop("Invalid entries in `coupling.predictors`. Allowed: 'default', 'variables'.")
      }
      if(!is.null(coupling.predictors$default)) {
        if(!coupling.predictors$default %in% valid_cp) {
          .stop(paste0("`coupling.predictors$default` must be one of: ", paste(valid_cp, collapse = ", ")))
        }
      }
      if(!is.null(coupling.predictors$variables)) {
        if (!is.list(coupling.predictors$variables)) .stop("`coupling.predictors$variables` must be a named list.")
        for(v in names(coupling.predictors$variables)) {
          mode_v <- coupling.predictors$variables[[v]]
          if(!mode_v %in% valid_cp) {
            .stop(paste0("Invalid coupling mode for variable '", v, "'. Allowed: ", paste(valid_cp, collapse = ", ")))
          }
        }
      }
    }
    else {
      .stop("`coupling.predictors` must be NULL, a character mode, or a list.")
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
          .stop(paste0("Variable '", v, "' requested 'bayesian_feedback' but is missing at global scale."))
        }
        if(identical(spec_re$model, "rw2")) {
          .stop(paste0("Variable '", v, "' cannot use 'bayesian_feedback' with 'rw2'. Supported for 'linear' only."))
        }
      # hierarchical joint methods
      } else if(cp_mode %in% c("ordered_hierarchical", "scale_decomposed")) {
        has_joint <- TRUE
        if(!(v %in% vg_check)) {
          .stop(paste0("Variable '", v, "' cannot use '", cp_mode, "'. Missing at global scale."))
        }
        if(!is.null(spec_gl) && spec_gl$model == "drop") {
          .stop(paste0("Variable '", v, "' cannot use '", cp_mode, "'. Global effect is dropped."))
        }
        if(spec_re$model == "drop") {
          .stop(paste0("Variable '", v, "' cannot use '", cp_mode, "'. Regional effect is dropped."))
        }
        if(cp_mode == "ordered_hierarchical" && !is.null(spec_gl) && spec_gl$model == "rw2") {
          .stop(paste0("Variable '", v, "' cannot use 'ordered_hierarchical'. Global effect uses RW2."))
        }
      }
    }
  }
  if(has_bf && has_joint) {
    .stop("Cannot mix 'bayesian_feedback' with 'ordered_hierarchical' or 'scale_decomposed'.")
  }
  #
    if(!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) .stop("'seed' must be a single numeric value.")
    set.seed(seed)
  }
  #
  valid_int_strategy <- c("eb", "ccd", "grid")
  if(!inla.int.strategy %in% valid_int_strategy) {
    .stop(paste0("`inla.int.strategy` must be one of: ", paste(valid_int_strategy, collapse = ", ")))
  }
  if(inla.int.strategy == "eb") {
    .warn("'inla.int.strategy = eb' (Empirical Bayes) is fast but underestimates uncertainty. Use 'ccd' for final results.")
  }
  #
  available_cores <- parallel::detectCores(logical = TRUE)
  if(!is.null(n.threads) && n.threads > available_cores) {
    .stop(paste0("Requested `n.threads` (", n.threads, ") exceeds available cores (", available_cores, ")."))
  }
  INLA::inla.setOption(num.threads = n.threads)
  #
  if(!is.null(background.weights) && !is.numeric(background.weights)) {
    if(!is.character(background.weights) || length(background.weights) != 1 || background.weights != "area_weighted") {
      .stop("`background.weights` must be NULL, 'area_weighted', or a numeric vector.")
    }
  }
  if(!is.null(background.weights) && fam == "cp") {
    .info("background.weights is ignored for family = 'cp' (handled via continuous samplers).")
    background.weights <- NULL
  }

  if(is.null(coupling.intercept)) {
    .check("Architecture: Regional-only (no global component)")
  } else {
    cp_pred_str <- if(is.list(coupling.predictors)) "variable-specific" else if(is.null(coupling.predictors)) "none" else coupling.predictors
    .check(paste0("Architecture: Joint model (Intercepts: ", coupling.intercept, " | Predictors: ", cp_pred_str, ")"))
  }
  fam_str <- if(fam == "cp") {
    "log-Gaussian Cox process (LGCP)"
  } else {
    paste0(fam, " (link: ", lnk, ")")
  }
  .check(paste0("Family: ", fam_str))


  ## Data preparation
  sp_covglo <- terra::unwrap(jmbm_obj$IndVar.Global.Selected)
  sp_covreg <- terra::unwrap(jmbm_obj$IndVar.Regional.Selected)
  crs <- sf::st_crs(sp_covglo)

  all_model_vars <- unique(c(names(sp_covglo), names(sp_covreg)))
  scale_params_glo <- list()
  scale_params_reg <- list()

  # scale decomposed pre-processing
  sd_vars <- shared_vr[vapply(shared_vr, function(X) .resolve_coupling_predictor(X, coupling.predictors, vg) == "scale_decomposed", logical(1))]
  if(length(sd_vars) > 0) {
    glo_resampled_stack <- terra::resample(sp_covglo[[sd_vars]], sp_covreg[[sd_vars[1]]])
    for(X in sd_vars) {
      sp_covreg[[paste0(X, "_reg_anom")]] <- sp_covreg[[X]] - glo_resampled_stack[[X]]
      sp_covreg[[paste0(X, "_glo_res")]] <- glo_resampled_stack[[X]]
    }
  }

  if(length(all_model_vars) > 0) {
    .info("Pre-processing environmental covariates...")
    if(length(sd_vars) > 0) .check("Applying scale decomposition (macro-trend vs micro-anomaly)")
    .check("Standardizing covariates (Z-score)")
  }

  n_cores_std <- min(n.threads, max(length(names(sp_covglo)), length(names(sp_covreg))))

  # global standarization 
  if (!is.null(coupling.intercept) && length(names(sp_covglo)) > 0) {
    result_glo <- .standardize_rasters(sp_covglo, n_cores = n.threads)
    sp_covglo <- result_glo$rast
    scale_params_glo <- result_glo$params
  }

  # regional standarization
  if (length(names(sp_covreg)) > 0) {
    result_reg <- .standardize_rasters(sp_covreg, n_cores = n_cores_std)
    sp_covreg <- result_reg$rast
    scale_params_reg <- result_reg$params
  }

  scale_params <- c(
    scale_params_glo,
    scale_params_reg[setdiff(names(scale_params_reg), names(scale_params_glo))]
  )

  # log stats
  if (length(names(sp_covreg)) > 0 && verbose) {
    for (v in names(sp_covreg)) {
      if (v %in% sd_vars) next
      # standardized values
      m_val <- terra::global(sp_covreg[[v]], "mean", na.rm = TRUE)[1, 1]
      s_val <- terra::global(sp_covreg[[v]], "sd", na.rm = TRUE)[1, 1]
      # original values
      original_mean <- scale_params_reg[[v]][["mean"]]
      original_sd <- scale_params_reg[[v]][["sd"]]
      .item(sprintf("%s (regional): Z-mean = %+.3f, Z-sd = %.3f | Original: mean = %+.2f, sd = %.2f", 
                    v, m_val, s_val, original_mean, original_sd))
    }
  }

  if (!is.null(coupling.intercept) && !is.null(sp_covglo) && length(names(sp_covglo)) > 0 && verbose) {
    for (v in names(sp_covglo)) {
      # standardized values
      m_val <- terra::global(sp_covglo[[v]], "mean", na.rm = TRUE)[1, 1]
      s_val <- terra::global(sp_covglo[[v]], "sd", na.rm = TRUE)[1, 1]
      # original values
      original_mean <- scale_params_glo[[v]][["mean"]]
      original_sd   <- scale_params_glo[[v]][["sd"]]
      .item(sprintf("%s (global): Z-mean = %+.3f, Z-sd = %.3f | Original: mean = %+.2f, sd = %.2f", 
                    v, m_val, s_val, original_mean, original_sd))
    }
  }

  #
  bg_or_abs_glo <- if(is.null(jmbm_obj$Absences.XY.Global)) jmbm_obj$Background.XY.Global else jmbm_obj$Absences.XY.Global
  bg_or_abs_reg <- if(is.null(jmbm_obj$Absences.XY.Regional)) jmbm_obj$Background.XY.Regional else jmbm_obj$Absences.XY.Regional

  pp_glo <- rbind(
    cbind(jmbm_obj$SpeciesData.XY.Global, resp = if(!is.null(jmbm_obj$Response.Global)) jmbm_obj$Response.Global else 1L),  #@@@JMB jmbm_obj$Response.Global y Regional habría que generarlos en sabinaNSDM input y arrastrar si hay algo
    cbind(jmbm_obj$bg_or_abs_glo, resp = 0L)     # si lo hacemo así poner algún check con stop/warning para que datos y family sean coherentes
  )
  pp_reg <- rbind(
    cbind(jmbm_obj$SpeciesData.XY.Regional, resp = if(!is.null(jmbm_obj$Response.Regional)) jmbm_obj$Response.Regional else 1L),  
    cbind(jmbm_obj$bg_or_abs_reg, resp = 0L)          
  )

  pp_glo <- sf::st_as_sf(pp_glo, coords = c("x","y"), crs = crs)
  pp_glo <- sf::st_transform(pp_glo, crs)
  pp_reg <- sf::st_as_sf(pp_reg, coords = c("x","y"), crs = crs) 
  pp_reg <- sf::st_transform(pp_reg, crs)

  pp_glo$region <- 0L
  pp_reg$region <- 1L

  # rm NAs
  if(length(vg) > 0) {
    keep_glo <- stats::complete.cases(terra::extract(sp_covglo, terra::vect(pp_glo), ID = FALSE))
    pp_glo <- pp_glo[keep_glo, ]
  }
  if(length(vr) > 0) {
    ext_reg <- terra::extract(sp_covreg, terra::vect(pp_reg), ID = FALSE)
    if(has_joint) {
      ext_reg <- cbind(ext_reg, terra::extract(sp_covglo, terra::vect(pp_reg), ID = FALSE))
    }
    pp_reg <- pp_reg[stats::complete.cases(ext_reg), ]
  }

  ## Background weights   #@@@JMB!! "area_weighted" sigue Warton & Shepherd 2010/Renner etal 2015 pero esta gente no contemplan dos escalas. Con area global >> area reg los pesos hay que equilibrarlos o algo así? ¿como hacemos esto para dos escalas?
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
      .stop(paste0("Length of background.weights (", length(background.weights), ") does not match observations."))
    }
  }

  ## SPDE domain definition
  pred_sf <- sf::st_as_sf(terra::as.points(sp_covreg, values = TRUE))
  pred_sf <- sf::st_transform(pred_sf, crs)

  #rm NAs
  ext_pred <- terra::extract(sp_covreg, pred_sf, ID = FALSE)
  if(has_joint) {
    ext_pred <- cbind(ext_pred, terra::extract(sp_covglo, pred_sf, ID = FALSE))
  }
  pred_sf <- pred_sf[stats::complete.cases(ext_pred), ]

  if(!is.null(coupling.intercept)) {
    # Use global raster bbox as domain. more conservative than convex hull over points,
    # avoids excluding raster cells that lack nearby presence/background records.
    bdy_glo <- sf::st_as_sfc(sf::st_bbox(sp_covglo))
    sf::st_crs(bdy_glo) <- crs
  } else {
    bdy_glo <- NULL
  }
  bdy_reg <- sf::st_union(sf::st_make_valid(sf::st_as_sf(terra::as.polygons(!is.na(sp_covreg[[1]]), dissolve = TRUE))))
  pred_sf$region <- 1L

  ## SPDE components
  if (has_Sre || has_Sshared) {
    .info("Building latent spatial fields (SPDE)...")
    if (has_Sre) .check("S_re (fine-scale residual spatial field) activated")
    if (has_Sshared) .check("S_shared (broad-scale spatial field) activated")
  }

  if(has_Sre)  {
    matern_loc <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = regional.pcprior.range,
      prior.sigma = regional.pcprior.sigma
    )
  }

  if(has_Sshared) {
    matern_shared <- INLA::inla.spde2.pcmatern(
      mesh = spde.mesh,
      prior.range = shared.pcprior.range,
      prior.sigma = shared.pcprior.sigma)
  } else {
    matern_shared <- NULL
  }


  ## Nested intercept
  # prior for intercepts iid. param = c(1, 0.01) -> P(σ > 1) = 0.01
  IID_PRIOR <- "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"
  #@@@JMB!! prior for beta_copy (regional deviation). beta_regional ~ Normal(mean = 1, sd = 0.5). INLA uses precisión, so sd = 0.5 -> tau = 1/(0.5^2) = 4  ¿está bien??
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
  jmbm_obj_fcov <- jmbm_obj
  if(is.null(coupling.intercept)) {
    jmbm_obj_fcov$Selected.Variables.Global <- character(0)
  }
  cmp_cov <- .fcov(obj = jmbm_obj_fcov,
                  spobjglo = "sp_covglo",
                  spobjreg = "sp_covreg",
                  sp_covglo = sp_covglo,
                  sp_covreg = sp_covreg,
                  covariate.effects = covariate.effects,
                  coupling.predictors = coupling.predictors,
                  pp_glo_sf = pp_glo,
                  pp_reg_sf = pp_reg)

  # component formula (intercept + spataial + Sshared + fcovs
  cmp <- c(
    base_intercepts,
    if(has_Sre) "Sre(geometry, model = matern_loc)" else NULL,
    if(has_Sshared) "Sshared(main = geometry, model = matern_shared)" else NULL,
    cmp_cov$cmp   # bio1GL(), bio1RE(),...
  )
  cmp <- paste(cmp[!is.na(cmp) & nzchar(cmp)], collapse = " + ")
  cmp <- as.formula(paste("~", cmp))

  dom <- if (has_Sre || has_Sshared) list(geometry = spde.mesh) else NULL
  f_Sre <- if(has_Sre) " + Sre" else ""
  f_Sshared <- if(has_Sshared) " + Sshared" else ""

  .opt_plus <- function(s) if(!is.null(s) && nzchar(s)) paste0(" + ", s) else ""
 
  rhs_glo <- if(is.null(coupling.intercept)) {
    NULL  # no for regional-only
  } else {
    paste0("IGlobal", f_Sshared, .opt_plus(cmp_cov$like$fglobal))
  }
  rhs_reg <- paste0("IRegional", f_Sshared, f_Sre, .opt_plus(cmp_cov$like$fregional))

  # likelihoods
  liks <- .build_likelihoods(fam, lnk, rhs_glo, rhs_reg,
                              pp_glo, pp_reg, bdy_glo, bdy_reg, dom,
                              coupling.intercept,
                              w_glo = w_glo, w_reg = w_reg)
  lik_glo <- liks$lik_glo
  lik_reg <- liks$lik_reg

  eta_terms <- c(
    "IRegional",
    if(has_Sshared) "Sshared" else NULL,
    if(has_Sre) "Sre" else NULL,
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
  
  .info(sprintf("Fitting Bayesian model (Integration strategy: '%s')...", inla.int.strategy))

  fit <- .fit_jmbm(
    cmp = cmp, 
    lik_list = lik_list, 
    coupling.intercept = coupling.intercept, 
    coupling.predictors = coupling.predictors,
    needs_feedback = needs_feedback,
    vr = vr,
    n.threads = n.threads,
    seed = seed,
    int.strategy = inla.int.strategy
  )


  ## K-fold CV (only fam no cp)
  if(cv.folds > 1) {      #@@@JMB!! la cv actual usa k-folds random. Deberíamos cosiderar spatial block CV??? 
    .info(sprintf("Performing spatial cross-validation (K = %d folds)...", cv.folds))
    if(fam == "cp") {
      .warn("Cross-validation (cv.folds > 1) is not implemented for family = 'cp'. CV results will be NULL.")
      cv_res <- NULL
    } else {
      # stratified k folds
      folds_g <- if(!is.null(coupling.intercept)) .make_stratified_kfolds(pp_glo$resp, cv.folds) else NULL
      folds_r <- .make_stratified_kfolds(pp_reg$resp, cv.folds)

      n_cores_cv <- min(n.threads, cv.folds)
      inla_threads_k <- max(1, floor(n.threads / n_cores_cv))

      if(n_cores_cv > 1) {
        furrr::plan(furrr::multisession, workers = n_cores_cv, quiet = TRUE)
      }

      A_glo_total <- if(is.character(background.weights) && background.weights == "area_weighted" && fam != "cp" && !is.null(coupling.intercept)) as.numeric(terra::expanse(sp_covglo[[1]], unit = "km")) else NULL
      A_reg_total <- if(is.character(background.weights) && background.weights == "area_weighted" && fam != "cp") as.numeric(terra::expanse(sp_covreg[[1]], unit = "km")) else NULL

      both_cov <- c(sp_covglo, sp_covreg)
      coords_all_r <- sf::st_coordinates(pp_reg)
      valid_pts_mask <- stats::complete.cases(terra::extract(both_cov, coords_all_r, ID = FALSE))

     cv_worker <- function(k) {
        train_g <- if(!is.null(coupling.intercept)) pp_glo[folds_g != k, ] else NULL
        train_r <- pp_reg[folds_r != k, ]
        test_r  <- pp_reg[folds_r == k, ]

        w_glo_k <- NULL
        w_reg_k <- NULL
        
        if(is.character(background.weights) && background.weights == "area_weighted" && fam != "cp") {
          if(!is.null(coupling.intercept) && !is.null(train_g)) {
            n_bg_glo_k <- sum(train_g$resp == 0L)
            w_glo_k <- ifelse(train_g$resp != 0L, 1, A_glo_total / n_bg_glo_k)
          }
          n_bg_reg_k <- sum(train_r$resp == 0L)
          w_reg_k <- ifelse(train_r$resp != 0L, 1, A_reg_total / n_bg_reg_k)
        } else if(is.numeric(background.weights)) {
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
        fit_k <- .fit_jmbm(
          cmp = cmp,
          lik_list = lik_list_k,
          coupling.intercept = coupling.intercept,
          coupling.predictors = coupling.predictors,
          needs_feedback = needs_feedback,
          vr = vr,
          n.threads = inla_threads_k,
          seed = seed,
          int.strategy = inla.int.strategy
        )

        # rm NAs
        keep_k <- valid_pts_mask[folds_r == k]
        test_r2 <- test_r[keep_k, ]

        if(nrow(test_r2) > 0) {
          pk <- predict(fit_k, test_r2, pred_formula)
          if(fam %in% c("binomial", "beta")) {
            if(length(unique(test_r2$resp)) > 1) {
              suppressMessages(as.numeric(pROC::auc(test_r2$resp, pk$mean)))
            } else { NA_real_ }
          } else {
            sqrt(mean((test_r2$resp - pk$mean)^2, na.rm = TRUE))
          }
        } else {
          NA_real_
        }
      }

      if(n_cores_cv > 1) {
        cv_metrics_list <- furrr::future_lapply(seq_len(cv.folds), cv_worker, .options = furrr::furrr_options(seed = TRUE))
        furrr::plan(furrr::sequential) # Reseteo de seguridad al terminar
      } else {
        cv_metrics_list <- lapply(seq_len(cv.folds), cv_worker)
      }
     
      cv_metrics <- unlist(cv_metrics_list)
      metric_name <- if(fam %in% c("binomial", "beta")) "AUC" else "RMSE"
   
      cv_res <- list(
        cv.folds = cv.folds,
        metric_name = metric_name,
        metric_mean = mean(cv_metrics, na.rm = TRUE),
        metric_sd = sd(cv_metrics, na.rm = TRUE)
      )
    }
  } else {
    cv_res <- NULL
  }


  ## Predictions
  .info("Generating spatial predictions...")
  .check("Current suitability map")

  pred_combined <- predict(fit, pred_sf, pred_formula)
  pred <- .pred_as_tif(pred_combined, sp_covreg)
  pred_Sre <- if(has_Sre) {
    .pred_as_tif(predict(fit, pred_sf, ~ Sre), sp_covreg)
  } else NULL
  pred_Sshared <- if(has_Sshared) {
    .pred_as_tif(predict(fit, pred_sf, ~ Sshared), sp_covreg)
  } else NULL

  rm(pred_combined, pred_sf)
  gc(verbose = FALSE)


  ## New scenarios
  proj_list <- list()
  if(proj.new.env && !is.null(jmbm_obj$Scenarios)) {
    .check(sprintf("Projecting to new scenarios: %s", paste(names(jmbm_obj$Scenarios), collapse = ", ")))

    sp_covglo_curr <- sp_covglo
    sp_covreg_curr <- sp_covreg
   
    for(sc in names(jmbm_obj$Scenarios)) {
      scen_rast <- terra::unwrap(jmbm_obj$Scenarios[[sc]])
           
      if(!identical(sf::st_crs(crs), terra::crs(scen_rast))) {
        scen_rast <- terra::project(scen_rast, terra::crs(sp_covreg_curr))
      }

      sp_covglo_fut <- sp_covglo_curr
      sp_covreg_fut <- sp_covreg_curr      

      vars_glo <- intersect(names(sp_covglo_fut), names(scen_rast))
      if(length(vars_glo) > 0) {
        scen_glo_resampled <- terra::resample(scen_rast[[vars_glo]], sp_covglo_fut)
      }
      
      vars_reg <- setdiff(names(sp_covreg_fut), grep("_glo_res$|_reg_anom$", names(sp_covreg_fut), value = TRUE))
      vars_reg <- intersect(vars_reg, names(scen_rast))
      if(length(vars_reg) > 0) {
        scen_reg_resampled <- terra::resample(scen_rast[[vars_reg]], sp_covreg_fut)
      }

      # global standarization fut
      if (!is.null(coupling.intercept) && length(vars_glo) > 0) {
        for(v in names(sp_covglo_fut)) {
          if(v %in% names(scen_rast)) {
            fut_layer <- scen_glo_resampled[[v]]
            std_fut_layer <- (fut_layer - scale_params[[v]]$mean) / scale_params[[v]]$sd
          
            if(!(v %in% names(sp_covreg_fut))) {
              m_pres <- mean(terra::values(sp_covglo_curr[[v]]), na.rm = TRUE)
              m_fut  <- mean(terra::values(std_fut_layer), na.rm = TRUE)
            }
            sp_covglo_fut[[v]] <- std_fut_layer
          } 
        }
      }

      # regional standarization fut
      for(v in names(sp_covreg_fut)) {
        if(grepl("_glo_res$|_reg_anom$", v)) next 
        
        if(v %in% names(scen_rast)) {
          fut_layer <- scen_reg_resampled[[v]]

          if (!is.null(scale_params[[v]])) {
            std_fut_layer <- (fut_layer - scale_params[[v]]$mean) / scale_params[[v]]$sd
         
            if(!(v %in% names(sp_covreg_fut))) {
              m_pres <- mean(terra::values(sp_covreg_curr[[v]]), na.rm = TRUE)
              m_fut  <- mean(terra::values(std_fut_layer), na.rm = TRUE)
            }
            sp_covreg_fut[[v]] <- std_fut_layer
          }
        } 
      }

      # scale decomposed fut
      if(length(sd_vars) > 0) {
        vars_sd <- intersect(sd_vars, names(scen_rast))
        
        if(length(vars_sd) > 0) {
          raw_glo_base <- scen_glo_resampled[[vars_sd]]
          glo_resampled_scenario <- terra::resample(raw_glo_base, sp_covreg_curr)
          
          for(X in vars_sd) {
            raw_reg <- scen_reg_resampled[[X]]
            raw_glo <- glo_resampled_scenario[[X]]
            
            # Standarization with historical parameters
            std_reg_anom <- ((raw_reg - raw_glo) - scale_params[[paste0(X, "_reg_anom")]]$mean) / scale_params[[paste0(X, "_reg_anom")]]$sd
            std_glo_res <- (raw_glo - scale_params[[paste0(X, "_glo_res")]]$mean) / scale_params[[paste0(X, "_glo_res")]]$sd
            
            scen_rast[[paste0(X, "_reg_anom")]] <- std_reg_anom
            scen_rast[[paste0(X, "_glo_res")]] <- std_glo_res
            
            sp_covreg_fut[[paste0(X, "_reg_anom")]] <- std_reg_anom
            sp_covreg_fut[[paste0(X, "_glo_res")]] <- std_glo_res

            shift_reg_anom <- mean(terra::values(std_reg_anom), na.rm=TRUE) - mean(terra::values(sp_covreg_curr[[paste0(X, "_reg_anom")]]), na.rm=TRUE)
            shift_glo_res <- mean(terra::values(std_glo_res), na.rm=TRUE) - mean(terra::values(sp_covreg_curr[[paste0(X, "_glo_res")]]), na.rm=TRUE)
          }
        }
      }
           
      sp_covglo <- sp_covglo_fut
      sp_covreg <- sp_covreg_fut
      
      scen_pts <- terra::as.points(scen_rast, values = TRUE, na.rm = TRUE)
      scen_df <- sf::st_as_sf(scen_pts)
      scen_df <- sf::st_transform(scen_df, crs)
      scen_df$region <- 1L

      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[sc]] <- .pred_as_tif(proj_pred, template = scen_rast)

      rm(scen_pts, scen_df, proj_pred)
      gc(verbose = FALSE)
    }
    
    sp_covglo <- sp_covglo_curr
    sp_covreg <- sp_covreg_curr
  }


  ## Diagnostics
  .info("Computing model diagnostics...")

  coords_reg <- sf::st_coordinates(pp_reg)
  data_used  <- data.frame(x = coords_reg[,1], y = coords_reg[,2], resp = pp_reg$resp)
  diag_block <- .jmbm_diagnostics(
                  fit,
                  fam,  
                  data_used,
                  n_glo = if(!is.null(coupling.intercept)) nrow(pp_glo) else 0L,
                  priors = list(regional.pcprior.range = regional.pcprior.range,
                                regional.pcprior.sigma = regional.pcprior.sigma,
                                shared.pcprior.range = shared.pcprior.range,
                                shared.pcprior.sigma = shared.pcprior.sigma),
                  pred_Sre = if(has_Sre) pred_Sre else NULL,
                  pred_Sshared = if(has_Sshared) pred_Sshared else NULL,
                  coupling.intercept = coupling.intercept,
                  scale_params_glo = scale_params_glo,
                  scale_params_reg = scale_params_reg)


  ## Save outputs
  species <- jmbm_obj$Species.Name
  if(save.output) {
    # directories
    values_path <- file.path("Results", "MBM", "Values")
    projections_path <- file.path("Results", "MBM", "Projections")
    fs::dir_create(values_path, recurse = TRUE)
    fs::dir_create(projections_path, recurse = TRUE)
    # fixed effects
    if(!is.null(fit$summary.fixed)) {
      write.csv(fit$summary.fixed, file = file.path(values_path, paste0(species, "_fixed_effects.csv")), row.names = TRUE)
    }
    # Sre random effects
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
    # CPO values (one per observation)   #@@@JMB!! useful for leave-one-out diagnostics or model comparison??
    write.csv(data.frame(CPO = fit$cpo$cpo), file = file.path(values_path, paste0(species, "_pointwise_CPO.csv")), row.names = FALSE)
    # full model object (fit)
    saveRDS(fit, file = file.path(values_path, paste0(species, "_model_fit.rds")))
    # pred current
    if(!is.null(pred)) {
      file_path <- file.path(projections_path, paste0(species, "_Current.tif"))
      terra::writeRaster(terra::unwrap(pred), file_path, overwrite = TRUE)   
    }
    # pred_sf
    if(!is.null(pred_Sre)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_Sre.tif"))
      terra::writeRaster(terra::unwrap(pred_Sre), file_path, overwrite = TRUE)
    }
    # Sshared contribution
    if(!is.null(pred_Sshared)) {
      file_path <- file.path(projections_path, paste0(species, "_Current_Sshared.tif"))
      terra::writeRaster(terra::unwrap(pred_Sshared), file_path, overwrite = TRUE)
    }
    # new scenarios
    if(length(proj_list) > 0 && !is.null(jmbm_obj$Scenarios)) {
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
      list(key = "Srefields", file = "_SPDEfields.png", w = 6, h = 5),
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


    .info("Results saved in the following local folder(s):")
    .item(paste0("Projections: ", projections_path))
    .item(paste0("Model values: ", values_path))
    .item(paste0("Full model object (.rds): ", file.path(values_path, paste0(species, "_model_fit.rds"))))
  }


  # summary          #@@@JMB pendiente revisar/completar...
  summary_df <- .jmbm_generate_summary(
    fit = fit, 
    species = species, 
    fam = fam, 
    lnk = lnk, 
    coupling.intercept = coupling.intercept, 
    coupling.predictors = coupling.predictors, 
    diag_block = diag_block, 
    cv_res = cv_res,
    vg = vg,
    vr = vr,
    scale_params = scale_params,
    scale_params_glo = scale_params_glo,
    scale_params_reg = scale_params_reg,
    has_spatial = !is.null(spde.mesh)
  )


  # return
  sabina <- list(
    Species.Name = species,
    args = list(
      family = fam,
      link = lnk,
      #spde.mesh = !is.null(spde.mesh),
      regional.pcprior.range = regional.pcprior.range,
      regional.pcprior.sigma = regional.pcprior.sigma,
      shared.pcprior.range = shared.pcprior.range,
      shared.pcprior.sigma = shared.pcprior.sigma,
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
    Selected.Variables.Global = jmbm_obj$Selected.Variables.Global,
    Selected.Variables.Regional = jmbm_obj$Selected.Variables.Regional,
    current.projections = list(
      pred = terra::wrap(pred),
      pred_Sre = if(!is.null(pred_Sre)) terra::wrap(pred_Sre) else NULL,
      pred_Sshared = if(!is.null(pred_Sshared)) terra::wrap(pred_Sshared) else NULL
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
 
  attr(sabina, "class") <- "jmbm.inlabru"

  if(verbose) message("\n  ✓ Model fitted successfully. Check 'summary()' for evaluation metrics.\n")

  return(sabina)

}


###----------------###




