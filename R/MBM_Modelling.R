#' @name MBM.Modelling
#'
#' @title Joint multiscale Bayesian species distribution modelling
#'
#' @description Fits a joint multiscale Bayesian species distribution model using spatial random fields via \code{INLA} and \code{inlabru}. Integrates global and regional scales through flexible coupling architectures with formal uncertainty propagation..
#'
#' @param jmbm_obj An object of class \code{nsdm.vinput}, resulting from \code{sabinaNSDM::NSDM.SelectCovariates()}.
#' @param spde.mesh An \code{inla.mesh} object created externally with \code{create_mesh()}. Required if any spatial field (S_shared or S_re) is used. If \code{NULL} (default), the model runs without spatial structure.
#' @param shared.pcprior.range Numeric vector of length 2. PC-prior for the range of the shared spatial field S_shared, e.g. \code{c(10, 0.01)} means P(range < 10) = 0.01. Must be substantially larger than \code{regional.pcprior.range} to ensure scale separation (Bakka et al., 2018). If \code{NULL} (and \code{regional.pcprior.range} is also \code{NULL}), no spatial fields are used.
#' @param shared.pcprior.sigma Numeric vector of length 2. PC-prior for the marginal standard deviation of S_shared, e.g. \code{c(1, 0.01)} means P(sigma > 1) = 0.01. For multi-scale separation, set \code{shared.pcprior.sigma[1]} substantially higher than \code{regional.pcprior.sigma[1]} (e.g., c(1, 0.01) vs c(0.5, 0.01)).
#' @param regional.pcprior.range Numeric vector of length 2. PC-prior for the range of the local residual field S_re (regional predictor only). If \code{NULL} (and \code{shared.pcprior.range} is also \code{NULL}), no spatial fields are used. Providing \code{shared.pcprior.range} without \code{regional.pcprior.range} raises an error.
#' @param regional.pcprior.sigma Numeric vector of length 2. PC-prior for the marginal standard deviation of S_re.
#' @param covariate.effects Optional named list controlling the functional form of covariate effects. If \code{NULL} (default), all effects are linear. Accepts entries \code{"global"}, \code{"regional"}, and/or \code{"default"}, each set to \code{"linear"}, \code{"drop"}, \code{"rw2"} (non-linear RW2 with default PC prior \code{u = 5}, \code{alpha = 0.01}), or \code{list(model = "rw2", u = ..., alpha = ...)} for custom PC prior parameters.
#' @param coupling.intercept Character. Controls how the regional intercept inherits information from the global intercept. Options:
#'   \itemize{
#'     \item \code{"unpooled"} (default): Independent intercepts (no borrowing of strength).
#'     \item \code{"ordered_hierarchical"}: Regional intercept borrows strength from the global (Zhou & Bradley, 2024).
#'     \item \code{"bayesian_feedback"}: Sequential updating using global posterior as regional prior via moment matching (mean and precision). A warning is issued if the global posterior is strongly skewed (|skewness| > 1), as moment matching may be unreliable in that case. See Figueira et al. (2024).
#'     \item \code{NULL}: Regional-only model (no global component).
#'   }
#' @param coupling.covariates Character or list. Controls how global and regional covariate effects are related. Options:
#'   \itemize{
#'     \item \code{"unpooled"} (default): Independent estimation of global and regional coefficients. No shared hyperparameter; not a hierarchical mechanism.
#'     \item \code{"nested_shrinkage"}: Additive hierarchical shrinkage model: \code{beta_RE = beta_GL + delta_RE}, where \code{beta_GL} is shared via INLA's \code{copy} mechanism (\code{beta} fixed at 1, i.e. an exact, non-estimated copy), and \code{delta_RE ~ N(0, sigma_delta^2)} with \code{sigma_delta} estimated via a PC-prior. The amount of cross-scale borrowing is learned from the data: \code{sigma_delta} shrinks toward zero when regional and global slopes agree, and grows when regional evidence supports a distinct slope, including a change of sign. Uses unified (global-based) Z-standardization so the shared coefficient has consistent meaning across scales (no raster resampling involved, only a shared scaling constant). Requires the variable to be present at both scales, a linear global effect, and a joint (non-sequential) model fit.
#'     \item \code{"ordered_hierarchical"}: Multiplicative hierarchical coupling: \code{beta_RE = beta_copy * beta_GL}, via INLA's \code{copy} mechanism with \code{beta} freely estimated (Krainski et al., 2018; conceptual precedent in shared-component models, Knorr-Held & Best, 2001). Uses unified (global-based) Z-standardization, as above. Requires the variable to be present at both scales and a linear global effect.
#'     \item \code{"scale_decomposed"}: NOT a hierarchical coupling mechanism (no shared hyperparameter between global and regional coefficients). Shared covariates are decomposed into a large-scale trend (\code{X_glo_res = X_glo resampled}) and a small-scale anomaly (\code{X_reg_anom = X_reg - X_glo}), each with its own independent coefficient, following Cressie & Wikle (2011) and Zhou & Bradley (2024). Serves the multiscale objective via covariate orthogonalization rather than coefficient partial-pooling. In the output, coefficients are named \code{XGL_glo_res} and \code{XRE_reg_anom}.
#'     \item \code{"bayesian_feedback"}: Sequential updating using global posteriors as regional priors via moment matching (Figueira et al., 2024).
#'     \item \code{NULL}: Regional-only covariate effects; global covariates are excluded from the regional predictor entirely. No shared spatial fields (S_shared).
#'   }
#' @param proj.new.env Logical. Whether to compute predictions under new environmental scenarios (default: \code{TRUE}).
#' @param cv.folds Integer. Number of k-folds for spatial cross-validation (default: 1 = no CV).
#' @param n.threads Integer. Number of parallel threads to be used by INLA/inlabru (default: 1).
#' @param inla.int.strategy Character. INLA hyperparameter integration strategy. One of \code{"eb"} (Empirical Bayes, default), \code{"ccd"} (Central Composite Design), or \code{"grid"} (full grid integration). \code{"eb"} is fast but underestimates uncertainty by fixing hyperparameters at their posterior mode. \code{"ccd"} integrates over hyperparameters, at higher computational cost. See Rue et al. (2009) and Simpson et al. (2017).
#' @param seed Optional integer to set the random seed for reproducibility.
#' @param save.output Logical. If \code{TRUE}, saves key model outputs to disk.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages during model construction and fitting.
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
#' \item{formula}{List with the model's linear predictor components: \code{components} (the full inlabru component formula, intercepts + spatial fields + covariate terms), \code{rhs_global} (right-hand side of the global likelihood, or \code{NULL} if \code{coupling.intercept = NULL}), and \code{rhs_regional} (right-hand side of the regional likelihood).}
#' \item{Summary}{Named list of \code{data.frame}s with: \code{Metadata} (model configuration), \code{Model fit} (DIC, WAIC, MLPD), \code{Hyperparameters} (posterior range and sigma of spatial fields), \code{Intercepts} (IGlobal, IRegional, beta_copy), \code{Fixed effects} (covariate coefficients with CIs and significance), \code{Predictive performance} (AUC, Brier, BSS, RMSE), \code{Diagnostics} (Moran's I, SSI, r2_fields, range ratio, field correlation).}
#'
#' @details
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
#' This option implements a hierarchical Bayesian nested structure in which the regional intercept borrows strength from the global one. 
#' The regional intercept is expressed as a scaled copy of the global one using an INLA copy structure (IRegional = β * IGlobal), with β given an informative prior centred at 1 (Normal(1, sd = 0.5)), so that the regional level is encouraged (but not forced) to follow the global signal.
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
#'     If a covariate is not listed under `global`/`regional`, it inherits from `default`.
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
#' coupling.covariates: control of cross-scale covariate effects
#'
#' bayesian_feedback for new scenarios:
#' When using \code{coupling.covariates = "bayesian_feedback"} with \code{proj.new.env = TRUE},
#' note that the regional priors are fixed based on the present-day global posteriors.
#' They are NOT updated for future scenarios. This is a current limitation of any sequential
#' (two-stage) fitting protocol. Joint (single-fit) alternatives such as
#' \code{coupling.covariates = "nested_shrinkage"}, \code{"ordered_hierarchical"}, or
#' \code{"scale_decomposed"} do not have this limitation, since there is no separate
#' prior-injection step to go stale under new environmental scenarios.
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
#' Krainski, E.T., Gomez-Rubio, V., Bakka, H., Lenzi, A., Castro-Camilo, D., Simpson, D.,
#' Lindgren, F. & Rue, H. (2018). \emph{Advanced Spatial Modeling with Stochastic Partial
#' Differential Equations Using R and INLA}. CRC Press. (Chapter 4: the \code{copy} feature
#' for sharing/scaling effects across likelihoods.)
#'
#' Knorr-Held, L. & Best, N.G. (2001). A shared component model for detecting joint and
#' selective clustering of two diseases. \emph{Journal of the Royal Statistical Society:
#' Series A}, 164(1), 73-85.
#'
#' @export
MBM.Modelling <- function(jmbm_obj,
                      covariate.effects = NULL,
                      coupling.intercept = "unpooled",
                      coupling.covariates = "unpooled",
                      spde.mesh = NULL,
                      regional.pcprior.range = NULL,
                      regional.pcprior.sigma = NULL,
                      shared.pcprior.range = NULL,
                      shared.pcprior.sigma = NULL,
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


  .info("sabinaMBM: Multiscale Species Distribution Model")

  ## Checks
  if(!inherits(jmbm_obj, "nsdm.vinput")) {
    .stop("The 'jmbm_obj' must be of class 'nsdm.vinput'. Please see sabinaNSDM::NSDM.SelectCovariates().")
  }
  #
  fam <- "binomial"
  lnk <- "logit"
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
  valid_cp <- c("NULL", "unpooled", "ordered_hierarchical", "scale_decomposed", "bayesian_feedback", "nested_shrinkage")
  if(!is.null(coupling.covariates)) {
    if(is.character(coupling.covariates)) {
      if(!coupling.covariates %in% valid_cp) {
        .stop(paste0("`coupling.covariates` must be one of: ", paste(valid_cp, collapse = ", ")))
      }
    }
    else if(is.list(coupling.covariates)) {
      if(any(!names(coupling.covariates) %in% c("default", "variables"))) {
       .stop("Invalid entries in `coupling.covariates`. Allowed: 'default', 'variables'.")
      }
      if(!is.null(coupling.covariates$default)) {
        if(!coupling.covariates$default %in% valid_cp) {
          .stop(paste0("`coupling.covariates$default` must be one of: ", paste(valid_cp, collapse = ", ")))
        }
      }
      if(!is.null(coupling.covariates$variables)) {
        if (!is.list(coupling.covariates$variables)) .stop("`coupling.covariates$variables` must be a named list.")
        for(v in names(coupling.covariates$variables)) {
          mode_v <- coupling.covariates$variables[[v]]
          if(!mode_v %in% valid_cp) {
            .stop(paste0("Invalid coupling mode for variable '", v, "'. Allowed: ", paste(valid_cp, collapse = ", ")))
          }
        }
      }
    }
    else {
      .stop("`coupling.covariates` must be NULL, a character mode, or a list.")
    }
  }
  # compatibility coupling.covariates x covariate.effects
  has_bf <- (!is.null(coupling.intercept) && coupling.intercept == "bayesian_feedback")
  has_joint <- (!is.null(coupling.intercept) && coupling.intercept %in% c("ordered_hierarchical", "scale_decomposed")) ||
    (length(vr) > 0 && any(vapply(vr, function(v) .resolve_coupling_predictor(v, coupling.covariates, vg, vr) %in% c("ordered_hierarchical", "scale_decomposed", "nested_shrinkage"), logical(1))))
  vg_check <- if(is.null(coupling.intercept)) character(0) else vg
  if(length(vr) > 0) {
    for(v in vr) {
      cp_mode <- .resolve_coupling_predictor(v, coupling.covariates, vg_check, vr)
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
      } else if(cp_mode %in% c("ordered_hierarchical", "scale_decomposed", "nested_shrinkage")) {
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
        if(cp_mode %in% c("ordered_hierarchical", "nested_shrinkage") && !is.null(spec_gl) && spec_gl$model == "rw2") {
          .stop(paste0("Variable '", v, "' cannot use '", cp_mode, "'. Global effect uses RW2; ", cp_mode, " requires a linear global effect (beta_GL is copied/shared as-is)."))
        }
      }
    }
  }
  if(has_bf) {
    intercept_ok  <- is.null(coupling.intercept) || coupling.intercept == "bayesian_feedback"
    predictors_ok <- length(shared_vr) == 0 || all(vapply(shared_vr, function(v) {
      .resolve_coupling_predictor(v, coupling.covariates, vg_check, vr) %in% c("bayesian_feedback", "NULL")
    }, logical(1)))
    if(!intercept_ok || !predictors_ok) {
      .stop("'bayesian_feedback' must be used for the intercept and every shared covariate simultaneously, or not at all. Mixing it with 'unpooled', 'ordered_hierarchical', 'scale_decomposed', or 'nested_shrinkage' in the same call would leave those components uninformed by the global likelihood in the final fit.")
    }
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
    .stop(paste0("Requested `n.threads` (", n.threads, ") exceeds available cores (", available_cores, "). Your optimism is inspiring, your CPU is begging for mercy."))
  }
  INLA::inla.setOption(num.threads = n.threads)
  #
  if(is.null(coupling.intercept)) {
    .check("Architecture: Regional-only (no global component)")
  } else {
    cp_pred_str <- if(is.list(coupling.covariates)) "variable-specific" else if(is.null(coupling.covariates)) "none" else coupling.covariates
    .check(paste0("Architecture: Coupled model (Intercepts: ", coupling.intercept, " | Covariates: ", cp_pred_str, ")"))
  }
  fam_str <- paste0(fam, " (link: ", lnk, ")")
  .check(paste0("Family: ", fam_str))


  ## Data preparation
  sp_covglo <- terra::unwrap(jmbm_obj$IndVar.Global.Selected)
  sp_covreg <- terra::unwrap(jmbm_obj$IndVar.Regional.Selected)
  crs <- sf::st_crs(sp_covglo)

  all_model_vars <- unique(c(names(sp_covglo), names(sp_covreg)))
  scale_params_glo <- list()
  scale_params_reg <- list()

  # Identify shared variables whose coupling demands unified Z-standardization.
  unified_vars <- shared_vr[vapply(shared_vr, function(X) {
    .resolve_coupling_predictor(X, coupling.covariates, vg, vr) %in% c("scale_decomposed", "bayesian_feedback", "nested_shrinkage", "ordered_hierarchical")
  }, logical(1))]

  # scale_decomposed (unified_vars subset): needs extra steps for macro/anomaly Z-score decomposition
  sd_vars <- shared_vr[vapply(shared_vr, function(X) {
    .resolve_coupling_predictor(X, coupling.covariates, vg, vr) %in% c("scale_decomposed")
  }, logical(1))]
  if(length(all_model_vars) > 0) {
    .info("Pre-processing environmental covariates...")
    if(length(unified_vars) > 0) {
      .check(sprintf("Unified Z-standardization (mu_GL, sigma_GL) for shared coupled variables: %s",
                     paste(unified_vars, collapse = ", ")))
    }
    .check("Standardizing covariates (Z-score)")
  }

  # Single-pass standardization: unified for shared/coupled vars, native otherwise.
  result_std <- .standardize_multiscale(
    sp_covglo = if (!is.null(coupling.intercept)) sp_covglo else NULL,
    sp_covreg = sp_covreg,
    unified_vars     = unified_vars,
    standardize_global = !is.null(coupling.intercept))

  if (!is.null(result_std$sp_covglo)) sp_covglo <- result_std$sp_covglo
  sp_covreg <- result_std$sp_covreg
  scale_params_glo <- result_std$scale_params_glo
  scale_params_reg <- result_std$scale_params_reg

  sp_covglo_backup <- sp_covglo
  sp_covreg_backup <- sp_covreg

  # scale_decomposed: Post-standardization macro/anomaly decomposition. 
  # Z_anom = (X_RE - X_GL_res) / sigma_GL; both components now Z-standardized for regional predictor. 
  if(length(sd_vars) > 0) {
    .check("Computing scale decomposition (macro-trend vs micro-anomaly) on unified Z-scores")
    glo_resampled_stack <- terra::resample(sp_covglo[[sd_vars]], sp_covreg[[sd_vars[1]]])
    for(X in sd_vars) {
      sp_covreg[[paste0(X, "_glo_res")]]  <- glo_resampled_stack[[X]]
      sp_covreg[[paste0(X, "_reg_anom")]] <- sp_covreg[[X]] - glo_resampled_stack[[X]]
      for(suffix in c("_glo_res", "_reg_anom")) {
        v_name <- paste0(X, suffix)
        stat <- terra::global(sp_covreg[[v_name]], fun = c("mean", "sd"), na.rm = TRUE)
        scale_params_reg[[v_name]] <- list(mean = stat[1,1], sd = stat[1,2], unified = TRUE)
      }
    }
  }

  scale_params <- c(
    scale_params_glo,
    scale_params_reg[setdiff(names(scale_params_reg), names(scale_params_glo))]
  )

  # log stats
  if (length(names(sp_covreg)) > 0 && verbose) {
    for (v in names(sp_covreg)) {
      # skip layers generated by scale_decomposed; their stats are derived
      if (length(sd_vars) > 0 && (v %in% c(paste0(sd_vars, "_glo_res"),
                                            paste0(sd_vars, "_reg_anom")))) next
      m_val <- terra::global(sp_covreg[[v]], "mean", na.rm = TRUE)[1, 1]
      s_val <- terra::global(sp_covreg[[v]], "sd",   na.rm = TRUE)[1, 1]
      original_mean <- scale_params_reg[[v]][["mean"]]
      original_sd   <- scale_params_reg[[v]][["sd"]]
      unified_flag  <- isTRUE(scale_params_reg[[v]][["unified"]])
      tag <- if (unified_flag) "regional, unified (sigma_GL)" else "regional"
      .item(sprintf("%s (%s): Z-mean = %+.3f, Z-sd = %.3f | Original: mean = %+.2f, sd = %.2f",
                    v, tag, m_val, s_val, original_mean, original_sd))
    }
  }

  if (!is.null(coupling.intercept) && !is.null(sp_covglo) && length(names(sp_covglo)) > 0 && verbose) {
    for (v in names(sp_covglo)) {
      m_val <- terra::global(sp_covglo[[v]], "mean", na.rm = TRUE)[1, 1]
      s_val <- terra::global(sp_covglo[[v]], "sd",   na.rm = TRUE)[1, 1]
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
    cbind(jmbm_obj$SpeciesData.XY.Global, resp = if(!is.null(jmbm_obj$Response.Global)) jmbm_obj$Response.Global else 1L),
    cbind(bg_or_abs_glo, resp = 0L)     # si lo hacemo asi poner algun check con stop/warning para que datos y family sean coherentes
  )
  pp_reg <- rbind(
    cbind(jmbm_obj$SpeciesData.XY.Regional, resp = if(!is.null(jmbm_obj$Response.Regional)) jmbm_obj$Response.Regional else 1L),  
    cbind(bg_or_abs_reg, resp = 0L)          
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

  # background vs real absences flag, used by .bf_intercept_offset
  bf_uses_background <- is.null(jmbm_obj$Absences.XY.Global) && is.null(jmbm_obj$Absences.XY.Regional)

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
    if (has_Sre) .check("S_re (fine-scale field specific to the regional model)")
    if (has_Sshared) .check("S_shared (broad-scale field shared across both scales)")
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


  ## Prior constants
  # prior for intercepts iid. param = c(1, 0.01) -> P(σ > 1) = 0.01
  IID_PRIOR <- "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"
  COPY_BETA_PRIOR <- "hyper = list(beta = list(prior = 'normal', param = c(1, 4)))"
  SLOPE_DELTA_PRIOR <- "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"

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
                  coupling.covariates = coupling.covariates,
                  slope_delta_prior = SLOPE_DELTA_PRIOR,
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

  # if coupling.intercept = NULL, desactivate bayesian feedback
  needs_feedback <- !is.null(coupling.intercept) && (
    coupling.intercept == "bayesian_feedback" ||
    any(vapply(vr, function(v) .resolve_coupling_predictor(v, coupling.covariates, vg, vr) == "bayesian_feedback", logical(1)))
  )

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
                              coupling.intercept)
  lik_glo <- liks$lik_glo
  lik_reg <- liks$lik_reg

  eta_terms <- c(
    "IRegional",
    if(has_Sshared) "Sshared" else NULL,
    if(has_Sre) "Sre" else NULL,
    cmp_cov$like$fregional
  )
  eta <- paste(stats::na.omit(eta_terms), collapse = " + ")

  pred_formula <- switch(
    lnk,
    "logit" = as.formula(paste0("~ 1 / (1 + exp(-(", eta, ")))")),
    "cloglog" = as.formula(paste0("~ 1 - exp(-exp(", eta, "))")),
    "log" = as.formula(paste0("~ exp(", eta, ")")),
    "identity" = as.formula(paste0("~ ", eta)))


  ## Assemble likelihoods and fit model
  lik_list <- list()
  if(!is.null(coupling.intercept)) lik_list <- c(lik_list, list(lik_glo))
  lik_list <- c(lik_list, list(lik_reg))

  .info(sprintf("Fitting Bayesian model (Integration strategy: '%s')...", inla.int.strategy))

  fit <- .fit_jmbm(
    cmp = cmp,
    lik_list = lik_list,
    coupling.intercept = coupling.intercept,
    coupling.covariates = coupling.covariates,
    needs_feedback = needs_feedback,
    vr = vr,
    n.threads = n.threads,
    seed = seed,
    int.strategy = inla.int.strategy,
    bf_delta_int = .bf_intercept_offset(pp_glo$resp, pp_reg$resp, bf_uses_background),
    has_Sshared = has_Sshared,
    spde.mesh = spde.mesh,
    verbose = verbose
  )


  ## K-fold CV
  if(cv.folds > 1) {
    .info(sprintf("Performing spatial cross-validation (K = %d folds)...", cv.folds))
    # stratified k folds
    folds_g <- if(!is.null(coupling.intercept)) .make_stratified_kfolds(pp_glo$resp, cv.folds) else NULL
    folds_r <- .make_stratified_kfolds(pp_reg$resp, cv.folds)

    n_cores_cv <- min(n.threads, cv.folds)
    inla_threads_k <- max(1, floor(n.threads / n_cores_cv))

    if(n_cores_cv > 1) {
      future::plan(future::multisession, workers = n_cores_cv, quiet = TRUE)
    }

    both_cov <- c(sp_covglo, sp_covreg)
    coords_all_r <- sf::st_coordinates(pp_reg)
    valid_pts_mask <- stats::complete.cases(terra::extract(both_cov, coords_all_r, ID = FALSE))

    cv_worker <- function(k) {
      train_g <- if(!is.null(coupling.intercept)) pp_glo[folds_g != k, ] else NULL
      train_r <- pp_reg[folds_r != k, ]
      test_r  <- pp_reg[folds_r == k, ]

      liks_k <- .build_likelihoods(fam, lnk, rhs_glo, rhs_reg,
                                    train_g, train_r,
                                    bdy_glo, bdy_reg, dom,
                                    coupling.intercept)
      lik_list_k <- Filter(Negate(is.null), list(liks_k$lik_glo, liks_k$lik_reg))

      # fit fold k
      fit_k <- .fit_jmbm(
        cmp = cmp,
        lik_list = lik_list_k,
        coupling.intercept = coupling.intercept,
        coupling.covariates = coupling.covariates,
        needs_feedback = needs_feedback,
        vr = vr,
        n.threads = inla_threads_k,
        seed = seed,
        int.strategy = inla.int.strategy,
        bf_delta_int = .bf_intercept_offset(train_g$resp, train_r$resp, bf_uses_background),
        has_Sshared = has_Sshared,
        spde.mesh = spde.mesh,
        verbose = verbose
      )

      # rm NAs
      keep_k <- valid_pts_mask[folds_r == k]
      test_r2 <- test_r[keep_k, ]

      if(nrow(test_r2) > 0) {
        test_r2 <- sf::st_transform(test_r2, terra::crs(sp_covreg))
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
      cv_metrics_list <- future.apply::future_lapply(seq_len(cv.folds), cv_worker, future.seed = TRUE)
      future::plan(future::sequential) # Reset
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


  #### New scenarios
  proj_list <- list()
  if(proj.new.env && !is.null(jmbm_obj$Scenarios)) {
    .check(sprintf("Projecting to new scenarios: %s", paste(names(jmbm_obj$Scenarios), collapse = ", ")))
    
    for(sc in names(jmbm_obj$Scenarios)) {
      scen_rast <- terra::unwrap(jmbm_obj$Scenarios[[sc]])
      
      if(!identical(terra::crs(scen_rast), terra::crs(sp_covreg_backup))) {
        scen_rast <- terra::project(scen_rast, terra::crs(sp_covreg_backup))
      }
      
      scen_reg_full <- terra::resample(scen_rast, sp_covreg_backup[[1]])
      pred_template <- scen_reg_full[[rep(1, terra::nlyr(sp_covreg_backup))]]
      names(pred_template) <- names(sp_covreg_backup)

      vars_glo <- intersect(names(sp_covglo_backup), names(scen_rast))
      vars_reg <- intersect(names(sp_covreg_backup), names(scen_rast))
   
      # resampling
      if(length(vars_glo) > 0) {
        fact <- max(1, round(terra::res(sp_covglo_backup)[1] / terra::res(scen_rast)[1]))
        scen_glo_raw <- if(fact >= 2) {
          terra::resample(terra::aggregate(scen_rast[[vars_glo]], fact = fact, fun = "mean", na.rm = TRUE),
                           sp_covglo_backup)
        } else {
          terra::resample(scen_rast[[vars_glo]], sp_covglo_backup)
        }
      } else {
        scen_glo_raw <- NULL
      }

      scen_reg_raw <- if(length(vars_reg) > 0) terra::resample(scen_rast[[vars_reg]], sp_covreg_backup) else NULL

      
      # standarization uni (vars simples)
      for(v in setdiff(vars_reg, sd_vars)) {
        if (!is.null(scale_params_reg[[v]])) {
          pred_template[[v]] <- (scen_reg_raw[[v]] - scale_params_reg[[v]]$mean) / scale_params_reg[[v]]$sd
        }
      }
      
      # scale-decomposed standarization uni
      if(length(sd_vars) > 0) {
        vars_sd <- intersect(sd_vars, names(scen_rast))
        if(length(vars_sd) > 0) {
          raw_glo_base <- scen_glo_raw[[vars_sd]]
          glo_resampled_scenario <- terra::resample(raw_glo_base, sp_covreg_backup)
          
          for(X in vars_sd) {
            raw_reg <- scen_reg_raw[[X]]
            raw_glo <- glo_resampled_scenario[[X]]
            
            p_unif <- scale_params_glo[[X]]
            pred_template[[paste0(X, "_glo_res")]]  <- (raw_glo - as.numeric(p_unif$mean)) / as.numeric(p_unif$sd)
            pred_template[[paste0(X, "_reg_anom")]] <- (raw_reg - raw_glo) / as.numeric(p_unif$sd)
          }
        }
      }
      
      sp_covreg <- pred_template 
      
      if(!is.null(scen_glo_raw)) {
        pred_template_glo <- terra::rast(sp_covglo_backup)
        for(v in vars_glo) {
          if(!is.null(scale_params_glo[[v]])) {
            pred_template_glo[[v]] <- (scen_glo_raw[[v]] - scale_params_glo[[v]]$mean) / scale_params_glo[[v]]$sd
          }
        }
        sp_covglo <- pred_template_glo
      }
      
      # prediction fut
      scen_pts <- terra::as.points(pred_template, values = TRUE, na.rm = TRUE)
      scen_df <- sf::st_as_sf(scen_pts)
      sf::st_crs(scen_df) <- crs
      scen_df$region <- 1L
      
      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[sc]] <- .pred_as_tif(proj_pred, template = pred_template)
      
      rm(scen_pts, scen_df, proj_pred, pred_template, scen_glo_raw, scen_reg_raw, scen_reg_full)
      if(exists("pred_template_glo")) rm(pred_template_glo)
      gc(verbose = FALSE)
    }
    
    sp_covglo <- sp_covglo_backup
    sp_covreg <- sp_covreg_backup
  }


  ## Diagnostics
  .info("Computing model diagnostics...")

  coords_reg <- sf::st_coordinates(pp_reg)
  data_used  <- data.frame(x = coords_reg[,1], y = coords_reg[,2], resp = pp_reg$resp)
  diag_block <- .jmbm_diagnostics(
                  fit,
                  fam,  
                  data_used,
                  n_glo = if(needs_feedback) 0L else if(!is.null(coupling.intercept)) nrow(pp_glo) else 0L,
                  # needs_feedback: fit is regional-only, no global block to skip
                  priors = list(regional.pcprior.range = regional.pcprior.range,
                                regional.pcprior.sigma = regional.pcprior.sigma,
                                shared.pcprior.range = shared.pcprior.range,
                                shared.pcprior.sigma = shared.pcprior.sigma),
                  pred_Sre = if(has_Sre) pred_Sre else NULL,
                  pred_Sshared = if(has_Sshared) pred_Sshared else NULL,
                  coupling.intercept = coupling.intercept,
                  scale_params_glo = scale_params_glo,
                  scale_params_reg = scale_params_reg)


  species <- jmbm_obj$Species.Name

  # summary
  summary_df <- .jmbm_generate_summary(
    fit = fit, 
    species = species, 
    fam = fam, 
    lnk = lnk, 
    coupling.intercept = coupling.intercept, 
    coupling.covariates = coupling.covariates, 
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
      coupling.covariates = coupling.covariates,
      covariate.effects = covariate.effects,
      proj.new.env = proj.new.env,
      cv.folds = cv.folds,
      n.threads = n.threads,
      inla.int.strategy = inla.int.strategy,
      seed = seed
    ),
    Selected.Variables.Global = jmbm_obj$Selected.Variables.Global,
    Selected.Variables.Regional = jmbm_obj$Selected.Variables.Regional,
    formula = list(
      components   = cmp,
      rhs_global   = rhs_glo,
      rhs_regional = rhs_reg
    ),
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
    scale_params_glo = scale_params_glo,
    scale_params_reg = scale_params_reg,
    pit_values = diag_block$calibration$pit_values,
    diagnostic_data = diag_block$diagnostic_data,
    Summary = summary_df
  )
 
  attr(sabina, "class") <- "jmbm.inlabru"


  ## Save outputs
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
    # CPO values (one per observation)
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
      list(key = "hyperparams",   file = "_hyperparams.png",       w = 6, h = 5),
      list(key = "correlogram",   file = "_correlogram.png",       w = 6, h = 5),
      list(key = "hist",          file = "_residualHistogram.png", w = 6, h = 5),
      list(key = "qq",            file = "_QQplot.png",            w = 6, h = 4),
      list(key = "Srefields",     file = "_SPDEfields.png",        w = 6, h = 5),
      list(key = "semivariogram", file = "_semivariogram.png",     w = 6, h = 4)
    )
    for(ps in plot_specs) {
      p <- tryCatch(plot(sabina, which = ps$key), error = function(e) NULL)
      if(!is.null(p)) {
        ggplot2::ggsave(
          filename = file.path(values_path, paste0(species, ps$file)),
          plot = p, width = ps$w, height = ps$h, dpi = 300, bg = "white")
      }
    }
    # predictive metrics
    pred_metrics <- data.frame(
      Metric = c("AUC", "Tjur_R2", "Brier", "BSS", "RMSE", "Cor_obs_pred",
                 "Calibration_slope", "Coverage_95", "MLPD", "PIT_KS_pvalue"),
      Value = c(diag_block$predictive$auc_full,
                diag_block$predictive$tjur_r2,
                diag_block$predictive$brier,
                diag_block$predictive$bss,
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


  if(verbose) message("\n  ✓ Model fitted successfully. Check 'summary()' for evaluation metrics.\n")

  return(sabina)

}


###----------------###




