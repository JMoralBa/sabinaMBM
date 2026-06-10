#' @name MBM.SuggestParams
#'
#' @title Suggest parameters for \code{create_mesh()} and \code{MBM.Modelling()}
#'
#' @description
#' Derives theory-based suggestions for the key numerical parameters required by
#' \code{create_mesh()} and \code{MBM.Modelling()} in the sabinaMBM package.
#' Specifically, it suggests:
#' \enumerate{
#'   \item \strong{Mesh geometry} (\code{edge}, \code{offset}): triangle size and
#'         boundary buffer for \code{create_mesh()}, derived from the expected
#'         spatial autocorrelation range of the fine-scale field (S_RE).
#'   \item \strong{SPDE PC priors} (\code{regional.pcprior.range},
#'         \code{regional.pcprior.sigma}, \code{shared.pcprior.range},
#'         \code{shared.pcprior.sigma}): penalized complexity priors for the
#'         fine-scale residual field (S_RE) and the broad-scale shared field
#'         (S_shared), derived sequentially from a single master spatial parameter
#'         following Fuglstad et al. (2019) and Bakka et al. (2018).
#'   \item \strong{Covariate functional form} (\code{covariate.effects}): whether
#'         each covariate at each scale is better represented as a linear effect or
#'         a non-linear rw2 spline, assessed via an auxiliary GAM test (Wood, 2017).
#' }
#'
#' All outputs are suggestions based on theoretical guidelines and exploratory
#' data analysis. The user retains full control over the final parameter values
#' passed to \code{create_mesh()} and \code{MBM.Modelling()}.
#'
#' @param jmbm_obj An object of class \code{nsdm.vinput} produced by
#'   \code{sabinaNSDM::NSDM.SelectCovariates()}.
#' @param expected.range.fraction Numeric scalar in (0, 1]. Assumed spatial
#'   autocorrelation range of the fine-scale process (S_RE) expressed as a
#'   fraction of the regional domain diameter. This is the master parameter from
#'   which all spatial suggestions are derived. Default: \code{0.20}, the INLA
#'   standard assumption that the effective range spans approximately one fifth of
#'   the study domain (Barber et al., 2016). Users should increase this value for
#'   species with broad continuous distributions (e.g., 0.40) and decrease it for
#'   highly fragmented or localised species (e.g., 0.05).
#' @param shared.range.ratio Numeric scalar >= 3. Assumed multiplicative factor by
#'   which the S_shared spatial process operates at a larger scale than S_RE.
#'   Controls scale separation between the two SPDE fields. Must be >= 3 to
#'   guarantee statistical identifiability (Bakka et al., 2018). Default:
#'   \code{3}, representing the minimum mathematical threshold required to prevent
#'   double-smoothing while separating local from broad-scale patterns.
#' @param shared.sigma.divisor Numeric scalar > 1. Divisor applied to
#'   \code{sigma_re} (1.5) to obtain \code{sigma_shared}. The shared field enters
#'   both likelihoods simultaneously and must be more tightly constrained to
#'   prevent variance cannibalization. Default: \code{5}, yielding
#'   \code{sigma_shared = 0.3}. On the logit scale, this restricts the shared 
#'   field's fluctuations to a narrow probability band, ensuring it acts as a 
#'   background modulator (Simpson et al., 2017).
#' @param range.alpha Numeric scalar in (0, 1). Probability for PC priors on
#'   spatial range: P(range < u) = \code{range.alpha}. Default: \code{0.05}.
#'   Following Fuglstad et al. (2019).
#' @param sigma.alpha Numeric scalar in (0, 1). Probability for PC priors on
#'   marginal standard deviation: P(sigma > u) = \code{sigma.alpha}.
#'   Default: \code{0.01}. Following Simpson et al. (2017).
#' @param rw2.alpha Numeric scalar in (0, 1). Significance threshold for the GAM
#'   non-linearity test (p-value of smooth term edf). Default: \code{0.01}
#'   (conservative, consistent with the PC prior philosophy of penalising
#'   complexity unless data provide strong evidence).
#' @param rw2.edf.min Numeric scalar > 1. Minimum effective degrees of freedom
#'   (edf) of the GAM smooth term required to suggest rw2, in addition to the
#'   significance criterion. Filters spurious significance due to large N.
#'   Default: \code{1.5}.
#' @param rw2.u Numeric scalar. Value of \code{u} in the PC prior for rw2
#'   covariate effects: P(sigma_spline > u) = \code{rw2.alpha.pc}. Default:
#'   \code{0.5}. Combined with a low alpha, this conservative prior strongly 
#'   penalizes structural deviations from linearity to prevent overfitting 
#'   (Simpson et al., 2017).
#' @param rw2.alpha.pc Numeric scalar. Value of \code{alpha} in the PC prior for
#'   rw2 covariate effects. Default: \code{0.01}.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted summary
#'   of all suggestions with their rationale to the console.
#'
#' @return A standard named list containing four main elements, structured so that 
#'   its contents can be passed directly to the package's core functions:
#'   \describe{
#'     \item{\code{create_mesh_args}}{Named list with suggested \code{edge} (length-2
#'       numeric vector) and \code{offset} (length-2 numeric vector), ready for \code{create_mesh()}.}
#'     \item{\code{MBM.Modelling_args}}{Named list containing all the suggested arguments 
#'       for \code{MBM.Modelling()}: \code{regional.pcprior.range}, \code{regional.pcprior.sigma}, 
#'       \code{shared.pcprior.range}, \code{shared.pcprior.sigma} (NULL for regional-only models), 
#'       and \code{covariate.effects} (which includes the specific non-linear/linear suggestions).}
#'     \item{\code{covariate_diagnostics}}{Named list with \code{regional_tests} and 
#'       \code{global_tests}. These contain data.frames with the GAM non-linearity test results 
#'       (effective degrees of freedom, p-values, and rationale) for each covariate.}
#'     \item{\code{spatial_metadata}}{Named list containing spatial context and intermediate 
#'       quantities: \code{crs_units} (degrees or metres), \code{regional_diameter}, 
#'       \code{range_re_expected} (the master ecological parameter), and \code{scale_separation_ratio}.}
#'   }
#'
#' @details
#' \strong{CRS detection and units.} The function detects whether the input
#' rasters use geographic coordinates (degrees, e.g. WGS84) or a projected CRS
#' (metres, e.g. UTM). All spatial suggestions are in native CRS units.
#'
#' \strong{Sequential derivation from a single master parameter.}
#' All spatial suggestions derive from a single ecologically meaningful quantity:
#' the expected autocorrelation range of the fine-scale process
#' (\code{range_re_expected = expected.range.fraction × regional_diameter}).
#' The default fraction of 0.20 is the INLA standard assumption (Barber et al.,
#' 2016) and should be adjusted based on species ecology. The full derivation
#' cascade is:
#' \enumerate{
#'   \item \strong{Master parameter:}
#'     \code{range_re_expected = expected.range.fraction × diam_re}
#'   \item \strong{S_RE range prior lower bound:}
#'     \code{rho_0_RE = range_re_expected / 10}
#'     (Fuglstad et al. 2019: optimal prior calibration sets the lower bound at
#'     1/10 of the expected range)
#'   \item \strong{S_shared range prior lower bound:}
#'     \code{rho_0_shared = (shared.range.ratio × range_re_expected) / 10}
#'     (same Fuglstad 2019 calibration applied to the shared field, whose
#'     expected range is \code{shared.range.ratio} times larger; Bakka et al. 2018)
#'   \item \strong{Mesh inner edge:}
#'     \code{edge[1] = range_re_expected / 5}
#'     (inner triangles must be smaller than range/5 to capture local spatial
#'     structure; SM_S1 sabinaMBM; Foster et al. 2024; Krainski et al. 2019)
#'   \item \strong{Mesh outer edge:}
#'     \code{edge[2] = 3 × edge[1]}
#'     (outer domain coarser to reduce computation; heuristic consistent with
#'     Anderson et al. 2022 and Lindgren et al. 2011)
#'   \item \strong{Inner offset:}
#'     \code{offset[1] = edge[1]}
#'     (minimum 1-triangle buffer to avoid mesh collapse at data boundary;
#'     no formal reference, practical heuristic)
#'   \item \strong{Outer offset:}
#'     \code{offset[2] = range_re_expected}
#'     (boundary buffer >= spatial range ensures Neumann boundary variance
#'     dissipates before reaching the study area; Krainski et al. 2019;
#'     Anderson et al. 2022)
#'   \item \strong{S_RE sigma prior upper bound:}
#'     \code{sigma_re = 1.5} (fixed)
#'     (on the logit scale, SD = 1.5 spans approximately the [0.05, 0.95]
#'     probability range, allowing S_RE to act as a residual spatial correction
#'     without overriding environmental predictors; Simpson et al. 2017, p.17)
#'   \item \strong{S_shared sigma prior upper bound:}
#'     \code{sigma_shared = sigma_re / shared.sigma.divisor}
#'     (S_shared enters both likelihoods and must be more tightly constrained to
#'     prevent variance cannibalization; divisor default of 5 gives sigma_shared
#'     = 0.3, consistent with case study S2; Simpson et al. 2017)
#' }
#'
#' \strong{Non-linearity test for rw2.}
#' For each covariate, a binomial GAM with thin-plate spline
#' (\code{mgcv::gam(resp ~ s(v, k = 5), family = binomial, method = "REML")})
#' is fitted to the presence/background data. rw2 is suggested only when
#' \emph{both}: (a) the smooth term p-value < \code{rw2.alpha}, AND (b) edf >
#' \code{rw2.edf.min}. The dual criterion prevents: (1) large-N spurious
#' significance (filtered by edf) and (2) overfitting to ecologically negligible
#' curvature (filtered by p). With k = 5 knots, only robust non-linearities are
#' detected, consistent with the PC prior philosophy (Wood, 2017).
#'
#' @references
#' Anderson, S.C. et al. (2022). sdmTMB: an R package for fast, flexible, and
#' user-friendly generalized linear mixed effects models with spatial and
#' spatiotemporal random fields.
#' \emph{bioRxiv}. \doi{10.1101/2022.03.24.485545}
#'
#' Bakka, H. et al. (2018). Spatial modelling with R-INLA: A review.
#' \emph{WIREs Computational Statistics}, 10, e1443.
#'
#' Barber, X. et al. (2016). Incorporating spatial structure into inclusion
#' probabilities for Bayesian variable selection in generalized linear models
#' with the spike-and-slab prior.
#' \emph{Journal of Applied Statistics}, 43, 2261–2279.
#'
#' Foster, S.D. et al. (2024). RISDM: species distribution modelling from
#' multiple data sources in R.
#' \emph{Ecography}, e06964.
#'
#' Fuglstad, G.A. et al. (2019). Constructing priors that penalize the complexity
#' of Gaussian random fields.
#' \emph{Journal of the American Statistical Association}, 114, 445–452.
#'
#' Krainski, E. et al. (2019). \emph{Advanced Spatial Modeling with Stochastic
#' Partial Differential Equations Using R and INLA}. CRC Press.
#'
#' Lindgren, F. & Rue, H. (2011). An explicit link between Gaussian fields and
#' Gaussian Markov random fields.
#' \emph{Journal of the Royal Statistical Society B}, 73, 423–498.
#'
#' Simpson, D. et al. (2017). Penalising model component complexity: a
#' principled, practical approach to constructing priors.
#' \emph{Statistical Science}, 32, 1–28.
#'
#' Wood, S.N. (2017). \emph{Generalized Additive Models: An Introduction with R}.
#' 2nd ed. CRC Press.
#'
#' @seealso \code{\link{create_mesh}}, \code{\link{MBM.Modelling}}
#'
#' @examples
#' \dontrun{
#' # After running the sabinaNSDM pipeline:
#' mySelvars <- sabinaNSDM::NSDM.SelectCovariates(myFormatting, ...)
#'
#' # Get parameter suggestions (defaults suitable for most regional SDMs)
#' params <- MBM.SuggestParams(mySelvars)
#'
#' # For a highly localised species, reduce expected.range.fraction:
#' params <- MBM.SuggestParams(mySelvars, expected.range.fraction = 0.05)
#'
#' # Use suggestions directly
#' myMesh <- create_mesh(
#'   nsdm_obj         = mySelvars,
#'   edge             = params$create_mesh$edge,
#'   offset           = params$create_mesh$offset,
#'   boundary.method  = "raster_mask"
#' )
#'
#' myMod <- MBM.Modelling(
#'   mbm_obj                = mySelvars,
#'   spde.mesh              = myMesh,
#'   regional.pcprior.range = params$regional.pcprior.range,
#'   regional.pcprior.sigma = params$regional.pcprior.sigma,
#'   shared.pcprior.range   = params$shared.pcprior.range,
#'   shared.pcprior.sigma   = params$shared.pcprior.sigma,
#'   covariate.effects      = params$covariate.effects
#' )
#' }
#'
#' @export
MBM.SuggestParams <- function(mbm_obj,
                               expected.range.fraction = 0.20,
                               shared.range.ratio      = 3,
                               shared.sigma.divisor    = 5,
                               range.alpha             = 0.05,
                               sigma.alpha             = 0.01,
                               rw2.alpha               = 0.01,
                               rw2.edf.min             = 1.5,
                               rw2.u                   = 0.5,
                               rw2.alpha.pc            = 0.01,
                               verbose                 = TRUE) {

  # ── 0. Input checks ────────────────────────────────────────────────────────

  if (!inherits(mbm_obj, "nsdm.vinput"))
    stop("'mbm_obj' must be of class 'nsdm.vinput'.\n",
         "  Please use sabinaNSDM::NSDM.SelectCovariates() to obtain it.")
  if (!is.numeric(expected.range.fraction) ||
      expected.range.fraction <= 0 || expected.range.fraction > 1)
    stop("'expected.range.fraction' must be in (0, 1]. ",
         "Default 0.20 is the INLA standard assumption (Barber et al., 2016).")
  if (!is.numeric(shared.range.ratio) || shared.range.ratio <= 3)
    stop("'shared.range.ratio' must be > 3 to guarantee scale separation ",
         "(Bakka et al., 2018).")
  if (!is.numeric(shared.sigma.divisor) || shared.sigma.divisor <= 1)
    stop("'shared.sigma.divisor' must be > 1: shared sigma must be more ",
         "constrained than regional sigma.")
  if (!is.numeric(range.alpha) || range.alpha <= 0 || range.alpha >= 1)
    stop("'range.alpha' must be in (0, 1).")
  if (!is.numeric(sigma.alpha) || sigma.alpha <= 0 || sigma.alpha >= 1)
    stop("'sigma.alpha' must be in (0, 1).")
  if (!is.numeric(rw2.alpha) || rw2.alpha <= 0 || rw2.alpha >= 1)
    stop("'rw2.alpha' must be in (0, 1).")
  if (!is.numeric(rw2.edf.min) || rw2.edf.min <= 1)
    stop("'rw2.edf.min' must be > 1 (a linear effect has edf = 1).")
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required for the non-linearity test. ",
         "Please install it.")

  # ── 1. Detect model type: joint vs regional-only ──────────────────────────

  has_global <- !is.null(mbm_obj$Selected.Variables.Global) &&
                length(mbm_obj$Selected.Variables.Global) > 0

  vars_gl <- if (has_global) mbm_obj$Selected.Variables.Global else character(0)
  vars_re <- mbm_obj$Selected.Variables.Regional

  # ── 2. Extract rasters and detect CRS ─────────────────────────────────────

  rast_re  <- terra::unwrap(mbm_obj$IndVar.Regional.Selected)
  rast_gl  <- if (has_global) terra::unwrap(mbm_obj$IndVar.Global.Selected) else NULL

  is_lonlat <- isTRUE(sf::st_is_longlat(sf::st_crs(rast_re)))
  crs_units <- if (is_lonlat) "degrees" else "metres (projected CRS)"

  # ── 3. Spatial resolutions (kept for diagnostics only) ────────────────────

  res_re <- mean(terra::res(rast_re[[1]]))
  res_gl <- if (has_global) mean(terra::res(rast_gl[[1]])) else NA_real_

  # ── 4. Domain diameters ────────────────────────────────────────────────────

  .diameter <- function(coords_df) {
    bb <- sf::st_bbox(
      sf::st_as_sf(as.data.frame(coords_df), coords = c("x", "y"))
    )
    sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  }

  coords_re <- rbind(
    mbm_obj$SpeciesData.XY.Regional[, c("x", "y")],
    mbm_obj$Background.XY.Regional[, c("x", "y")]
  )
  diam_re <- .diameter(coords_re)

  diam_gl <- if (has_global) {
    coords_gl <- rbind(
      mbm_obj$SpeciesData.XY.Global[, c("x", "y")],
      mbm_obj$Background.XY.Global[, c("x", "y")]
    )
    .diameter(coords_gl)
  } else NA_real_

  # ── 5 & 6. Sequential spatial parameter derivation ────────────────────────
  #
  # Master parameter: expected autocorrelation range of S_RE.
  # All spatial suggestions derive from this single quantity.
  #
  # Cascade:
  #   range_re_expected  = expected.range.fraction * diam_re
  #                        (INLA default: 20% of domain; Barber et al. 2016)
  #   rho_0_RE           = range_re_expected / 10
  #                        (Fuglstad et al. 2019: lower bound = 1/10 of expected range)
  #   rho_0_shared       = (shared.range.ratio * range_re_expected) / 10
  #                        (same Fuglstad calibration for shared field; Bakka et al. 2018)
  #   edge[1]            = range_re_expected / 5
  #                        (SM_S1 sabinaMBM: edge < range/5; Foster et al. 2024)
  #   edge[2]            = 3 * edge[1]
  #                        (heuristic; Anderson et al. 2022; Lindgren et al. 2011)
  #   offset[1]          = edge[1]
  #                        (practical heuristic: 1-triangle buffer, no formal ref)
  #   offset[2]          = range_re_expected
  #                        (>= range ensures Neumann variance dissipates;
  #                         Krainski et al. 2019; Anderson et al. 2022)
  #   sigma_re           = 1.5 (fixed)
  #                        (logit scale: SD=1.5 spans ~[0.05,0.95]; Simpson et al. 2017)
  #   sigma_shared       = sigma_re / shared.sigma.divisor  [default: 1.5/5 = 0.3]
  #                        (shared field more constrained; Simpson et al. 2017;
  #                         case study S2 Morales-Barbero et al.)

  range_re_expected <- expected.range.fraction * diam_re

  # PC prior lower bounds (Fuglstad et al. 2019)
  range_re     <- range_re_expected / 10
  range_shared <- if (has_global)
    (shared.range.ratio * range_re_expected) / 10 else NULL

  # Sigma prior upper bounds
  sigma_re     <- 1.5
  sigma_shared <- if (has_global) round(sigma_re / shared.sigma.divisor, 2) else NULL

  # Mesh geometry
  edge_inner   <- range_re_expected / 5
  edge_outer   <- 3 * edge_inner
  offset_inner <- edge_inner
  offset_outer <- range_re_expected

  # ── 7. Non-linearity test (GAM auxiliary) ─────────────────────────────────
  # Binomial GAM with thin-plate spline (k=5, conservative).
  # rw2 suggested when BOTH: p < rw2.alpha AND edf > rw2.edf.min.
  # Dual criterion avoids: (a) large-N spurious significance; (b) negligible curvature.

  .build_data <- function(xy_pres, xy_bg, rast) {
    coords <- rbind(
      cbind(xy_pres[, c("x", "y")], resp = 1L),
      cbind(xy_bg[,   c("x", "y")], resp = 0L)
    )
    vals <- terra::extract(rast, coords[, c("x", "y")])[, -1, drop = FALSE]
    cbind(data.frame(resp = coords$resp), vals)
  }

  .test_rw2 <- function(data, var_name, rw2.alpha, rw2.edf.min) {
    keep <- stats::complete.cases(data[, c("resp", var_name)])
    d    <- data[keep, , drop = FALSE]
    if (nrow(d) < 20)
      return(list(suggest = "linear", p = NA, edf = NA,
                  reason = "insufficient data (< 20 complete cases)"))
    v <- d[[var_name]]
    if (diff(range(v, na.rm = TRUE)) == 0)
      return(list(suggest = "linear", p = NA, edf = NA, reason = "no variation"))
    df_gam <- data.frame(resp = d$resp, v = v)
    fit <- tryCatch(
      mgcv::gam(resp ~ s(v, k = 5), data = df_gam,
                family = stats::binomial(link = "logit"),
                method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit))
      return(list(suggest = "linear", p = NA, edf = NA,
                  reason = "GAM failed to converge"))
    sm   <- summary(fit)
    edf  <- sm$s.table[1, "edf"]
    pval <- sm$s.table[1, "p-value"]
    suggest <- if (!is.na(pval) && pval < rw2.alpha && edf > rw2.edf.min)
      "rw2" else "linear"
    reason <- if (suggest == "rw2") {
      sprintf("p=%.4f < %.2f and edf=%.2f > %.1f", pval, rw2.alpha, edf, rw2.edf.min)
    } else if (!is.na(pval) && pval < rw2.alpha) {
      sprintf("p=%.4f significant but edf=%.2f <= %.1f (curvature too weak)",
              pval, edf, rw2.edf.min)
    } else {
      sprintf("p=%.4f >= %.2f (no evidence of non-linearity)", pval, rw2.alpha)
    }
    list(suggest = suggest, p = pval, edf = edf, reason = reason)
  }

  .rw2_entry <- function(suggest, rw2.u, rw2.alpha.pc) {
    if (suggest == "rw2") list(model = "rw2", u = rw2.u, alpha = rw2.alpha.pc)
    else "linear"
  }

  data_re <- .build_data(mbm_obj$SpeciesData.XY.Regional,
                         mbm_obj$Background.XY.Regional, rast_re)
  data_gl <- if (has_global)
    .build_data(mbm_obj$SpeciesData.XY.Global,
                mbm_obj$Background.XY.Global, rast_gl) else NULL

  tests_re <- lapply(vars_re, .test_rw2, data = data_re,
                     rw2.alpha = rw2.alpha, rw2.edf.min = rw2.edf.min)
  names(tests_re) <- vars_re

  tests_gl <- if (has_global) {
    t <- lapply(vars_gl, .test_rw2, data = data_gl,
                rw2.alpha = rw2.alpha, rw2.edf.min = rw2.edf.min)
    names(t) <- vars_gl; t
  } else NULL

  cov_effects_re <- lapply(tests_re,
                           function(x) .rw2_entry(x$suggest, rw2.u, rw2.alpha.pc))
  cov_effects_gl <- if (has_global)
    lapply(tests_gl, function(x) .rw2_entry(x$suggest, rw2.u, rw2.alpha.pc)) else NULL

  covariate.effects <- list(regional = cov_effects_re, default = "linear")
  if (has_global)
    covariate.effects <- c(list(global = cov_effects_gl), covariate.effects)

  rw2_tests_re <- data.frame(
    scale = "regional", variable = vars_re,
    edf     = sapply(tests_re, `[[`, "edf"),
    p_value = sapply(tests_re, `[[`, "p"),
    suggest = sapply(tests_re, `[[`, "suggest"),
    reason  = sapply(tests_re, `[[`, "reason"),
    stringsAsFactors = FALSE
  )
  rw2_tests_gl <- if (has_global) data.frame(
    scale = "global", variable = vars_gl,
    edf     = sapply(tests_gl, `[[`, "edf"),
    p_value = sapply(tests_gl, `[[`, "p"),
    suggest = sapply(tests_gl, `[[`, "suggest"),
    reason  = sapply(tests_gl, `[[`, "reason"),
    stringsAsFactors = FALSE
  ) else NULL
  rw2_tests <- rbind(rw2_tests_gl, rw2_tests_re)
  rownames(rw2_tests) <- NULL

  # ── 8. Assemble output ─────────────────────────────────────────────────────

  out <- list(
    create_mesh = list(
      edge   = c(edge_inner,   edge_outer),
      offset = c(offset_inner, offset_outer)
    ),
    MBM.Modelling_args = list(
      regional.pcprior.range = c(range_re, range.alpha),
      regional.pcprior.sigma = c(sigma_re, sigma.alpha),
      shared.pcprior.range   = if (has_global) c(range_shared, range.alpha) else NULL,
      shared.pcprior.sigma   = if (has_global) c(sigma_shared, sigma.alpha) else NULL,
      covariate.effects      = covariate.effects
    ),
    covariate_diagnostics = list(
      regional_tests = if(length(tests_re) > 0) do.call(rbind, tests_re) else NULL,
      global_tests   = if(length(tests_gl) > 0) do.call(rbind, tests_gl) else NULL
    ),
    spatial_metadata = list(
      crs_units              = crs_units,
      regional_diameter      = diam_re,
      range_re_expected      = range_re_expected,
      scale_separation_ratio = if (has_global) shared.range.ratio else NA_real_
    )
  )

  # ── 9. Verbose output ──────────────────────────────────────────────────────

  if (verbose) {
    cat("\n=================================================================\n")
    cat("  MBM.SuggestParams \u2014 Automatic Configuration Tool\n")
    cat("=================================================================\n\n")
    
    cat(">> SPATIAL UNITS DETECTED:", toupper(crs_units), "\n")
    cat("   Domain diameter:", round(diam_re, 2), crs_units, "\n")
    cat("   Expected regional range (20%):", round(range_re_expected, 2), crs_units, "\n\n")
    
    cat(">> FOR create_mesh():\n")
    cat("   edge   = c(", round(edge_inner, 3), ", ", round(edge_outer, 3), ")\n", sep="")
    cat("   offset = c(", round(offset_inner, 3), ", ", round(offset_outer, 3), ")\n\n", sep="")
    
    cat(">> FOR MBM.Modelling():\n")
    cat("   regional.pcprior.range = c(", round(range_re, 3), ", ", range.alpha, ")\n", sep="")
    cat("   regional.pcprior.sigma = c(", sigma_re, ", ", sigma.alpha, ")\n", sep="")
    if (has_global) {
      cat("   shared.pcprior.range   = c(", round(range_shared, 3), ", ", range.alpha, ")\n", sep="")
      cat("   shared.pcprior.sigma   = c(", sigma_shared, ", ", sigma.alpha, ")\n", sep="")
    }
    
    cat("\n>> COVARIATE RECOMMENDATIONS (covariate.effects):\n")
    cat("   * A GAM test (p-value <", rw2.alpha, "& EDF >", rw2.edf.min, ") was used to detect non-linearities.\n")
    
    vars_rw2 <- names(covariate.effects$regional)[sapply(covariate.effects$regional, function(x) is.list(x) && x$model == "rw2")]
    if(length(vars_rw2) > 0) {
      cat("   Suggested 'rw2' (non-linear):", paste(vars_rw2, collapse=", "), "\n")
    } else {
      cat("   All regional covariates suggested as 'linear'.\n")
    }
    cat("=================================================================\n")
  }
  
  invisible(out)
}
