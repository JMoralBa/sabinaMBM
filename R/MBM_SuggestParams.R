#' @name MBM.SuggestParams
#' @noRd
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
#'         (S_shared), derived sequentially from a single master spatial parameter.
#'   \item \strong{Covariate functional form} (\code{covariate.effects}): whether
#'         each covariate at each scale is better represented as a linear effect or
#'         a non-linear rw2 spline, assessed via an auxiliary GAM test.
#' }
#'
#' All outputs are suggestions based on theoretical guidelines and exploratory
#' data analysis. The user retains full control over the final parameter values
#' passed to \code{create_mesh()} and \code{MBM.Modelling()}.
#'
#' @param jmbm_obj An object of class \code{nsdm.vinput} produced by
#'   \code{sabinaNSDM::NSDM.SelectCovariates()}.
#' @param custom.range.Sre Numeric. Optional. The expected spatial autocorrelation
#'   range of the local fine-scale process (S_RE). \strong{Values must be provided
#'   in the exact same units as the native CRS of your environmental rasters.}
#'   For example, if your raster is in a metric projection (e.g., UTM) and you 
#'   expect a spatial dependence of 50 km, you must input \code{50000}. If your 
#'   raster is in an unprojected CRS (e.g., WGS84 degrees), you must input the 
#'   range in decimal degrees (e.g., \code{0.45}). If provided, this value overrides
#'   the automatic derivation based on \code{prior.range.fraction}.
#' @param prior.range.fraction Numeric scalar in (0, 1]. Assumed spatial
#'   autocorrelation range of the fine-scale process (S_RE) expressed as a
#'   fraction of the regional domain diameter. Default: \code{0.20} (i.e., the 
#'   effective range spans approximately one fifth of the study domain).
#' @param shared.range.ratio Numeric scalar >= 3. Assumed multiplicative factor by
#'   which the S_shared spatial process operates at a larger scale than S_RE.
#'   Default: \code{3}.
#' @param regional.sigma.u Numeric scalar. Upper bound for the PC prior on the
#'   marginal standard deviation of the regional field. Default: \code{1.5}.
#' @param shared.sigma.divisor Numeric scalar > 1. Divisor applied to
#'   \code{regional.sigma.u} to obtain the upper bound for the shared field.
#'   Default: \code{5} (yields a shared upper bound of 0.3).
#' @param range.alpha Numeric scalar in (0, 1). Probability for PC priors on
#'   spatial range: P(range < u) = \code{range.alpha}. Default: \code{0.05}.
#' @param sigma.alpha Numeric scalar in (0, 1). Probability for PC priors on
#'   marginal standard deviation: P(sigma > u) = \code{sigma.alpha}. Default: \code{0.01}.
#' @param rw2.alpha Numeric scalar in (0, 1). Significance threshold for the GAM
#'   non-linearity test (p-value of smooth term edf). Default: \code{0.01}.
#' @param rw2.edf.min Numeric scalar > 1. Minimum effective degrees of freedom
#'   (edf) of the GAM smooth term required to suggest rw2. Default: \code{1.5}.
#' @param rw2.u Numeric scalar. Value of \code{u} in the PC prior for rw2
#'   covariate effects: P(sigma_spline > u) = \code{rw2.alpha.pc}. Default: \code{0.5}.
#' @param rw2.alpha.pc Numeric scalar. Value of \code{alpha} in the PC prior for
#'   rw2 covariate effects. Default: \code{0.01}.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted summary
#'   of all suggestions with their rationale to the console.
#'
#' @return A standard named list containing four main elements:
#'   \describe{
#'     \item{\code{create_mesh}}{Named list with suggested \code{edge} and \code{offset}.}
#'     \item{\code{MBM.Modelling_args}}{Named list containing all the suggested arguments
#'       for \code{MBM.Modelling()}: \code{regional.pcprior.range}, \code{regional.pcprior.sigma},
#'       \code{shared.pcprior.range}, \code{shared.pcprior.sigma}, and \code{covariate.effects}.}
#'     \item{\code{covariate_diagnostics}}{Named list with \code{regional_tests} and
#'       \code{global_tests} containing GAM non-linearity test results.}
#'     \item{\code{spatial_metadata}}{Named list containing spatial context.}
#'   }
#'
#' @details
#' \strong{Sequential derivation of spatial parameters:}
#' To prevent user misspecification and double-smoothing in joint models, all spatial 
#' suggestions in this function are derived from a single master parameter: the expected 
#' autocorrelation range of the fine-scale process (\code{range_re_expected}). If not 
#' provided via \code{custom.range.Sre}, it is estimated as \code{prior.range.fraction} 
#' times the diameter of the regional bounding box. The full derivation cascade follows 
#' established guidelines in spatial statistics and INLA literature:
#' \enumerate{
#'   \item \strong{Expected regional range:} \code{range_re_expected = prior.range.fraction * diam_re}.
#'     The default fraction (0.20) assumes the effective range spans approximately one fifth 
#'     of the study domain, a standard INLA heuristic (Barber et al., 2016).
#'   \item \strong{S_RE range prior lower bound:} \code{rho_0_RE = range_re_expected / 10}.
#'     Optimal PC prior calibration establishes the lower bound at 1/10 of the expected 
#'     true range (Fuglstad et al., 2019).
#'   \item \strong{S_shared range prior lower bound:} \code{rho_0_shared = (shared.range.ratio * range_re_expected) / 10}.
#'     Applies the same calibration (Fuglstad et al., 2019) while guaranteeing statistical 
#'     identifiability by forcing the global field to operate at a significantly broader scale (Bakka et al., 2018).
#'   \item \strong{Mesh inner edge:} \code{edge[1] = range_re_expected / 5}.
#'     Inner triangles must be smaller than the spatial range (specifically < range/5) to capture 
#'     local spatial structure (Foster et al., 2024).
#'   \item \strong{Mesh outer edge:} \code{edge[2] = 3 * edge[1]}.
#'     The outer domain is made coarser to reduce computational burden (Lindgren et al., 2011; Anderson et al., 2022).
#'   \item \strong{Inner offset:} \code{offset[1] = edge[1]}.
#'     A practical heuristic creating a minimal 1-triangle buffer to avoid mesh collapse at the data boundary.
#'   \item \strong{Outer offset:} \code{offset[2] = range_re_expected}.
#'     A boundary buffer greater than or equal to the spatial range ensures that Neumann boundary variance 
#'     dissipates before reaching the study area (Krainski et al., 2019; Anderson et al., 2022).
#'   \item \strong{Regional sigma upper bound:} \code{sigma_re = 1.5} (fixed).
#'     On the logit scale, a standard deviation of 1.5 spans the approximate [0.05, 0.95] 
#'     probability range (Simpson et al., 2017).
#'   \item \strong{Shared sigma upper bound:} \code{sigma_shared = sigma_re / shared.sigma.divisor}.
#'     The shared field must be more tightly constrained to prevent variance cannibalization 
#'     (Simpson et al., 2017).
#' }
#'
#' \strong{Non-linearity test for rw2:}
#' For each covariate, an independent binomial GAM with a thin-plate spline is fitted 
#' (Wood, 2017). A non-linear random walk (rw2) is suggested only if the smooth term 
#' is statistically significant (p-value < \code{rw2.alpha}) AND the Estimated Degrees 
#' of Freedom exhibit true ecological curvature (edf > \code{rw2.edf.min}). This dual 
#' criterion prevents false positives in large datasets where negligible curvature might 
#' yield significant p-values.
#'
#' @references
#' Anderson, S. C., et al. (2022). sdmTMB: An R package for fast, flexible, and 
#' user-friendly generalized linear mixed effects models with spatial and spatiotemporal 
#' random fields. \emph{bioRxiv}.
#' 
#' Bakka, H., et al. (2018). Spatial modelling with R-INLA: A review. 
#' \emph{WIREs Computational Statistics}, 10, e1443.
#' 
#' Barber, X., et al. (2016). Incorporating spatial structure into inclusion probabilities 
#' for Bayesian variable selection in generalized linear models with the spike-and-slab prior. 
#' \emph{Journal of Applied Statistics}, 43, 2261–2279.
#'
#' Foster, S. D., et al. (2024). RISDM: species distribution modelling from multiple data 
#' sources in R. \emph{Ecography}, e06964.
#'
#' Fuglstad, G.-A., et al. (2019). Constructing priors that penalize the complexity 
#' of Gaussian random fields. \emph{Journal of the American Statistical Association}, 
#' 114, 445–452.
#'
#' Krainski, E., et al. (2019). \emph{Advanced Spatial Modeling with Stochastic 
#' Partial Differential Equations Using R and INLA}. Chapman and Hall/CRC.
#'
#' Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields 
#' and Gaussian Markov random fields: the stochastic partial differential equation approach. 
#' \emph{Journal of the Royal Statistical Society: Series B}, 73, 423–498.
#'
#' Simpson, D., et al. (2017). Penalising model component complexity: A principled, 
#' practical approach to constructing priors. \emph{Statistical Science}, 32, 1–28.
#'
#' Wood, S. N. (2017). \emph{Generalized Additive Models: An Introduction with R}. 
#' 2nd Edition. Chapman and Hall/CRC.
#'
#' @seealso \code{\link{create_mesh}}, \code{\link{MBM.Modelling}}
#'
#' @examples
#' \dontrun{
#' myParams <- MBM.SuggestParams(jmbm_obj = mySelvars)
#' 
#' myMesh <- create_mesh(
#'   nsdm_obj        = mySelvars, 
#'   edge            = myParams$create_mesh$edge, 
#'   offset          = myParams$create_mesh$offset, 
#'   boundary.method = "raster_mask"
#' )
#' 
#' myMod <- MBM.Modelling(
#'   jmbm_obj               = mySelvars,
#'   spde.mesh              = myMesh,
#'   regional.pcprior.range = myParams$MBM.Modelling_args$regional.pcprior.range,
#'   regional.pcprior.sigma = myParams$MBM.Modelling_args$regional.pcprior.sigma,
#'   shared.pcprior.range   = myParams$MBM.Modelling_args$shared.pcprior.range,
#'   shared.pcprior.sigma   = myParams$MBM.Modelling_args$shared.pcprior.sigma,
#'   covariate.effects      = myParams$MBM.Modelling_args$covariate.effects
#' )
#' }
#'
MBM.SuggestParams <- function(jmbm_obj,
                              custom.range.Sre = NULL,
                              prior.range.fraction = 0.20,
                              shared.range.ratio = 3,
                              regional.sigma.u = 1.5,
                              shared.sigma.divisor = 5,
                              range.alpha = 0.05,
                              sigma.alpha = 0.01,
                              rw2.alpha = 0.01,
                              rw2.edf.min = 1.5,
                              rw2.u = 0.5,
                              rw2.alpha.pc = 0.01,
                              verbose = TRUE) {

  # checks
  if (!inherits(jmbm_obj, "nsdm.vinput"))
    .stop("'jmbm_obj' must be of class 'nsdm.vinput'.
          Please use sabinaNSDM::NSDM.SelectCovariates() to obtain it.\n")
  if (!is.numeric(prior.range.fraction) || prior.range.fraction <= 0 || prior.range.fraction > 1) 
    .stop("'prior.range.fraction' must be in (0, 1]. 
          Default 0.20 is the INLA standard assumption.")
  if (!is.numeric(shared.range.ratio) || shared.range.ratio <= 1) 
    .stop("'shared.range.ratio' must be > 1 to maintain structural scale hierarchy.
          Values below 3 risk double-smoothing and weak scale separation")
  if (shared.range.ratio < 3) 
    .warn("'shared.range.ratio' < 3 risks double-smoothing and weak scale separation (Bakka et al., 2018).")
  if (!is.numeric(shared.sigma.divisor) || shared.sigma.divisor <= 1) 
    .stop("'shared.sigma.divisor' must be > 1: shared sigma must be more constrained than regional sigma.")
  #
  if (!is.numeric(range.alpha) || range.alpha <= 0 || range.alpha >= 1) 
    .stop("'range.alpha' must be a valid probability in (0, 1).")
  if (!is.numeric(sigma.alpha) || sigma.alpha <= 0 || sigma.alpha >= 1) 
    .stop("'sigma.alpha' must be a valid probability in (0, 1).")
  if (!is.numeric(rw2.alpha) || rw2.alpha <= 0 || rw2.alpha >= 1) 
    .stop("'rw2.alpha' must be a valid probability in (0, 1).")
  #
  if (!is.numeric(rw2.edf.min) || rw2.edf.min <= 1) 
    .stop("'rw2.edf.min' must be > 1 (a linear effect has edf = 1).")


  # model type
  has_global <- !is.null(jmbm_obj$Selected.Variables.Global) &&
                length(jmbm_obj$Selected.Variables.Global) > 0

  vars_gl <- if (has_global) jmbm_obj$Selected.Variables.Global else character(0)
  vars_re <- jmbm_obj$Selected.Variables.Regional


  # detect CRS
  rast_re <- terra::unwrap(jmbm_obj$IndVar.Regional.Selected)
  rast_gl <- if (has_global) terra::unwrap(jmbm_obj$IndVar.Global.Selected) else NULL
  
  is_lonlat <- isTRUE(sf::st_is_longlat(sf::st_crs(rast_re)))
  if(is_lonlat) {
    .warn(paste0("The input rasters use an unprojected geographic CRS (longitude/latitude).\n",
                 "   Distances and mesh geometry will be computed in degrees, which distorts true spatial relationships.\n",
                 "   sabinaMBM will proceed, but projecting your covariates to a metric CRS (e.g., UTM) is strongly recommended."))
    crs_units <- "degrees (unprojected)"
  } else {
    crs_units <- "native CRS units (projected)"
  }


  # Domain diameters 
  .diameter <- function(coords_df) {
    bb <- sf::st_bbox( sf::st_as_sf(as.data.frame(coords_df), coords = c("x", "y")) )
    sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  }

  coords_re <- rbind(
    jmbm_obj$SpeciesData.XY.Regional[, c("x", "y")],
    jmbm_obj$Background.XY.Regional[, c("x", "y")]
  )
  diam_re <- .diameter(coords_re)

  diam_gl <- if (has_global) {
    coords_gl <- rbind(
      jmbm_obj$SpeciesData.XY.Global[, c("x", "y")],
      jmbm_obj$Background.XY.Global[, c("x", "y")]
    )
    .diameter(coords_gl)
  } else NA_real_


  # sequential spatial parameter derivation
  if (!is.null(custom.range.Sre)) {
    range_re_expected <- custom.range.Sre
  } else {
    range_re_expected <- prior.range.fraction * diam_re
  }


  # PC prior lower bounds (Fuglstad et al. 2019)
  range_re <- range_re_expected / 10
  range_shared <- if (has_global) (shared.range.ratio * range_re_expected) / 10 else NULL


  # Sigma prior upper bounds
  sigma_re <- regional.sigma.u
  sigma_shared <- if (has_global) round(sigma_re / shared.sigma.divisor, 2) else NULL


  # Mesh geometry
  edge_inner <- range_re_expected / 5
  edge_outer <- 3 * edge_inner
  offset_inner <- edge_inner
  offset_outer <- range_re_expected


  # Non-linearity test
  # Binomial GAM + thin-plate spline (k=5)   #@@@JMB conservative???
  # rw2 suggested when p < rw2.alpha AND edf > rw2.edf.min.
  # Dual filter avoids large-N spurious significance and negligible curvature
  .build_data <- function(xy_pres, xy_bg, rast) {
    coords <- rbind(
      cbind(xy_pres[, c("x", "y")], resp = 1L),
      cbind(xy_bg[, c("x", "y")], resp = 0L)
    )
    vals <- terra::extract(rast, coords[, c("x", "y")])[, -1, drop = FALSE]
    cbind(data.frame(resp = coords$resp), vals)
  }

  .test_rw2 <- function(data, var_name, rw2.alpha, rw2.edf.min) {
    keep <- stats::complete.cases(data[, c("resp", var_name)])
    d <- data[keep, , drop = FALSE]
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
    sm <- summary(fit)
    edf <- sm$s.table[1, "edf"]
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

  data_re <- .build_data(jmbm_obj$SpeciesData.XY.Regional,
                         jmbm_obj$Background.XY.Regional, rast_re)
  data_gl <- if (has_global)
    .build_data(jmbm_obj$SpeciesData.XY.Global,
                jmbm_obj$Background.XY.Global, rast_gl) else NULL

  tests_re <- lapply(vars_re, .test_rw2, data = data_re,
                     rw2.alpha = rw2.alpha, rw2.edf.min = rw2.edf.min)
  names(tests_re) <- vars_re

  tests_gl <- if (has_global) {
    t <- lapply(vars_gl, .test_rw2, data = data_gl,
                rw2.alpha = rw2.alpha, rw2.edf.min = rw2.edf.min)
    names(t) <- vars_gl; t
  } else NULL

  rw2_tests_re <- data.frame(
    scale = "regional", variable = vars_re,
    edf = sapply(tests_re, `[[`, "edf"),
    p_value = sapply(tests_re, `[[`, "p"),
    suggest = sapply(tests_re, `[[`, "suggest"),
    reason = sapply(tests_re, `[[`, "reason"),
    stringsAsFactors = FALSE
  )
  rw2_tests_gl <- if (has_global) data.frame(
    scale = "global", variable = vars_gl,
    edf = sapply(tests_gl, `[[`, "edf"),
    p_value = sapply(tests_gl, `[[`, "p"),
    suggest = sapply(tests_gl, `[[`, "suggest"),
    reason  = sapply(tests_gl, `[[`, "reason"),
    stringsAsFactors = FALSE
  ) else NULL
  rw2_tests <- rbind(rw2_tests_gl, rw2_tests_re)
  rownames(rw2_tests) <- NULL

  covariate.effects <- list(
    regional = lapply(tests_re, function(res) .rw2_entry(res$suggest, rw2.u, rw2.alpha.pc)),
    global = if(has_global) lapply(tests_gl, function(res) .rw2_entry(res$suggest, rw2.u, rw2.alpha.pc)) else NULL,
    default  = "linear"
  )

  # output
  out <- list(
    create_mesh = list(
      edge = c(edge_inner, edge_outer),
      offset = c(offset_inner, offset_outer)
    ),
    MBM.Modelling_args = list(
      regional.pcprior.range = c(range_re, range.alpha),
      regional.pcprior.sigma = c(sigma_re, sigma.alpha),
      shared.pcprior.range = if (has_global) c(range_shared, range.alpha) else NULL,
      shared.pcprior.sigma = if (has_global) c(sigma_shared, sigma.alpha) else NULL,
      covariate.effects = covariate.effects
    ),
    covariate_diagnostics = list(
      regional_tests = if(nrow(rw2_tests_re) > 0) rw2_tests_re else NULL,
      global_tests = if(!is.null(rw2_tests_gl) && nrow(rw2_tests_gl) > 0) rw2_tests_gl else NULL
    ),
    spatial_metadata = list(
      crs_units = crs_units,
      regional_diameter = diam_re,
      range_re_expected = range_re_expected,
      scale_separation_ratio = if (has_global) shared.range.ratio else NA_real_
    )
  )

  # verbose output
  if (verbose) {
    cat(paste0(
      "\n=================================================================\n",
      "  MBM.SuggestParams // Suggested configuration\n",
      "=================================================================\n\n",
      
      ">> SPATIAL UNITS DETECTED: ", toupper(crs_units), "\n",
      "   Domain diameter: ", round(diam_re, 2), " ", crs_units, "\n",
      sprintf("   Expected regional range (%g%%): ", prior.range.fraction * 100),
      round(range_re_expected, 2), " ", crs_units, "\n\n",
      
      ">> FOR create_mesh():\n",
      "   edge   = c(", round(edge_inner, 3), ", ", round(edge_outer, 3), ")\n",
      "   offset = c(", round(offset_inner, 3), ", ", round(offset_outer, 3), ")\n\n",
      
      ">> FOR MBM.Modelling():\n",
      "   regional.pcprior.range = c(", round(range_re, 3), ", ", range.alpha, ")\n",
      "   regional.pcprior.sigma = c(", sigma_re, ", ", sigma.alpha, ")\n",
      if (has_global) paste0(
        "   shared.pcprior.range   = c(", round(range_shared, 3), ", ", range.alpha, ")\n",
        "   shared.pcprior.sigma   = c(", sigma_shared, ", ", sigma.alpha, ")\n"),
      
      "\n>> COVARIATE RECOMMENDATIONS (covariate.effects):\n",
      "   * A GAM test (p-value < ", rw2.alpha, " & EDF > ", rw2.edf.min, ") was used to detect non-linearities.\n"
    ))

    if (has_global && !is.null(covariate.effects$global)) {
      vars_rw2_gl <- names(covariate.effects$global)[sapply(covariate.effects$global, function(x) is.list(x) && x$model == "rw2")]
      if (length(vars_rw2_gl) > 0) {
        cat("   Global scale — suggested 'rw2' (non-linear):", paste(vars_rw2_gl, collapse = ", "), "\n")
      } else {
        cat("   Global scale — all covariates suggested as 'linear'.\n")
      }
    }

    vars_rw2_re <- names(covariate.effects$regional)[sapply(covariate.effects$regional, function(x) is.list(x) && x$model == "rw2")]
    if (length(vars_rw2_re) > 0) {
      cat("   Regional scale — suggested 'rw2' (non-linear):", paste(vars_rw2_re, collapse = ", "), "\n")
    } else {
      cat("   Regional scale — all covariates suggested as 'linear'.\n")
    }
    
    cat("=================================================================\n")
    cat("Parameters locked. Ready to initialize MBM.Modelling().\n\n")
  }
  
  invisible(out)
}

