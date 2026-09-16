#' @name plot.jmbm.inlabru
#'
#' @title Plot predictions and spatial fields from a fitted MBM model
#'
#' @description Generates diagnostic plots and spatial prediction maps for objects of class \code{jmbm.inlabru}. Supports current and future suitability maps, spatial field visualizations, posterior marginals, and calibration diagnostics.
#'
#' @param x An object of class "jmbm.inlabru".
#' @param which Component to plot. One of:
#'   \itemize{
#'     \item \code{"pred"}: current suitability prediction (default).
#'     \item \code{"pred_Sre"}: Sre field (requires SPDE).
#'     \item \code{"pred_Sshared"}: Sshared field (requires SPDE).
#'     \item \code{"hyperparams"}: posterior marginals of SPDE hyperparameters with PC-priors overlaid.
#'     \item \code{"intercepts"}: posterior marginals of IGlobal and IRegional intercepts.
#'     \item \code{"fixed"}: posterior marginals of all fixed effect coefficients, back-transformed to original covariate scale.
#'     \item \code{"correlogram"}: empirical correlogram of regional model residuals against pairwise distance, used to check for residual spatial autocorrelation not captured by the spatial field(s). Requires Sre or Sshared.
#'     \item \code{"semivariogram"}: empirical semivariogram of regional model residuals, complementary diagnostic to the correlogram for residual spatial structure. Requires Sre or Sshared.
#'     \item \code{"Srefields"}: posterior mean maps of the fitted spatial field(s) (Sre and/or Sshared, whichever are present in the model).
#'     \item \code{"hist"}: histogram of regional model residuals, centred at 0 (or at the mean residual when \code{coupling.intercept = "ordered_hierarchical"}, since the soft-constraint intercept does not force a zero-centred residual distribution).
#'     \item \code{"qq"}: normal QQ-plot of regional model residuals, for visual assessment of normality.
#'     \item \code{"pit"}: histogram of Probability Integral Transform values. A uniform distribution indicates good calibration.
#'     \item \code{"new.projections"}(first scenario), or scenario name as string (e.g. \code{"MRI_ESM2_0_2070_SSP585"}): future/alternative scenario prediction.
#'   }
#' @param layer Raster layer to display (default = "mean"). Available: "mean", "sd", "q0.025", "q0.5", "q0.975", "median", "sd.mc_std_err", "mean.mc_std_err".
#' @param palette Color palette passed to ggplot2::scale_fill_distiller(). Default = "Spectral".
#' @param title Optional plot title. If NULL, a default title is generated as "Species - type - scope - layer".
#' @param legend_title Optional legend title. If NULL, it is automatically inferred from the model family and the layer.
#' @param ... Additional graphical arguments passed to methods.
#'
#' @return A ggplot object showing the selected MBM output layer.
#'
#' @examples
#' \dontrun{
#' # Default: current prediction (mean layer)
#' plot(myModel)
#'
#' ## Sre field
#' plot(myModel, which = "pred_Sre", layer = "sd")
#'
#' ## Scenario by index or name
#' plot(myModel, which = "new.projections[[1]]")
#' plot(myModel, which = "scenario1", layer = "q0.975")
#' }
#'
#' @seealso \code{\link{MBM.Modelling}}, \code{\link{summary.jmbm.inlabru}}
#'
#' @export
#' @method plot jmbm.inlabru
plot.jmbm.inlabru <- function(x,
                              which = "pred",
                              layer = "mean",
                              palette = "Spectral",
                              title = NULL,
                              legend_title = NULL,
                              ...) {

  stopifnot(inherits(x, "jmbm.inlabru"))

  scope_label <- NULL
  r <- NULL

  if (identical(which, "pred")) {
    # current
    r <- x$current.projections$pred
    scope_label <- "Current"

  } else if(identical(which, "pred_Sre")) {
    r <- x$current.projections$pred_Sre
    if(is.null(r)) .stop("No Sre field ('pred_Sre') in this model. Refit with `regional.pcprior.range` and `regional.pcprior.sigma`.")
    scope_label <- "Sre field"

  } else if(identical(which, "pred_Sshared")) {
    r <- x$current.projections$pred_Sshared
    if(is.null(r)) .stop("No Sshared field ('pred_Sshared') in this model. Refit with `shared.pcprior.range` and `shared.pcprior.sigma`.")
    scope_label <- "Sshared field"

  } else if(identical(which, "hyperparams")) {
    hp <- x$diagnostic_data$hyperparams
    if(is.null(hp) || is.null(hp$posteriors) || nrow(hp$posteriors) == 0)
      .stop("No hyperparameter marginals in this model. Refit with SPDE priors.")
    species <- gsub("\\.", " ", x$Species.Name)
    return(.plot_hyperparams_jmbm(hp$posteriors, hp$prior_ticks,
                                   title = if(!is.null(title)) title else paste(species, "| Hyperparameters: posterior marginals")))

  } else if(identical(which, "intercepts")) {
    marg_r <- x$marginals$random
    if(is.null(marg_r) || length(marg_r) == 0)
      .stop("No intercept marginals in this model.")
    species <- gsub("\\.", " ", x$Species.Name)
    coupling_label <- if(!is.null(x$args$coupling.intercept)) x$args$coupling.intercept else "regional only"
    int_names <- intersect(c("IGlobal", "IRegional"), names(marg_r))
    if(length(int_names) == 0)
      .stop("IGlobal/IRegional not found in marginals$random.")
    df_list <- lapply(int_names, function(nm) {
      m <- marg_r[[nm]][[1]]  # iid has one level
      sm <- INLA::inla.smarginal(m)
      data.frame(x = sm$x, y = sm$y, intercept = nm)
    })
    df <- do.call(rbind, df_list)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = intercept)) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::scale_colour_manual(
        values = c(IGlobal = "#E74C3C", IRegional = "#2E75B6"),
        name = NULL) +
      ggplot2::labs(
        title = if(!is.null(title)) title else paste(species, "| Intercept posteriors"),
        subtitle = paste("coupling.intercept =", coupling_label),
        x = "Value", y = "Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = c(0.88, 0.88),
        legend.background = ggplot2::element_rect(fill = "white", colour = NA),
        legend.key.size = ggplot2::unit(0.4, "cm"),
        legend.text = ggplot2::element_text(size = 9))

    return(p)

  } else if(identical(which, "fixed")) {
    marg_f <- x$marginals$fixed
    sp_glo <- if (!is.null(x$scale_params_glo)) x$scale_params_glo else x$scale_params
    sp_reg <- if (!is.null(x$scale_params_reg)) x$scale_params_reg else x$scale_params
    if(is.null(marg_f) || length(marg_f) == 0)
      .stop("No fixed effect marginals in this model.")
    species <- gsub("\\.", " ", x$Species.Name)
    cp_label <- if(!is.null(x$args$coupling.covariates)) x$args$coupling.covariates else "NULL"
    # back-transform marginals to original covariate scale
    df_list <- lapply(names(marg_f), function(nm) {
      base_var <- gsub("GL$|RE$|RE_oh$|GL_glo_res$|RE_reg_anom$", "", nm)
      is_global <- grepl("GL$|GL_glo_res$", nm)
      sp <- if(is_global) sp_glo else sp_reg
      m <- marg_f[[nm]]
      if(!is.null(sp) && base_var %in% names(sp) && isTRUE(sp[[base_var]]$sd > 0)) {
        sd_x <- sp[[base_var]]$sd
        m <- INLA::inla.tmarginal(function(x) x / sd_x, m)
      }
      sm <- INLA::inla.smarginal(m)
      data.frame(x = sm$x, y = sm$y, coef = nm)
    })
    df <- do.call(rbind, df_list)
    # significance: CI do not cross 0
    signif_coefs <- names(marg_f)[vapply(names(marg_f), function(nm) {
      base_var <- gsub("GL$|RE$|RE_oh$|GL_glo_res$|RE_reg_anom$", "", nm)
      is_global <- grepl("GL$|GL_glo_res$", nm)
      sp <- if(is_global) sp_glo else sp_reg
      m <- marg_f[[nm]]
      if(!is.null(sp) && base_var %in% names(sp) && isTRUE(sp[[base_var]]$sd > 0))
        m <- INLA::inla.tmarginal(function(x) x / sp[[base_var]]$sd, m)
      q025 <- INLA::inla.qmarginal(0.025, m)
      q975 <- INLA::inla.qmarginal(0.975, m)
      q025 * q975 > 0
    }, logical(1))]
    df$significant <- df$coef %in% signif_coefs

    # label non-significant panels
    ns_labels <- df[!df$significant & !duplicated(df$coef), c("coef", "significant")]
    ns_labels <- ns_labels[!ns_labels$significant, , drop = FALSE]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = significant)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_colour_manual(values = c("TRUE" = "#2E75B6", "FALSE" = "#95A5A6"),
                                   guide = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::facet_wrap(~ coef, scales = "free") +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title = if(!is.null(title)) title else paste(species, "| Fixed effects posteriors (original scale)"),
        subtitle = paste0("coupling.covariates = ", cp_label,
                         " | Coefficients in original covariate units (back-transformed)",
                         " | Grey = non-significant (95% CI)"),
        x = "Coefficient value", y = "Density") +
      ggplot2::theme_minimal()

    # add NS annotation to non-significant panels
    if(nrow(ns_labels) > 0) {
      p <- p + ggplot2::geom_text(
        data = ns_labels,
        ggplot2::aes(label = "NS", x = -Inf, y = Inf),
        hjust = -0.2, vjust = 1.5,
        colour = "#95A5A6", size = 3.5, fontface = "italic",
        inherit.aes = FALSE)
    }
    return(p)

  } else if(identical(which, "correlogram")) {
    cor_df <- x$diagnostic_data$correlogram
    if(is.null(cor_df) || nrow(cor_df) == 0)
      .stop("No correlogram available. Refit with Sre or Sshared SPDE.")
    species <- gsub("\\.", " ", x$Species.Name)
    range_eff <- attr(cor_df, "range_eff")
    return(.plot_correlogram_jmbm(cor_df, range_eff,
                                   title = if(!is.null(title)) title else paste(species, "| Residual correlogram")))

  } else if(identical(which, "semivariogram")) {
    sv_df <- x$diagnostic_data$semivariogram
    if(is.null(sv_df) || nrow(sv_df) == 0)
      .stop("No semivariogram available. Refit with Sre or Sshared SPDE.")
    species <- gsub("\\.", " ", x$Species.Name)
    range_res_mean <- attr(sv_df, "range_res_mean")
    range_lat_mean <- attr(sv_df, "range_lat_mean")
    return(.plot_semivariogram_jmbm(sv_df, range_res_mean, range_lat_mean,
                                     title = if(!is.null(title)) title else paste(species, "| Empirical semivariogram")))

  } else if(identical(which, "Srefields")) {
    species <- gsub("\\.", " ", x$Species.Name)
    ordered_hierarchical <- !is.null(x$args$coupling.intercept) && x$args$coupling.intercept == "ordered_hierarchical"
    return(.plot_srefields_jmbm(x$current.projections$pred_Sre, x$current.projections$pred_Sshared,
                                 ordered_hierarchical,
                                 title = if(!is.null(title)) title else paste(species, "| Sre fields (posterior mean)")))

  } else if(identical(which, "hist")) {
    df_r <- x$diagnostic_data$hist
    if(is.null(df_r) || nrow(df_r) == 0)
      .stop("No residual data available for histogram.")
    species <- gsub("\\.", " ", x$Species.Name)
    return(.plot_hist_jmbm(df_r, title = if(!is.null(title)) title else paste(species, "| Residual histogram")))

  } else if(identical(which, "qq")) {
    df_qq <- x$diagnostic_data$qq
    if(is.null(df_qq) || nrow(df_qq) == 0)
      .stop("No residual data available for QQ-plot.")
    return(.plot_qq_jmbm(df_qq, title = title))

  } else if(identical(which, "pit")) {
    pit <- x$pit_values
    if(is.null(pit) || length(pit) == 0)
      .stop("No PIT values in this model. Refit with control.compute = list(cpo = TRUE).")
    species <- gsub("\\.", " ", x$Species.Name)
    ks_p <- tryCatch(
      suppressWarnings(ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_)
    ks_label <- if(is.finite(ks_p)) paste0("KS p-value = ", round(ks_p, 3)) else ""
    df_pit <- data.frame(pit = pit)
    p <- ggplot2::ggplot(df_pit, ggplot2::aes(x = pit)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 20, fill = "#2E75B6", colour = "white", alpha = 0.8) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "#E74C3C", linewidth = 0.8) +
      ggplot2::annotate("text", x = 0.95, y = Inf, label = ks_label,
                        hjust = 1, vjust = 1.5, size = 3.5, colour = "grey40") +
      ggplot2::labs(
        title = if(!is.null(title)) title else paste(species, "| PIT calibration histogram"),
        subtitle = "Dashed line = uniform (perfect calibration)",
        x = "PIT value", y = "Density") +
      ggplot2::theme_minimal()
    return(p)

  } else {
    # new.projections
    np <- x$new.projections
    if(is.null(np) || !length(np)) {
      .stop("No 'new.projections' in the object.")
    }

    get_np_by <- function(id) {
      nms <- names(np)
      if(is.numeric(id)) {
        id <- as.integer(id)
        if(id < 1 || id > length(np)) {
          .stop("Index out of range in new.projections.")
        }
        list(obj = np[[id]], name = nms[id])
      } else {
        if(!(id %in% nms)) {
          .stop(sprintf("Scenario '%s' not found in new.projections.", id))
        }
        list(obj = np[[id]], name = id)
      }
    }

    # first list obj, [[index]] / [[name]], or direct name
    if(identical(which, "new.projections")) {
      pick <- get_np_by(1)
    } else if(grepl("^new\\.projections\\[\\[.*\\]\\]$", which)) {
      inside <- sub("^new\\.projections\\[\\[(.*)\\]\\]$", "\\1", which)
      if(grepl("^[0-9]+$", inside)) inside <- as.integer(inside)
      pick <- get_np_by(inside)
    } else if(which %in% names(np)) {
      pick <- get_np_by(which)
    } else {
      stop("`which` must be 'pred', 'pred_Sre', 'new.projections', 'new.projections[[...]]', or an exact scenario name in new.projections.")
    }

    r <- pick$obj
    scope_label <- pick$name
  }

  rr <- if(inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r
  stopifnot(inherits(rr, "SpatRaster"))

  if(!layer %in% names(rr)) {
    .stop(paste0("Layer '", layer, "' does not exist. Available layers: ",
                paste(names(rr), collapse = ", ")))
  }
  r_show <- rr[[layer]]

  # legend title
  fam <- tolower(x$args$family)
  ln  <- tolower(layer)

  legend_prob <- function(ln) switch(ln,
    "mean"="Occurrence probability",
    "median"="Median probability",
    "q0.5"="Median probability",
    "sd"="Uncertainty (SD)",
    "q0.025"="Lower 95% quantile",
    "q0.975"="Upper 95% quantile",
    "sd.mc_std_err"="MC SE of SD",
    "mean.mc_std_err"="MC SE of mean",
    "Suitability")

  legend_eta <- function(ln) switch(ln,
    "mean"="Linear predictor (η)",
    "median"="Median η",
    "q0.5"="Median η",
    "sd"="Uncertainty of η (SD)",
    "q0.025"="Lower 95% quantile of η",
    "q0.975"="Upper 95% quantile of η",
    "sd.mc_std_err"="MC SE of SD (η)",
    "mean.mc_std_err"="MC SE of mean (η)",
    "Linear predictor (η)")

  legend_int <- function(ln) switch(ln,
    "mean"="Intensity (λ)",
    "median"="Median intensity (λ)",
    "q0.5"="Median intensity (λ)",
    "sd"="Uncertainty of intensity (SD)",
    "q0.025"="Lower 95% quantile of intensity",
    "q0.975"="Upper 95% quantile of intensity",
    "sd.mc_std_err"="MC SE of SD (λ)",
    "mean.mc_std_err"="MC SE of mean (λ)",
    "Intensity (λ)")

  auto_legend <- if(scope_label %in% c("Sre field", "Sshared field")) {
    legend_eta(ln)
  } else if(grepl("binomial", fam)) {
    legend_prob(ln)
  } else if(grepl("poisson|cp", fam)) {
    legend_int(ln)
  } else if(grepl("gaussian", fam)) {
    switch(ln,
      "mean"="Expected value",
      "median"="Median value",
      "q0.5"="Median value",
      "sd"="Uncertainty (SD)",
      "q0.025"="Lower 95% quantile",
      "q0.975"="Upper 95% quantile",
      "sd.mc_std_err"="MC SE of SD",
      "mean.mc_std_err"="MC SE of mean",
      "Expected value")
  } else {
    legend_prob(ln)
  }

  wrap_legend_title <- function(s, width = 22) {
    s <- as.character(s)
    paste(strwrap(s, width = width), collapse = "\n")
  }
  lt <- wrap_legend_title(if(is.null(legend_title)) auto_legend else legend_title)

  # plot title 
  species <- x$Species.Name
  plot_title <- if(is.null(title)) {
    paste(gsub("\\.", " ", species), "|", scope_label, "|", layer)
  } else {
    title
  }

  # plot
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r_show) +
    ggplot2::scale_fill_distiller(palette = palette, name = lt, na.value = "transparent") +
    ggplot2::labs(x = "Longitude", y = "Latitude", title = plot_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10))
}



# Auxiliary plot builders
# -----------------------

.plot_hyperparams_jmbm <- function(post_df, prior_ticks, title = NULL) {
  param_labels <- c(
    "Range for Sre" = "Sre field: Range",
    "Stdev for Sre" = "Sre field: Sigma",
    "Range for Sshared" = "Sshared field: Range",
    "Stdev for Sshared" = "Sshared field: Sigma",
    "Precision for IGlobal" = "IGlobal: Precision",
    "Precision for IRegional" = "IRegional: Precision",
    "Beta for IRegional" = "IRegional: Beta (copy)"
  )
  relabel <- function(x) ifelse(x %in% names(param_labels), param_labels[x], x)
  post_df$par <- relabel(post_df$par)
  if(!is.null(prior_ticks) && nrow(prior_ticks) > 0) {
    prior_ticks$par <- relabel(prior_ticks$par)
  }

  p <- ggplot2::ggplot(post_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(linewidth = 0.6, colour = "#1a5276") +
    ggplot2::facet_wrap(~par, scales = "free", ncol = 2)
  if(!is.null(prior_ticks) && nrow(prior_ticks) > 0) {
    p <- p +
      ggplot2::geom_vline(data = prior_ticks, ggplot2::aes(xintercept = x),
                          linetype = "dashed", linewidth = 0.5, colour = "#c0392b", alpha = 0.7) +
      ggplot2::geom_text(data = prior_ticks,
                         ggplot2::aes(x = x, y = y, label = paste0("PC prior (u) = ", round(x, 2))),
                         vjust = -0.4, hjust = 1, size = 2.8, colour = "#c0392b", angle = 90)
  }
  p + ggplot2::labs(
        title = if(!is.null(title)) title else "Hyperparameters: posterior marginals",
        subtitle = "Red dashed lines = PC prior (u)",
        y = "Density", x = "Value") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      strip.text = ggplot2::element_text(face = "bold", size = 9, colour = "#2c3e50"),
      axis.title = ggplot2::element_text(size = 9, colour = "#2c3e50"),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#d5d8dc"))
}

.plot_srefields_jmbm <- function(pred_Sre, pred_Sshared, ordered_hierarchical = FALSE, title = NULL) {
  ras_to_df <- function(r, nm) {
    rr <- terra::unwrap(r)[["mean"]]
    df <- terra::as.data.frame(rr, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", "mean")
    df$which <- nm
    df
  }
  maps_df <- data.frame()
  if(!is.null(pred_Sre)) maps_df <- rbind(maps_df, ras_to_df(pred_Sre, "Sre field"))
  if(!is.null(pred_Sshared)) maps_df <- rbind(maps_df, ras_to_df(pred_Sshared, "Sshared field"))

  if(nrow(maps_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = if(!is.null(title)) title else "Sre fields (posterior mean)",
          subtitle = "No Sre or Sshared fields present in the model") +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
          plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"))
    )
  }

  zlim <- range(maps_df$mean, na.rm = TRUE)

  ggplot2::ggplot(maps_df, ggplot2::aes(x = x, y = y, fill = mean)) +
    ggplot2::geom_raster(na.rm = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "#c0392b", mid = "white", high = "#1a5276",
      midpoint = 0, limits = zlim, na.value = "white",
      oob = scales::squish) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::facet_wrap(~which, ncol = 2, scales = "fixed") +
    ggplot2::labs(
      title = if(!is.null(title)) title else "Sre fields (posterior mean)",
      subtitle = "Red = below average, Blue = above average (centered at zero)",
      fill = "Posterior\nmean") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      strip.text = ggplot2::element_text(face = "bold", size = 9, colour = "#2c3e50"),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(), plot.background = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(0.3, "cm"), legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.title = ggplot2::element_text(size = 9, colour = "#2c3e50"),
      legend.text = ggplot2::element_text(size = 8))
}

.plot_correlogram_jmbm <- function(cor_df, range_eff = NULL, title = NULL) {
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = dist_mid, y = rho)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(na.rm = TRUE, size = 1.2, colour = "#1a5276") +
    ggplot2::geom_line(na.rm = TRUE, colour = "#1a5276", linewidth = 0.6)
  if(!is.null(range_eff) && is.finite(range_eff)) {
    p <- p +
      ggplot2::geom_vline(xintercept = range_eff, linetype = "dashed", colour = "#c0392b", alpha = 0.7) +
      ggplot2::annotate("text", x = range_eff, y = max(cor_df$rho, na.rm = TRUE),
                       label = "Model range", angle = 90, vjust = -0.8, hjust = 0.9,
                       colour = "#c0392b", size = 3)
  }
  p + ggplot2::labs(
        title = if(!is.null(title)) title else "Residual correlogram (residuals: obs \u2212 fitted mean)",
        subtitle = if(!is.null(range_eff) && is.finite(range_eff))
          paste0("Model range \u2248 ", round(range_eff, 3), " (map units)\n(distance where correlation vanishes)")
        else "No Sre/Sshared field: full extent shown",
        x = "Distance (map units)", y = "Residual correlation") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      strip.text = ggplot2::element_text(face = "bold", size = 9, colour = "#2c3e50"),
      axis.title.x = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, colour = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#d5d8dc"))
}

.plot_hist_jmbm <- function(df_r, title = NULL) {
  center_label <- attr(df_r, "center_label")
  line_x <- attr(df_r, "line_x")
  ggplot2::ggplot(df_r, ggplot2::aes(x = resid)) +
    ggplot2::geom_histogram(bins = 30, fill = "#1a5276", colour = "white", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#c0392b", linewidth = 0.4) +
    ggplot2::annotate("text", x = line_x, y = Inf, label = center_label,
                     angle = 90, vjust = -0.8, hjust = 1.2, size = 2.8, colour = "#c0392b") +
    ggplot2::labs(
      title = if(!is.null(title)) title else "Residual histogram (residuals: obs \u2212 fitted mean)",
      subtitle = paste0("Distribution of residuals\n(dashed line = ", center_label, ")"),
      x = "Residuals", y = "Frequency") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, colour = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#d5d8dc"),
      plot.background = ggplot2::element_blank())
}

.plot_qq_jmbm <- function(df_qq, title = NULL) {
  ggplot2::ggplot(df_qq, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point(colour = "#1a5276", size = 1.3, alpha = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#c0392b", linewidth = 0.4) +
    ggplot2::labs(
      title = if(!is.null(title)) title else "Residual QQ-plot",
      subtitle = "Residuals vs. theoretical quantiles\n(dashed = normal expectation)",
      x = "Theoretical quantiles (Normal)", y = "Sample quantiles (residuals)") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, colour = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#d5d8dc"),
      plot.background = ggplot2::element_blank())
}

.plot_semivariogram_jmbm <- function(sv_df, range_res_mean = NULL, range_lat_mean = NULL, title = NULL) {
  p <- ggplot2::ggplot(sv_df, ggplot2::aes(x = dist, y = gamma)) +
    ggplot2::geom_point(colour = "#1a5276", size = 1.5, alpha = 0.8) +
    ggplot2::geom_line(colour = "#1a5276", linewidth = 0.6, alpha = 0.8) +
    ggplot2::labs(
      title = if(!is.null(title)) title else "Empirical semivariogram (residuals: obs \u2212 fitted mean)",
      subtitle = " ",
      x = "Distance (map units)", y = "Semivariance \u03b3(h)") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, colour = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9, colour = "#2c3e50", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8, colour = "#2c3e50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#d5d8dc"))
  if(!is.null(range_res_mean) && is.finite(range_res_mean)) {
    p <- p + ggplot2::geom_vline(xintercept = range_res_mean, linetype = "dashed", colour = "#c0392b", linewidth = 0.5) +
      ggplot2::geom_text(
        data = data.frame(x = range_res_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Sre range"),
        ggplot2::aes(x = x, y = y, label = label),
        colour = "#c0392b", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }
  if(!is.null(range_lat_mean) && is.finite(range_lat_mean)) {
    p <- p + ggplot2::geom_vline(xintercept = range_lat_mean, linetype = "dashed", colour = "#2980b9", linewidth = 0.5) +
      ggplot2::geom_text(
        data = data.frame(x = range_lat_mean, y = max(sv_df$gamma, na.rm = TRUE), label = "Sshared range"),
        ggplot2::aes(x = x, y = y, label = label),
        colour = "#2980b9", angle = 90, hjust = 1, vjust = -0.5, size = 3)
  }
  p
}

