#' @name plot.nsbm.inlabru
#'
#' @title Plot predictions and spatial fields from a fitted NSBM model
#'
#' @description Prints a structured summary for objects of class \code{nsbm.inlabru}, including model metadata, Bayesian fit criteria, hyperparameters, intercepts, fixed effects, predictive performance, calibration, and diagnostics. Returns the summary list invisibly.
#'
#' @param x An object of class "nsbm.inlabru".
#' @param which Component to plot. One of:
#'   \itemize{
#'     \item \code{"pred"}: current suitability prediction (default).
#'     \item \code{"pred_local"}: Sloc field (requires SPDE).
#'     \item \code{"pred_shared"}: Sshared field (requires SPDE).
#'     \item \code{"hyperparams"}: posterior marginals of SPDE hyperparameters with PC-priors overlaid.
#'     \item \code{"intercepts"}: posterior marginals of IGlobal and IRegional intercepts.
#'     \item \code{"fixed"}: posterior marginals of all fixed effect coefficients, back-transformed to original covariate scale.
#'     \item \code{"pit"}: histogram of Probability Integral Transform values. A uniform distribution indicates good calibration.
#'     \item \code{"new.projections"}, \code{"new.projections[[1]]"}, or scenario name: future/alternative scenario.
#'   }
#' @param layer Raster layer to display (default = "mean"). Available: "mean", "sd", "q0.025", "q0.5", "q0.975", "median", "sd.mc_std_err", "mean.mc_std_err".
#' @param palette Color palette passed to ggplot2::scale_fill_distiller(). Default = "Spectral".
#' @param title Optional plot title. If NULL, a default title is generated as "Species - type - scope - layer".
#' @param legend_title Optional legend title. If NULL, it is automatically inferred from the model family and the layer.
#'
#' @return A ggplot object showing the selected NSBM output layer.
#'
#' @examples
#' # Default: current prediction (mean layer)
#' plot(myPred.pure)
#'
#' ## Sloc field
#' # plot(myPred.pure, which = "pred_local", layer = "sd")
#'
#' ## Scenario by index or name
#' # plot(myPred.pure, which = "new.projections[[1]]")
#' # plot(myPred.pure, which = "scenario1", layer = "q0.975")
#'
#' @seealso \code{\link{NSBM.pure}}, \code{\link{summary.nsbm.inlabru}}
#'
#' @export
#' @method plot nsbm.inlabru
plot.nsbm.inlabru <- function(x,
                              which = "pred",
                              layer = "mean",
                              palette = "Spectral",
                              title = NULL,
                              legend_title = NULL) {

  stopifnot(inherits(x, "nsbm.inlabru"))

  scope_label <- NULL
  r <- NULL

  if (identical(which, "pred")) {
    # current
    r <- x$current.projections$pred
    scope_label <- "Current"

  } else if(identical(which, "pred_local")) {
    r <- x$current.projections$pred_local
    if(is.null(r)) stop("❌ No Sloc field ('pred_local') in this model. Refit with `local.pcprior.range` and `local.pcprior.sigma`.\n")
    scope_label <- "Sloc field"

  } else if(identical(which, "pred_shared")) {
    r <- x$current.projections$pred_shared
    if(is.null(r)) stop("❌ No Sshared field ('pred_shared') in this model. Refit with `shared.pcprior.range` and `shared.pcprior.sigma`.\n")
    scope_label <- "Sshared field"

  } else if(identical(which, "hyperparams")) {
    marg <- x$marginals$hyperpar
    if(is.null(marg) || length(marg) == 0)
      stop("❌ No hyperparameter marginals in this model. Refit with SPDE priors.\n")
    species <- gsub("\\.", " ", x$Species.Name)
    param_labels <- c(
      "Range for Sloc" = "Sloc field: Range",
      "Stdev for Sloc" = "Sloc field: Sigma",
      "Range for GLspde" = "Sshared field: Range",
      "Stdev for GLspde" = "Sshared field: Sigma",
      "Precision for IGlobal" = "IGlobal: Precision",
      "Precision for IRegional" = "IRegional: Precision",
      "Beta for IRegional" = "IRegional: Beta (copy)"
    )
    df_list <- lapply(names(marg), function(nm) {
      m  <- INLA::inla.smarginal(marg[[nm]])
      label <- if(nm %in% names(param_labels)) param_labels[[nm]] else nm
      data.frame(x = m$x, y = m$y, param = label)
    })
    df <- do.call(rbind, df_list)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(colour = "#2E75B6", linewidth = 0.8) +
      ggplot2::facet_wrap(~ param, scales = "free") +
      ggplot2::labs(
        title = if(!is.null(title)) title else paste(species, "| Hyperparameter posteriors"),
        x = "Value", y = "Density") +
      ggplot2::theme_minimal()
    return(p)

  } else if(identical(which, "intercepts")) {
    marg_r <- x$marginals$random
    if(is.null(marg_r) || length(marg_r) == 0)
      stop("❌ No intercept marginals in this model.\n")
    species <- gsub("\\.", " ", x$Species.Name)
    coupling_label <- if(!is.null(x$args$coupling.intercept)) x$args$coupling.intercept else "regional only"
    int_names <- intersect(c("IGlobal", "IRegional"), names(marg_r))
    if(length(int_names) == 0)
      stop("❌ IGlobal/IRegional not found in marginals$random.\n")
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
        legend.position      = c(0.88, 0.88),
        legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
        legend.key.size      = ggplot2::unit(0.4, "cm"),
        legend.text          = ggplot2::element_text(size = 9))

    return(p)

  } else if(identical(which, "fixed")) {
    marg_f <- x$marginals$fixed
    sp     <- x$scale_params
    if(is.null(marg_f) || length(marg_f) == 0)
      stop("❌ No fixed effect marginals in this model.\n")
    species <- gsub("\\.", " ", x$Species.Name)
    cp_label <- if(!is.null(x$args$coupling.predictors)) x$args$coupling.predictors else "NULL"
    # back-transform marginals to original covariate scale
    df_list <- lapply(names(marg_f), function(nm) {
      base_var <- gsub("GL$|RE$|RE_oh$|GL_ls$|RE_ss$", "", nm)
      m <- marg_f[[nm]]
      if(!is.null(sp) && base_var %in% names(sp) && sp[[base_var]]$sd > 0) {
        sd_x <- sp[[base_var]]$sd
        m <- INLA::inla.tmarginal(function(x) x / sd_x, m)
      }
      sm <- INLA::inla.smarginal(m)
      data.frame(x = sm$x, y = sm$y, coef = nm)
    })
    df <- do.call(rbind, df_list)
    # significance: CI do not cross 0
    signif_coefs <- names(marg_f)[vapply(names(marg_f), function(nm) {
      base_var <- gsub("GL$|RE$|RE_oh$|GL_ls$|RE_ss$", "", nm)
      m <- marg_f[[nm]]
      if(!is.null(sp) && base_var %in% names(sp) && sp[[base_var]]$sd > 0)
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
        subtitle = paste("coupling.predictors =", cp_label, "| Grey = non-significant (95% CI)"),
        x = "Coefficient value", y = "Density") +
      ggplot2::theme_minimal()

    # add "NS" annotation to non-significant panels
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
    p <- x$diagnostic_plots$correlogram
    if(is.null(p))
      stop("❌ No correlogram available. Refit with Sloc or Sshared SPDE.\n")
    if(!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)

  } else if(identical(which, "semivariogram")) {
    p <- x$diagnostic_plots$semivariogram
    if(is.null(p))
      stop("❌ No semivariogram available. Refit with Sloc or Sshared SPDE.\n")
    if(!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)

  } else if(identical(which, "pit")) {
    pit <- x$pit_values
    if(is.null(pit) || length(pit) == 0)
      stop("❌ No PIT values in this model. Refit with control.compute = list(cpo = TRUE).\n")
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
      stop("No 'new.projections' in the object.")
    }

    get_np_by <- function(id) {
      nms <- names(np)
      if(is.numeric(id)) {
        id <- as.integer(id)
        if(id < 1 || id > length(np)) {
          stop("Index out of range in new.projections.")
        }
        list(obj = np[[id]], name = nms[id])
      } else {
        if(!(id %in% nms)) {
          stop(sprintf("Scenario '%s' not found in new.projections.", id))
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
      stop("`which` must be 'pred', 'pred_local', 'new.projections', 'new.projections[[...]]', or an exact scenario name in new.projections.")
    }

    r <- pick$obj
    scope_label <- pick$name
  }

  rr <- if(inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r
  stopifnot(inherits(rr, "SpatRaster"))

  if(!layer %in% names(rr)) {
    stop(paste0("Layer '", layer, "' does not exist. Available layers: ",
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

  auto_legend <- if(scope_label %in% c("Sloc field", "Sshared field")) {
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
    legend_prob(ln) #@@@JMB arreglar para v2 con nbinomial, beta, tweedie...
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

