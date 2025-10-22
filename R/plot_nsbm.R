#' @name plot.nsbm.inlabru
#'
#' @title Plot for nsbm.inlabru...
#'
#' @description This function plots a single layer from a fitted NSBM model object ("nsbm.inlabru"). It can display the current prediction ("pred"), the spatial field ("pred_sp"), or any future/environmental scenario stored in "new.projections".
#'
#' @param x An object of class "nsbm.inlabru".
#' @param which Component to plot. One of: "pred", "pred_sp", "new.projections", "new.projections[[1]]", "new.projections[[NAME]]", or directly the scenario stored in x$new.projections.
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
#' ## Spatial field
#' # plot(myPred.pure, which = "pred_sp", layer = "sd")
#'
#' ## Scenario by index or name
#' # plot(myPred.pure, which = "new.projections[[1]]")
#' # plot(myPred.pure, which = "scenario1", layer = "q0.975")
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
    if (!is.null(x$current.projections$pred)) {
      r <- x$current.projections$pred
    } else if (!is.null(x$current.projections$pred.multiply)) {
      r <- x$current.projections$pred.multiply
    } else {
      r <- x$current.projections$pred.covariate
    }
    scope_label <- "Current"

  } else if (identical(which, "pred_sp")) {
    # spatial field
    r <- x$current.projections$pred_sp
    scope_label <- "Spatial field"

  } else {
    # new.projections
    np <- x$new.projections
    if (is.null(np) || !length(np)) {
      stop("No 'new.projections' in the object.")
    }

    get_np_by <- function(id) {
      nms <- names(np)
      if (is.numeric(id)) {
        id <- as.integer(id)
        if (id < 1 || id > length(np)) {
          stop("Index out of range in new.projections.")
        }
        list(obj = np[[id]], name = nms[id])
      } else {
        if (!(id %in% nms)) {
          stop(sprintf("Scenario '%s' not found in new.projections.", id))
        }
        list(obj = np[[id]], name = id)
      }
    }

    # first, [[index]] / [[name]], or direct name
    if (identical(which, "new.projections")) {
      pick <- get_np_by(1)
    } else if (grepl("^new\\.projections\\[\\[.*\\]\\]$", which)) {
      inside <- sub("^new\\.projections\\[\\[(.*)\\]\\]$", "\\1", which)
      if (grepl("^[0-9]+$", inside)) inside <- as.integer(inside)
      pick <- get_np_by(inside)
    } else if (which %in% names(np)) {
      pick <- get_np_by(which)
    } else {
      stop("`which` must be 'pred', 'pred_sp', 'new.projections', 'new.projections[[...]]', or an exact scenario name in new.projections.")
    }

    r <- pick$obj
    scope_label <- pick$name
  }

  rr <- if(inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r
  stopifnot(inherits(rr, "SpatRaster"))

  if (!layer %in% names(rr)) {
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

  auto_legend <- if (grepl("binomial", fam)) {
    if (identical(scope_label, "Spatial field")) legend_eta(ln) else legend_prob(ln)
  } else if (grepl("poisson", fam)) {
    if (identical(scope_label, "Spatial field")) legend_eta(ln) else legend_int(ln)
  } else if (grepl("gaussian", fam)) {
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
  lt <- wrap_legend_title(if (is.null(legend_title)) auto_legend else legend_title)

  # plot title 
  species   <- x$Species.Name
  type_lab  <- sub(".*NSBM\\.|\\s*\\(.*\\)", "", x$Summary$Value[grep("Model type", x$Summary$Field)][1])
  plot_title <- if (is.null(title)) {
    paste(species, "|", type_lab, "|", scope_label, "|", layer)
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

