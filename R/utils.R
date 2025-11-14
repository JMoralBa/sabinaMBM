#' prepare summary
#' @noRd
generate_summary_nsbm <- function(fit, species=species, lcpo_val) {
  
  spatial_local <- "spatial" %in% names(fit$summary.random)
  latent_global <- "GLspde" %in% names(fit$summary.random)

  # Fixed effects and hyperpar
  summary_fixed <- fit$summary.fixed
  hyper <- fit$summary.hyperpar

  # Species and model type
  species_name <- gsub("\\.", " ", species)
  model_type <- paste0("NSBM",
   if(spatial_local) " + spatial" else "",
   if(latent_global) " + latent" else ""
  )

  # Model fit metrics
  dic_val <- if(!is.null(fit$dic$dic) && !is.na(fit$dic$dic)) round(fit$dic$dic, 2) else "Not computed"
  waic_val <- if(!is.null(fit$waic$waic) && !is.na(fit$waic$waic)) round(fit$waic$waic, 2) else "Not computed"
  mlik_val <- if(!is.null(fit$mlik) && !is.na(fit$mlik[1,1])) round(fit$mlik[1, 1], 2) else "Not computed"
  

  # Hyperparameters spde
  range_res <- if(!is.null(hyper) && "Range for spatial" %in% rownames(hyper))
    paste0(round(hyper["Range for spatial","mean"],2), " ± ", round(hyper["Range for spatial","sd"],2)) else "—"
  sigma_res <- if(!is.null(hyper) && "Stdev for spatial" %in% rownames(hyper))
    paste0(round(hyper["Stdev for spatial","mean"],2), " ± ", round(hyper["Stdev for spatial","sd"],2)) else "—"
  range_lat <- if(!is.null(hyper) && "Range for GLspde" %in% rownames(hyper))
    paste0(round(hyper["Range for GLspde","mean"],2), " ± ", round(hyper["Range for GLspde","sd"],2)) else "—"
  sigma_lat <- if(!is.null(hyper) && "Stdev for GLspde" %in% rownames(hyper))
    paste0(round(hyper["Stdev for GLspde","mean"],2), " ± ", round(hyper["Stdev for GLspde","sd"],2)) else "—"


  # Significant vars
  valid_vars <- summary_fixed[!is.na(summary_fixed$mean), ]
  valid_vars$type <- ifelse(
    rownames(valid_vars) %in% c("IGlobal", "IRegional"), "Intercept",
    ifelse(grepl("GL$", rownames(valid_vars)), "Global",
           ifelse(grepl("RE$", rownames(valid_vars)), "Regional", "Unclassified"))
  )
  # Filter significant vars (CI does not cross 0 and mean is relevant)
  signif_vars <- valid_vars[
    valid_vars[,"0.025quant"] * valid_vars[,"0.975quant"] > 0 &
    abs(valid_vars[,"mean"]) > 0.05,
  ]
  # Order by abs mean (importance?)
  signif_vars <- signif_vars[order(-abs(signif_vars[,"mean"])), ]
  var_labels <- paste0(
    rownames(signif_vars), " ",
    ifelse(signif_vars$mean > 0, "(+)", "(–)")
  )

  #
  summary_df <- data.frame(
    Field = c(
      "Species name:",
      "Model type:", 
      "",
      "-- Model fit (Bayesian criteria) --",
      "  DIC",
      "  WAIC",
      "  Marginal log-likelihood (log ML)",
      # log sum of conditional predictive ordinates: a bayesian metric for model validation (lower values indicate better fit)
      "  LCPO (sum of log-CPO values)",
      "",
      "--------- Model parameters --------",   
      if(spatial_local) "  Spatial/residual local field:",
      if(spatial_local) "      Range (posterior mean ± SD)" else NULL,
      if(spatial_local) "      Standard deviation (σ) (posterior mean ± SD)" else NULL,
      "",
      if(latent_global) "  Latent global field:" else NULL,
      if(latent_global) "      Range (posterior mean ± SD)" else NULL,
      if(latent_global) "      Standard deviation (σ) (posterior mean ± SD)" else NULL,
      "",
      "  Significant covariates (ordered):"
    ),
    Value = c(
      species_name, model_type,
      "", "", 
      dic_val, waic_val, mlik_val, lcpo_val,
      "", "",
      if(spatial_local) "" else NULL,
      if(spatial_local) range_res else NULL,
      if(spatial_local) sigma_res else NULL,
      "",
      if(latent_global) "" else NULL,
      if(latent_global) range_lat else NULL,
      if(latent_global) sigma_lat else NULL,
      "",
      if(length(var_labels) > 0) paste(var_labels, collapse = ", ") else "None"
    ),
    stringsAsFactors = FALSE
  )

  return(summary_df)
}


###----------------###

#' from sf to tif
#' @noRd
pred_as_tif <- function(pred, template, vars_to_export = c("mean", 
                                                           "sd", 
                                                           "q0.025", 
                                                           "q0.5", 
                                                           "q0.975",
                                                           "median",
                                                           "sd.mc_std_err",
                                                           "mean.mc_std_err")) {

  stopifnot(inherits(pred, c("bru_prediction", "sf")))
  stopifnot(inherits(template, "SpatRaster"))

  r_stack <- lapply(vars_to_export, function(var) {
    if(!is.null(pred[[var]])) {
      df <- as.data.frame(cbind(sf::st_coordinates(pred), value = pred[[var]]))
      sf_pts <- sf::st_as_sf(df, coords = c("X", "Y"), crs = sf::st_crs(pred))
      sf_pts <- sf::st_transform(sf_pts, terra::crs(template))
      r <- terra::rasterize(sf_pts, template, field = "value")
      names(r) <- var
      return(r)
    } else {
      return(NULL)
    }
  })

  r_stack <- Filter(function(x) inherits(x, "SpatRaster"), r_stack)
  if(length(r_stack) == 0) return(NULL)

  r_pred <- do.call(c, r_stack)
  return(r_pred)
}

###----------------###

#' random distribution of data in k-folds
#' @noRd
make_stratified_kfolds <- function(y, K) {
  # y: vector de 1/0
  # K: número de folds
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  f <- integer(length(y))
  f[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  f[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  f
}
