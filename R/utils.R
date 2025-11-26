#' prepare significant covariates
#' @noRd
signif_vars <- function(fit) {
  sf <- fit$summary.fixed
  
  # p-value aprox with marginal posterior
  tailprob <- function(term){
    if (!is.null(fit$marginals.fixed) && !is.null(fit$marginals.fixed[[term]])) {
      p0 <- INLA::inla.pmarginal(0, fit$marginals.fixed[[term]])
      2 * min(p0, 1 - p0)
    } else {
      # fallback Normal
      2 * stats::pnorm(-abs(sf[term,"mean"] / sf[term,"sd"]))
    }
  }

  # 95% CI no cruza 0
  signif95 <- sf[,"0.025quant"] * sf[,"0.975quant"] > 0
  
  stars <- ifelse(signif95, "***", "")
  
  out <- data.frame(
    coef = rownames(sf),
    estimate = sf$mean,
    sd = sf$sd,
    `2.5%` = sf[,"0.025quant"],
    `97.5%` = sf[,"0.975quant"],
    signif = stars,
    tail_p = vapply(rownames(sf), tailprob, numeric(1)),
    row.names = NULL,
    check.names = FALSE
  )
  
  out <- out[order(-abs(out$estimate)), ]

  out$estimate <- round(out$estimate, 6)
  out$sd <- round(out$sd, 6)
  out$`2.5%` <- round(out$`2.5%`, 6)
  out$`97.5%` <- round(out$`97.5%`, 6)
  
  return(out)
}

#' plot covariates importance
#' @noRd
nsbm_vars_importance <- function(fit) {

  df <- fit$summary.fixed
  df$var <- rownames(df)

  # exclude intercepts and latent
  drop <- c("IGlobal", "IRegional", "beta_GL")
  df <- df[!(df$var %in% drop), , drop = FALSE]

  df$sign <- ifelse(df$mean > 0, "Positive", "Negative")

  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = reorder(var, abs(mean)),
                 y = abs(mean),
                 fill = sign)) +    ggplot2::geom_col(alpha = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Positive" = "#1a5276","Negative" = "#c0392b")) +
    ggplot2::labs(
      title = "F) Variable importance (|β|)",
      subtitle = "Higher bars indicate stronger effect magnitude",
      x = "Covariate",
      y = "|Effect size|"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold",size = 11,color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 9,color = "#5d6d7e"),
      axis.title.x = ggplot2::element_text(size = 9,margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(size = 9,margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(size = 8,color = "#2c3e50"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2,color = "#d5d8dc")
    )


  return(p)
}


#' prepare summary
#' @noRd
nsbm_generate_summary <- function(fit, 
                                 species, 
                                 fam, lnk, 
                                 diag_block, 
                                 cv_res=NULL) {
  
  fmt_val <- function(x, digits = 3) {
    if(is.null(x) || length(x) == 0) return("—")
    if(is.character(x)) return(x)
    if(!is.finite(x)) return("—")
    round(x, digits)
  }

  fmt_pm <- function(mean, sd, digits = 3) {
    if(!is.finite(mean) || !is.finite(sd)) return(NULL)
    paste0(round(mean, digits), " ± ", round(sd, digits))
  }

  spatial_local <- "spatial" %in% names(fit$summary.random)
  latent_global <- "GLspde" %in% names(fit$summary.random)
  has_covariates <- nrow(fit$summary.fixed) > 0

  # metadata
  species_name <- gsub("\\.", " ", species)
  model_type <- paste0("NSBM",
   if(spatial_local) " + spatial" else "",
   if(latent_global) " + latent" else "",
   if(has_covariates) " + covariates" else ""
  )

  tbl_metadata <- data.frame(
    Field = c("Species name:", "Model type:", "Family | Link:"),
    Value = c(species_name,
              model_type,
              paste0(fam, " | ", ifelse(is.null(lnk), "—", lnk))),
    stringsAsFactors = FALSE
  )

  # bayesian model fit
  tbl_fit <- data.frame(
    Metric = c(
      "DIC",
      "WAIC",
      "Marginal log-likelihood (log ML)",
      "LCPO (sum log-CPO)",
      "MLPD (mean log predictive density)"
    ),
    Value = c(diag_block$bayes_fit$dic_val, 
              diag_block$bayes_fit$waic_val,
              diag_block$bayes_fit$mlik_val, 
              diag_block$bayes_fit$lcpo_val,
              diag_block$bayes_fit$mlpd_val),
    stringsAsFactors = FALSE
  )

  # hyperparametres
  hyp_block <- diag_block$hiperpars
  params <- character(0)
  values <- character(0)

  if(spatial_local) {
    params <- c(
      params,
      "Spatial local field – Range (posterior mean ± SD)",
      "Spatial local field – Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_res,
      hyp_block$sigma_res
    )
  }
  if(latent_global) {
    params <- c(
      params,
      "Latent global field – Range (posterior mean ± SD)",
      "Latent global field – Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_lat,
      hyp_block$sigma_lat
    )
  }
  params <- c(params, "Precision for IGlobal (mean ± SD)")
  values <- c(values, hyp_block$prec_IGlobal)

  tbl_hyper <- data.frame(
    Parameter = params,
    Value = values,
    stringsAsFactors = FALSE
  )

  # intercepts
  int_block <- diag_block$intercepts
  rows <- list()
  # IGlobal exists only for coupling.intercept = "additive"/"hierarchical"
  if(!is.null(int_block$iGlobal)) {
    rows[["IGlobal (mean ± SD, CI95%)"]] <- int_block$iGlobal
  }
  # IRegional always exists
  if(!is.null(int_block$iRegional)) {
    if(!is.null(int_block$copy_beta)) {
      rows[["IRegional (copy from IGlobal, mean ± SD, CI95%)"]] <- int_block$iRegional
    } else {
      rows[["IRegional (mean ± SD, CI95%)"]] <- int_block$iRegional
    }
  }
  # copy_beta only exists  coupling.intercept = "hierarchical
  if(!is.null(int_block$copy_beta)) {
    rows[["Copy β (IGlobal -> IRegional) (mean ± SD)"]] <- int_block$copy_beta
  }

  tbl_intercepts <- data.frame(
    Term = names(rows),
    Value = unname(unlist(rows)),
    stringsAsFactors = FALSE
  )


  # fixed covariates
  tbl_fixed <- diag_block$fixed_covariates
  if(is.null(tbl_fixed) || !nrow(tbl_fixed)) {
    tbl_fixed <- data.frame(Term = "No significant covariates detected")
  }

  # predictive performance
  tbl_pred <- data.frame(
    Metric = c("AUC (full model)", 
               "Brier score", 
               "RMSE", 
               "Observed-predicted correlation (r)"),
    Value  = c(diag_block$predictive$auc_full,
               fmt_val(diag_block$predictive$brier),
               fmt_val(diag_block$predictive$rmse),
               fmt_val(diag_block$predictive$corr_obs_pred)),
    stringsAsFactors = FALSE
  )

  # calibration & coverage
  cal_block <- diag_block$calibration
  cov_block <- diag_block$coverage

  tbl_cal <- data.frame(
    Metric = c("Calibration slope",
               "PIT KS p-value",
               "Coverage (central 50%)",
               "Coverage (central 95%)"),
    Value = c(fmt_val(cal_block$slope),
              fmt_val(cal_block$ks_pit),
              fmt_val(cov_block$cov50),
              fmt_val(cov_block$cov95)),
    stringsAsFactors = FALSE
  )

  # diagnostics
  tbl_diag <- data.frame(
    Metric = c(
      "Residual Moran's I",
      "Max CI/median ratio (hyperparameters)",
      "Scale-separation ratio (range_latent / range_spatial)",
      "Variance ratio (sigma_latent / sigma_spatial)",
      "Latent–residual field correlation (r)"
    ),
    Value = c(
      fmt_val(diag_block$diagnostics$moran_I),
      fmt_val(diag_block$diagnostics$max_CIratio),
      fmt_val(diag_block$diagnostics$range_ratio),
      fmt_val(diag_block$diagnostics$sigma_ratio),
      fmt_val(diag_block$diagnostics$field_correlation)
    ),
    stringsAsFactors = FALSE
  )


  # cv
  tbl_cv <- NULL
  if(!is.null(cv_res)) {
    tbl_cv <- data.frame(
      Metric = c("CV folds", "CV AUC (mean ± sd)"),
      Value = c(cv_res$cv.folds,
                sprintf("%.3f ± %.3f", cv_res$auc_mean, cv_res$auc_sd)),
      stringsAsFactors = FALSE
    )
  }


  drop_null_rows <- function(df) {
    ok <- !sapply(df$Value, function(x) is.null(x) || identical(x, "—"))
    df[ok, , drop = FALSE]
  }

  tbl_fit <- drop_null_rows(tbl_fit)
  tbl_hyper <- drop_null_rows(tbl_hyper)
  tbl_intercepts <- drop_null_rows(tbl_intercepts)
  tbl_pred <- drop_null_rows(tbl_pred)
  tbl_cal <- drop_null_rows(tbl_cal)
  tbl_diag <- drop_null_rows(tbl_diag)

  # list of tables
  out <- list(
    Metadata = tbl_metadata,
    `Model fit` = tbl_fit,
    Hyperparameters = tbl_hyper,
    Intercepts = tbl_intercepts,
    `Fixed effects` = tbl_fixed,
    `Predictive performance` = tbl_pred,
    `Calibration & coverage` = tbl_cal,
    Diagnostics = tbl_diag
  )

  if(!is.null(tbl_cv))
    out[["Cross-validation"]] <- tbl_cv

  return(out)
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

