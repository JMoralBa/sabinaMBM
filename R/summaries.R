#' @title Summary for nsbm.inlabru Objects
#' @description Provides a summary for objects of class \code{nsbm.inlabru}.
#' @param object An object of class \code{nsbm.inlabru}.
#' @param ... Additional arguments (not used).
#' @return A summary in data frame format.
#' @seealso \code{\link{NSBM.inlabru}}
#' @export
#' @method summary nsbm.inlabru
summary.nsbm.inlabru <- function(object, ...) {

  S <- object$Summary

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(" Summary for:", object$Species.Name, "\n")
  cat(strrep("=", 70), "\n\n", sep = "")

  print_block <- function(title, df) {
    if (is.null(df) || !nrow(df)) return()
    cat(title, "\n")
    print(df, row.names = FALSE)
    cat("\n")
  }

  print_block("------ Model metadata ------", S$Metadata)
  print_block("------ Model fit (Bayesian criteria) ------", S$`Model fit`)
  print_block("------ Hyperparameters ------", S$Hyperparameters)
  print_block("------ Intercepts ------", S$Intercepts)
  print_block("------ Fixed effects ------", S$`Fixed effects`)
  print_block("------ Predictive performance ------", S$`Predictive performance`)
  print_block("------ Calibration & coverage ------", S$`Calibration & coverage`)
  print_block("------ Diagnostics ------", S$Diagnostics)
  #print_block("------ Variance Decomposition ------", S$`Variance Decomposition`)

  if (!is.null(S$`Cross-validation`)) {
    print_block("------ Cross-validation ------", S$`Cross-validation`)
  }

  invisible(S)
}


# -----------------------------


#' prepare summary
#' @noRd
.nsbm_generate_summary <- function(fit, species, fam, lnk, coupling.intercept, coupling.predictors, diag_block, cv_res=NULL, vg=NULL, vr=NULL) {
  
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

  Sloc_local <- "Sloc" %in% names(fit$summary.random)
  Sshared_global <- "Sshared" %in% names(fit$summary.random)
  has_covariates <- nrow(fit$summary.fixed) > 0

  # metadata
  species_name <- gsub("\\.", " ", species)
  model_type <- paste0("NSBM",
   if(Sloc_local) " + Sloc" else "",
   if(Sshared_global) " + Sshared" else "",
   if(has_covariates) " + covariates" else ""
  )

  cp_int <- if(is.null(coupling.intercept)) "NULL" else coupling.intercept
  cp_pred <- if(is.list(coupling.predictors)) "custom list" else if(is.null(coupling.predictors)) "NULL" else coupling.predictors
  
  tbl_metadata <- data.frame(
    Field = c("Species name:", "Model type:", "Family | Link:", 
              "Coupling (Intercept):", "Coupling (Predictors):"),
    Value = c(species_name,
              model_type,
              paste0(fam, " | ", ifelse(is.null(lnk), "—", lnk)),
              cp_int,
              cp_pred),
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

  if(Sloc_local) {
    params <- c(
      params,
      "Sloc field: Range (posterior mean ± SD)",
      "Sloc field: Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_res,
      hyp_block$sigma_res
    )
  }
  if(Sshared_global) {
    params <- c(
      params,
      "Sshared field: Range (posterior mean ± SD)",
      "Sshared field: Sigma (posterior mean ± SD)"
    )
    values <- c(
      values,
      hyp_block$range_Sshared,
      hyp_block$sigma_Sshared
    )
  }
  params <- c(params, "IGlobal: Precision (mean ± SD)")
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
    if(!is.null(coupling.intercept) && coupling.intercept == "ordered_hierarchical") {
      rows[["IRegional (copy from IGlobal, mean ± SD, CI95%)"]] <- int_block$iRegional
    } else {
      rows[["IRegional (mean ± SD, CI95%)"]] <- int_block$iRegional
    }
  }
  # copy_beta only exists  coupling.intercept = "hierarchical
  if(!is.null(int_block$copy_beta)) {
    rows[["IRegional: Beta (copy) (mean ± SD)"]] <- int_block$copy_beta
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
  } else {
    # var name
    base_vars <- gsub("GL$|GL_ls$|RE$|RE_ss$|RE_oh$", "", tbl_fixed$coef)
    
    # which coupling?
    tbl_fixed$Coupling <- vapply(base_vars, function(v) {
      if(v %in% c(vg, vr)) {
        .resolve_coupling_predictor(v, coupling.predictors, vg)
      } else {
        "—"
      }
    }, character(1))
    
    cols_order <- c("coef", "coupling", "estimate", "sd", "2.5%", "97.5%", "signif", "tail_p")
    cols_order <- intersect(cols_order, names(tbl_fixed))
    tbl_fixed <- tbl_fixed[, cols_order]
    attr(tbl_fixed, "note") <- "Coefficients back-transformed to original covariate scale (standardized internally for numerical stability)."
  }

  # predictive performance
  tbl_pred <- data.frame(
    Metric = c("AUC (full model)",
               "Tjur R\u00b2 (discrimination coefficient)",
               "Brier score",
               "RMSE",
               "Observed-predicted correlation (r)"),
    Value  = c(diag_block$predictive$auc_full,
               fmt_val(diag_block$predictive$tjur_r2),
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
      "Residual Moran's I (obs - fitted mean)",
      "Max CI/median ratio (hyperparameters)",
      "Scale-separation ratio (range_Sshared / range_Sloc)",
      "Variance ratio (sigma_Sshared / sigma_Sloc)",
      "Sshared–Sloc field correlation (r)"
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
    metric_label <- paste0("CV ", cv_res$metric_name, " (mean ± sd)")
    tbl_cv <- data.frame(
      Metric = c("CV folds", metric_label),
      Value = c(cv_res$cv.folds,
                sprintf("%.3f ± %.3f", cv_res$metric_mean, cv_res$metric_sd)),
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
