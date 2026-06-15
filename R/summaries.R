#' @title Summary for jmbm.inlabru objects
#' @description Prints a structured summary for objects of class \code{jmbm.inlabru}, including model metadata, Bayesian fit criteria, hyperparameters, intercepts, fixed effects, predictive performance, and spatial diagnostics.
#' @param object An object of class \code{jmbm.inlabru}.
#' @param ... Additional arguments (not used).
#' @return A summary in data frame format.
#' @seealso \code{\link{MBM.Modelling}}, \code{\link{plot.jmbm.inlabru}}
#' @export
#' @method summary jmbm.inlabru
summary.jmbm.inlabru <- function(object, ...) {

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

  fe <- S$`Fixed effects`
  if (!is.null(fe) && nrow(fe) > 0 && "coef" %in% names(fe)) {
    cat("------ Fixed effects ------\n")
    is_regional <- grepl("RE($|_oh$|_reg_anom$)", fe$coef)
    fe_gl <- fe[!is_regional, , drop = FALSE]
    fe_re <- fe[ is_regional, , drop = FALSE]
    if (nrow(fe_gl) > 0) {
      cat("  [Global scale]\n")
      print(fe_gl, row.names = FALSE)
    }
    if (nrow(fe_gl) > 0 && nrow(fe_re) > 0) cat("\n")
    if (nrow(fe_re) > 0) {
      cat("  [Regional scale]\n")
      print(fe_re, row.names = FALSE)
    }
    # only shown when non-trivial suffixes are present
    has_oh <- any(grepl("_oh$",       fe$coef))
    has_glo_res <- any(grepl("_glo_res$",  fe$coef))
    has_reg_ano <- any(grepl("_reg_anom$", fe$coef))
    cat("\n")
    if (has_oh || has_glo_res || has_reg_ano) {
      cat("  Note: GL = global scale; RE = regional scale;",
          "_oh = ordered-hierarchical constraint (beta_RE ~ N(1, 0.5^2), native Z-scale);",
          "_glo_res = large-scale macro-trend (scale-decomposed, unified Z-scale);",
          "_reg_anom = fine-scale anomaly (scale-decomposed, unified Z-scale).",
          "\n  Coefficients are back-transformed to original covariate units (effect per unit of X).",
          "\n  Variables with scale_decomposed or bayesian_feedback coupling use global sigma (sigma_GL)",
          "as the standardization denominator; all other variables use their native sigma.\n")
    } else {
      cat("  Note: GL = global scale; RE = regional scale.",
          "\n  Coefficients are back-transformed to original covariate units (effect per unit of X).\n")
    }
    cat("\n")
  } else {
    print_block("------ Fixed effects ------", fe)
  }

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


