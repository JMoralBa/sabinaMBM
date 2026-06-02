#' @title Summary for jmbm.inlabru Objects
#' @description Provides a summary for objects of class \code{jmbm.inlabru}.
#' @param object An object of class \code{jmbm.inlabru}.
#' @param ... Additional arguments (not used).
#' @return A summary in data frame format.
#' @seealso \code{\link{NSBM.inlabru}}
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


