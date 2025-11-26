#' @title Summary for nsbm.inlabru Objects
#' @description Provides a summary for objects of class \code{nsbm.inlabru}.
#' @param object An object of class \code{nsbm.inlabru}.
#' @param ... Additional arguments (not used).
#' @return A summary in data frame format.
#' @seealso \code{\link{NSBM.inlabru}}
#' @export
#' @method summary nsbm.inlabru
#summary.nsbm.inlabru <- function(object, ...) {
#  summary_df <- as.data.frame(object$Summary)
#  return(summary_df)
#}

summary.nsbm.inlabru <- function(object, ...) {

  S <- object$Summary
  
  if (!is.list(S)) {
    # compatibilitat amb versions antigues
    print(as.data.frame(S))
    return(invisible(S))
  }

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(" Summary for:", object$Species.Name, "\n")
  cat(strrep("=", 70), "\n\n", sep = "")

  print_block <- function(title, df) {
    if (is.null(df) || !nrow(df)) return()
    cat(title, "\n")
    print(df, row.names = FALSE)
    cat("\n")
  }

  # ordre dels blocs basat en nsbm_generate_summary()
  print_block("------ Model metadata ------", S$Metadata)
  print_block("------ Model fit (Bayesian criteria) ------", S$`Model fit`)
  print_block("------ Hyperparameters ------", S$Hyperparameters)
  print_block("------ Intercepts ------", S$Intercepts)
  print_block("------ Fixed effects ------", S$`Fixed effects`)
  print_block("------ Predictive performance ------", S$`Predictive performance`)
  print_block("------ Calibration & coverage ------", S$`Calibration & coverage`)
  print_block("------ Diagnostics ------", S$Diagnostics)
  print_block("------ Variance Decomposition ------", S$`Variance Decomposition`)

  if (!is.null(S$`Cross-validation`)) {
    print_block("------ Cross-validation ------", S$`Cross-validation`)
  }

  invisible(S)
}