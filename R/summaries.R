#' @title Summary for nsbm.inlabru Objects
#' @description Provides a summary for objects of class \code{nsbm.inlabru}.
#' @param object An object of class \code{nsbm.inlabru}.
#' @param ... Additional arguments (not used).
#' @return A summary in data frame format.
#' @seealso \code{\link{NSBM.inlabru}}
#' @export
#' @method summary nsbm.inlabru
summary.nsbm.inlabru <- function(object, ...) {
  summary_df <- as.data.frame(object$Summary)
  return(summary_df)
}
