.new_prism_result <- function(x) {
  structure(x, class = "prism_result")
}

#' Print a PRISM result
#' @param x A `prism_result`.
#' @param ... Unused.
#' @export
print.prism_result <- function(x, ...) {
  cat("PRISM covariance analysis\n")
  cat("  features: ", x$diagnostics$data$D, "\n", sep = "")
  cat("  samples: ", x$diagnostics$data$N, "\n", sep = "")
  cat("  composition estimator: ", x$parameters$composition_estimator, "\n", sep = "")
  cat("  bootstrap: ", x$parameters$bootstrap, "\n", sep = "")
  cat("  draws: ", x$parameters$S_effective, "\n", sep = "")
  invisible(x)
}
