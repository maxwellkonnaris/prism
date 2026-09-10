.new_prism_result <- function(x) {
  structure(x, class = c("prism_result", "list"))
}

#' Print a PRISM result
#'
#' @param x A `prism_result` object.
#' @param ... Unused.
#' @export
print.prism_result <- function(x, ...) {
  dimensions <- x$diagnostics$data[c("D", "N")]
  cat("PRISM covariance analysis\n")
  if (length(dimensions) == 2L) {
    cat("  features: ", dimensions[["D"]], "\n", sep = "")
    cat("  samples: ", dimensions[["N"]], "\n", sep = "")
  }
  cat("  composition estimator: ", x$parameters$composition_estimator, "\n", sep = "")
  cat("  Monte Carlo draws: ", x$parameters$S, "\n", sep = "")
  cat("  streamed: ", if (isTRUE(x$parameters$stream_active)) "yes" else "no", "\n", sep = "")
  invisible(x)
}

.finalize_prism_result <- function(x, estimator) {
  diagnostics <- list(
    data = x$data_diagnostics,
    scale = x$params$scale_summary,
    composition = list(
      name = estimator$name,
      fit_summary = x$mln_fit_summary,
      model = x$mln_diagnostics
    ),
    sampling = x$sampling_diagnostics,
    high_resolution = x$high_resolution_diagnostics,
    simulation = list(
      identification_width_median = x$id_region_width_median,
      lower_endpoint_iqr = x$lower_iqr,
      upper_endpoint_iqr = x$upper_iqr
    )
  )
  x$diagnostics <- diagnostics
  x$parameters <- x$params
  x$parameters$composition_estimator <- estimator$name
  x$params <- x$parameters
  x$data_diagnostics <- NULL
  x$mln_fit_summary <- NULL
  x$mln_diagnostics <- NULL
  x$sampling_diagnostics <- NULL
  x$high_resolution_diagnostics <- NULL
  x$id_region_width_median <- NULL
  x$lower_iqr <- NULL
  x$upper_iqr <- NULL
  .new_prism_result(x)
}
