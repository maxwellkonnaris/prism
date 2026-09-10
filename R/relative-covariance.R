#' Closure and pseudocount for compositional counts
#'
#' Adds a pseudocount and renormalizes each sample (column) to sum to one.
#' A zero pseudocount is only accepted when no count is exactly zero, since
#' a zero proportion has no logarithm.
#'
#' @param counts Non-negative features-by-samples count matrix.
#' @param pseudocount Non-negative scalar added to counts before closure.
#' @return A finite, strictly positive features-by-samples proportion matrix
#'   whose columns sum to one.
#' @keywords internal
prism_closure <- function(counts, pseudocount = 0.5) {
  Y <- as.matrix(counts)
  if (!is.numeric(Y) || length(Y) == 0L) {
    stop("counts must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (nrow(Y) < 2L || ncol(Y) < 2L) {
    stop("counts must have at least 2 features (rows) and 2 samples (columns).", call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(Y < 0)) {
    stop("counts must be finite and non-negative.", call. = FALSE)
  }
  pseudocount <- .assert_scalar_finite(pseudocount, "pseudocount", lower = 0)
  if (pseudocount == 0 && any(Y == 0)) {
    stop(
      "pseudocount = 0 is invalid when counts contains zeros: ",
      "a zero proportion has no logarithm.",
      call. = FALSE
    )
  }
  Yp <- Y + pseudocount
  totals <- colSums(Yp)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Each sample must have a positive, finite total after adding the pseudocount.", call. = FALSE)
  }
  sweep(Yp, 2, totals, "/")
}

#' Log relative composition
#'
#' @inheritParams prism_closure
#' @return A finite features-by-samples matrix of log proportions.
#' @keywords internal
prism_log_composition <- function(counts, pseudocount = 0.5) {
  log(prism_closure(counts, pseudocount))
}

#' Relative log-covariance
#'
#' Sample covariance of the log-composition across samples, i.e.
#' \code{Sigma_rel = Cov(log P)} for features-by-samples log proportions
#' \code{log P}.
#'
#' @param log_props A finite features-by-samples matrix of log proportions.
#' @return A features-by-features covariance matrix.
#' @keywords internal
prism_relative_covariance <- function(log_props) {
  X <- as.matrix(log_props)
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("log_props must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("At least 2 samples are required to estimate a covariance.", call. = FALSE)
  }
  stats::cov(t(X))
}
