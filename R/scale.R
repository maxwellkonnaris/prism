#' Validate raw sample-level scale input
#'
#' `scale` is either one raw measurement per sample (a length-N vector) or
#' a matrix of replicate measurements per sample (N rows, one column per
#' replicate). Both are normalized to an N x R positive matrix so the rest
#' of the package only has to handle one shape.
#'
#' @param scale Positive numeric vector (length N) or matrix (N rows).
#' @param N Expected number of samples.
#' @return A finite, strictly positive N x R matrix.
#' @keywords internal
.validate_scale_input <- function(scale, N) {
  if (is.null(dim(scale))) {
    if (length(scale) != N) {
      stop("scale must have length N (one value per sample) or N rows.", call. = FALSE)
    }
    scale <- matrix(as.numeric(scale), nrow = N, ncol = 1L)
  } else {
    scale <- as.matrix(scale)
    if (nrow(scale) != N) {
      stop("scale must have length N (one value per sample) or N rows.", call. = FALSE)
    }
  }
  if (!is.numeric(scale) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("scale must be finite and strictly positive.", call. = FALSE)
  }
  scale
}

#' Select one raw scale value per sample and log it
#'
#' For samples with a single measurement this is just `log(scale)`. For
#' samples with replicate measurements, `replicate_index` (one draw per
#' sample, produced by the caller so all RNG state lives in the bootstrap
#' loop) selects which replicate is used, propagating measurement
#' uncertainty across bootstrap draws by resampling replicates the same
#' way samples are resampled.
#'
#' @param scale_mat N x R matrix from `.validate_scale_input()`.
#' @param sample_index Integer sample indices for this draw (length N,
#'   with replacement under the bootstrap).
#' @param replicate_index Integer replicate indices, one per entry of
#'   `sample_index`, in `1:ncol(scale_mat)`. Ignored when `scale_mat` has
#'   one column.
#' @return A numeric vector of log-scale values, length `length(sample_index)`.
#' @keywords internal
.select_scale_log <- function(scale_mat, sample_index, replicate_index = NULL) {
  R <- ncol(scale_mat)
  ri <- if (R == 1L) rep(1L, length(sample_index)) else replicate_index
  log(scale_mat[cbind(sample_index, ri)])
}

#' Estimate scale SD and taxon-scale correlation bounds from scale_log
#'
#' Fits `sigma = SD(log scale)` and, for each feature `d`,
#' `rho_d = Cor(log P_d, log scale)` from one sample of paired
#' (log-composition, log-scale) data, with an optional confidence interval
#' for each. The point estimate is always returned as `rho_witness`; when
#' an interval is requested, each feature's interval is a chi-square
#' interval for sigma and a Fisher-z interval for rho, both widened, never
#' narrowed, to contain the point estimate.
#'
#' @param log_proportions Features-by-samples log-composition matrix.
#' @param scale_log Sample-aligned finite log-scale vector.
#' @param ci_level Interval coverage in `[0, 1)`. `0` returns point
#'   estimates for both endpoints.
#' @param estimate_rho If `FALSE`, `rho_L`/`rho_U` are `NULL` and only
#'   sigma is estimated.
#' @param lower_zero If `TRUE`, force the lower sigma endpoint to 0
#'   instead of the chi-square lower quantile. Sigma cannot be negative,
#'   so the ordinary two-sided lower quantile asserts a positive floor
#'   that is an artifact of splitting the interval evenly across both
#'   tails, not a real constraint.
#' @return A list with `sigma_L`, `sigma_U`, `rho_L`, `rho_U`, `rho_witness`.
#' @keywords internal
.estimate_scale_log_bounds <- function(log_proportions, scale_log, ci_level = 0,
                                        estimate_rho = TRUE, lower_zero = FALSE) {
  X <- as.matrix(log_proportions)
  u <- as.numeric(scale_log)
  n <- length(u)
  if (ncol(X) != n || n < 2L || any(!is.finite(X)) || any(!is.finite(u))) {
    stop("log_proportions and scale_log must be finite and sample-aligned.", call. = FALSE)
  }
  ci_level <- .assert_scalar_finite(ci_level, "ci_level", lower = 0, upper = 1, upper_inclusive = FALSE)

  sigma_hat <- max(stats::sd(u), 0)
  rho_hat <- if (estimate_rho) .feature_correlations(X, u) else NULL

  if (ci_level == 0) {
    return(list(
      sigma_L = sigma_hat, sigma_U = sigma_hat,
      rho_L = rho_hat, rho_U = rho_hat, rho_witness = rho_hat
    ))
  }

  alpha <- 1 - ci_level
  df <- n - 1L
  sigma_L <- if (lower_zero) 0 else sqrt(df * sigma_hat^2 / stats::qchisq(1 - alpha / 2, df = df))
  sigma_U <- sqrt(df * sigma_hat^2 / stats::qchisq(alpha / 2, df = df))

  rho_L <- rho_U <- NULL
  if (estimate_rho) {
    if (n < 4L) {
      rho_L <- rho_U <- rho_hat
    } else {
      eps <- sqrt(.Machine$double.eps)
      z <- atanh(pmax(-1 + eps, pmin(1 - eps, rho_hat)))
      margin <- stats::qnorm(1 - alpha / 2) / sqrt(n - 3L)
      rho_L <- pmin(pmax(-1, pmin(1, tanh(z - margin))), rho_hat)
      rho_U <- pmax(pmax(-1, pmin(1, tanh(z + margin))), rho_hat)
    }
  }

  list(sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U, rho_witness = rho_hat)
}

#' Per-feature Pearson correlation with a shared vector, zero-variance safe
#' @keywords internal
.feature_correlations <- function(X, u) {
  n <- length(u)
  u_c <- u - mean(u)
  sd_u <- sqrt(sum(u_c^2) / (n - 1L))
  X_c <- X - rowMeans(X)
  sd_x <- sqrt(rowSums(X_c^2) / (n - 1L))
  if (!is.finite(sd_u) || sd_u <= 0) {
    return(rep(0, nrow(X)))
  }
  covariance <- rowSums(X_c * rep(u_c, each = nrow(X))) / (n - 1L)
  rho <- covariance / (sd_x * sd_u)
  rho[!is.finite(rho)] <- 0
  pmax(-1, pmin(1, rho))
}
