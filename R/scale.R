#' Fixed scale-SD and correlation bounds
#'
#' Specifies `sigma`/`rho` bounds directly (no per-sample data), for
#' `scale = prism_scale_bounds(...)` in [prism()]. Correlation bounds
#' require scale-SD bounds; both members of a pair must be supplied
#' together.
#'
#' @param sigma_L,sigma_U Optional scale-SD bounds, `0 <= sigma_L <= sigma_U`.
#' @param rho_L,rho_U Optional scalar or feature-level correlation bounds
#'   in `[-1, 1]`, `rho_L <= rho_U`.
#' @return A `prism_scale_bounds` object.
#' @export
prism_scale_bounds <- function(sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL) {
  has_sigma <- !is.null(sigma_L) || !is.null(sigma_U)
  has_rho <- !is.null(rho_L) || !is.null(rho_U)
  if (has_rho && !has_sigma) {
    stop("Correlation bounds require scale-SD bounds (sigma_L, sigma_U).", call. = FALSE)
  }
  if (has_sigma) {
    if (is.null(sigma_L) || is.null(sigma_U)) {
      stop("sigma_L and sigma_U must be supplied together.", call. = FALSE)
    }
    sigma_L <- .assert_scalar_finite(sigma_L, "sigma_L", lower = 0)
    sigma_U <- .assert_scalar_finite(sigma_U, "sigma_U", lower = 0)
    if (sigma_L > sigma_U) stop("sigma_L must be <= sigma_U.", call. = FALSE)
  }
  if (has_rho) {
    if (is.null(rho_L) || is.null(rho_U)) {
      stop("rho_L and rho_U must be supplied together.", call. = FALSE)
    }
    ok <- is.numeric(rho_L) && is.numeric(rho_U) && all(is.finite(rho_L)) && all(is.finite(rho_U)) &&
      all(rho_L >= -1) && all(rho_L <= 1) && all(rho_U >= -1) && all(rho_U <= 1) && all(rho_L <= rho_U)
    if (!ok) stop("rho_L/rho_U must be finite, within [-1, 1], and rho_L <= rho_U.", call. = FALSE)
  }
  structure(
    list(sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U),
    class = "prism_scale_bounds"
  )
}

#' Configure sample-level log-scale measurements
#'
#' PRISM requires sample-level scale measurements to be supplied on the log
#' scale. Bare numeric vectors and matrices passed to [prism()] are interpreted
#' as already log transformed. This constructor is only needed to change the
#' within-draw estimation settings.
#'
#' @param values Finite numeric vector (length N) or matrix (N rows, one
#'   column per technical replicate), already on the log scale.
#' @param ci_level Interval coverage in `[0, 1)` for the scale-SD/
#'   correlation confidence interval estimated within each draw. The default
#'   `0` uses point estimates. Positive values are experimental and require
#'   `experimental_ci = TRUE`.
#' @param estimate_rho Estimate taxon-scale correlation. If `FALSE`, only
#'   scale-SD is estimated.
#' @param lower_zero Force the lower scale-SD endpoint to `0` instead of
#'   the chi-square lower quantile (see [prism()]).
#' @param experimental_ci Explicitly enable analytic chi-square/Fisher-z
#'   intervals inside each bootstrap draw. This sensitivity-analysis mode is
#'   not the default bootstrap procedure.
#' @return A `prism_scale_log` object.
#' @export
prism_scale_log <- function(values, ci_level = 0, estimate_rho = TRUE,
                            lower_zero = TRUE, experimental_ci = FALSE) {
  .new_prism_scale_data(
    values, ci_level, estimate_rho, lower_zero, experimental_ci
  )
}

.new_prism_scale_data <- function(values, ci_level, estimate_rho, lower_zero,
                                  experimental_ci) {
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("values must be a finite numeric vector or matrix.", call. = FALSE)
  }
  ci_level <- .assert_scalar_finite(ci_level, "ci_level", lower = 0, upper = 1, upper_inclusive = FALSE)
  if (!is.logical(estimate_rho) || length(estimate_rho) != 1L || is.na(estimate_rho)) {
    stop("estimate_rho must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(lower_zero) || length(lower_zero) != 1L || is.na(lower_zero)) {
    stop("lower_zero must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(experimental_ci) || length(experimental_ci) != 1L || is.na(experimental_ci)) {
    stop("experimental_ci must be TRUE or FALSE.", call. = FALSE)
  }
  if (ci_level > 0 && !experimental_ci) {
    stop(
      "ci_level > 0 enables experimental within-draw analytic intervals; ",
      "set experimental_ci = TRUE to use this sensitivity analysis.",
      call. = FALSE
    )
  }
  structure(
    list(
      values = values, ci_level = ci_level, estimate_rho = estimate_rho,
      lower_zero = lower_zero, experimental_ci = experimental_ci
    ),
    class = "prism_scale_log"
  )
}

#' Validate sample-level log-scale input
#'
#' `scale` is either one log-scale measurement per sample (a length-N vector) or
#' a matrix of replicate measurements per sample (N rows, one column per
#' replicate). Both are normalized to an N x R finite matrix so the rest
#' of the package only has to handle one shape.
#'
#' @param scale Finite numeric vector (length N) or matrix (N rows), already
#'   on the log scale.
#' @param N Expected number of samples.
#' @return A finite N x R matrix.
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
  if (!is.numeric(scale) || any(!is.finite(scale))) {
    stop("scale must be finite.", call. = FALSE)
  }
  scale
}

#' Bootstrap the subject-level mean log scale
#'
#' For samples with one measurement this returns that value on the log scale.
#' For technical replicates, it resamples all replicates with replacement
#' within each selected subject and averages them on the log scale. This
#' propagates uncertainty in the subject-level mean without discarding all but
#' one replicate.
#'
#' @param scale_mat N x R matrix from `.validate_scale_input()`.
#' @param sample_index Integer sample indices for this draw (length N,
#'   with replacement under the bootstrap).
#' @return A numeric vector of log-scale values, length `length(sample_index)`.
#' @keywords internal
.select_scale_log <- function(scale_mat, sample_index) {
  R <- ncol(scale_mat)
  values <- scale_mat[sample_index, , drop = FALSE]
  if (R == 1L) return(as.numeric(values[, 1L]))
  vapply(
    seq_len(nrow(values)),
    function(i) mean(values[i, sample.int(R, R, replace = TRUE)]),
    numeric(1)
  )
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
#' If the scale measurements have zero or numerically zero variance, the
#' function returns `sigma_L = sigma_U = 0` and leaves all rho fields `NULL`:
#' correlation with a constant variable is undefined but irrelevant when
#' multiplied by zero scale SD.
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
#'   The rho fields are `NULL` when `estimate_rho = FALSE` or scale variance
#'   is zero.
#' @keywords internal
.estimate_scale_log_bounds <- function(log_proportions, scale_log, ci_level = 0,
                                        estimate_rho = TRUE, lower_zero = TRUE) {
  X <- as.matrix(log_proportions)
  u <- as.numeric(scale_log)
  n <- length(u)
  if (ncol(X) != n || n < 2L || any(!is.finite(X)) || any(!is.finite(u))) {
    stop("log_proportions and scale_log must be finite and sample-aligned.", call. = FALSE)
  }
  ci_level <- .assert_scalar_finite(ci_level, "ci_level", lower = 0, upper = 1, upper_inclusive = FALSE)

  sigma_hat <- stats::sd(u)
  if (.numerically_zero_variance(u, sigma_hat^2)) {
    return(list(
      sigma_L = 0, sigma_U = 0,
      rho_L = NULL, rho_U = NULL, rho_witness = NULL
    ))
  }
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

#' Per-feature Pearson correlation with a shared vector
#'
#' A correlation is undefined when either variable has zero variance. This
#' function therefore errors instead of substituting an arbitrary correlation.
#' @keywords internal
.feature_correlations <- function(X, u) {
  n <- length(u)
  u_c <- u - mean(u)
  sd_u <- sqrt(sum(u_c^2) / (n - 1L))
  X_c <- X - rowMeans(X)
  sd_x <- sqrt(rowSums(X_c^2) / (n - 1L))
  if (.numerically_zero_variance(u, sd_u^2)) {
    stop(
      "The log-scale measurements have zero or numerically zero marginal variance; ",
      "correlations are undefined.",
      call. = FALSE
    )
  }
  bad_x <- vapply(
    seq_len(nrow(X)),
    function(i) .numerically_zero_variance(X[i, ], sd_x[i]^2),
    logical(1)
  )
  if (any(bad_x)) {
    stop(
      "Log-composition feature(s) ", paste(which(bad_x), collapse = ", "),
      " have zero or numerically zero marginal variance; correlations are undefined.",
      call. = FALSE
    )
  }
  covariance <- rowSums(X_c * rep(u_c, each = nrow(X))) / (n - 1L)
  rho <- covariance / (sd_x * sd_u)
  if (any(!is.finite(rho))) {
    stop("Feature correlations could not be estimated as finite values.", call. = FALSE)
  }
  pmax(-1, pmin(1, rho))
}
