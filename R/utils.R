`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

log_props_from_counts <- function(counts, pseudocount = 0.5) {
  Y <- as.matrix(counts)
  if (any(Y < 0)) stop("counts must be nonnegative.")
  Yp <- Y + as.numeric(pseudocount)
  colsum <- colSums(Yp)
  if (any(colsum <= 0)) stop("Each sample must have positive total after pseudocount.")
  P <- sweep(Yp, 2, colsum, "/")
  log(P)
}

rel_cov_from_log_props <- function(log_props, ddof = 1) {
  X <- as.matrix(log_props)
  # X is D x N (variables x samples); cov expects N x D
  stats::cov(t(X))
}

.estimate_scale_log_bounds <- function(log_proportions, scale_log, ci_level = 0) {
  X <- as.matrix(log_proportions)
  u <- as.numeric(scale_log)
  n <- length(u)
  if (ncol(X) != n || n < 2L || any(!is.finite(X)) || any(!is.finite(u))) {
    stop("log_proportions and scale_log must be finite and sample-aligned.", call. = FALSE)
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level < 0 || ci_level >= 1) {
    stop("ci_level must be a single finite number in [0, 1).", call. = FALSE)
  }

  sigma_hat <- stats::sd(u)
  if (!is.finite(sigma_hat) || sigma_hat < 0) sigma_hat <- 0
  u_centered <- u - mean(u)
  X_centered <- X - rowMeans(X)
  sd_u <- sqrt(sum(u_centered^2) / (n - 1L))
  sd_x <- sqrt(rowSums(X_centered^2) / (n - 1L))
  if (!is.finite(sd_u) || sd_u <= 0 || n < 3L) {
    rho_hat <- rep(0, nrow(X))
  } else {
    covariance <- rowSums(X_centered * rep(u_centered, each = nrow(X))) / (n - 1L)
    rho_hat <- covariance / (sd_x * sd_u)
    rho_hat[!is.finite(rho_hat)] <- 0
    rho_hat <- pmax(-1, pmin(1, rho_hat))
  }

  if (ci_level == 0) {
    return(list(
      sigma_L = sigma_hat,
      sigma_U = sigma_hat,
      rho_L = rho_hat,
      rho_U = rho_hat,
      rho_witness = rho_hat
    ))
  }

  alpha <- 1 - ci_level
  df <- n - 1L
  sigma_L <- sqrt(df * sigma_hat^2 / stats::qchisq(1 - alpha / 2, df = df))
  sigma_U <- sqrt(df * sigma_hat^2 / stats::qchisq(alpha / 2, df = df))
  if (n < 4L) {
    rho_L <- rho_U <- rho_hat
  } else {
    epsilon <- sqrt(.Machine$double.eps)
    z <- atanh(pmax(-1 + epsilon, pmin(1 - epsilon, rho_hat)))
    margin <- stats::qnorm(1 - alpha / 2) / sqrt(n - 3L)
    rho_L <- pmax(-1, pmin(1, tanh(z - margin)))
    rho_U <- pmax(-1, pmin(1, tanh(z + margin)))
    rho_L <- pmin(rho_L, rho_hat)
    rho_U <- pmax(rho_U, rho_hat)
  }
  list(
    sigma_L = sigma_L,
    sigma_U = sigma_U,
    rho_L = rho_L,
    rho_U = rho_U,
    rho_witness = rho_hat
  )
}
