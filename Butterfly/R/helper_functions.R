#' Multivariate Normal Random Variable Generator
#'
#' @param n Number of samples
#' @param mu Mean vector
#' @param Sigma Covariance matrix
#' @return A matrix of random samples
rmvnorm <- function(n, mu, Sigma) {
  p <- length(mu)
  r <- matrix(rnorm(n * p), p, n)
  L <- t(chol(Sigma))
  r <- L %*% r
  sweep(r, 1, mu, FUN = `+`)
}

#' Check if value is within range
#'
#' @param x Value to check
#' @param l Lower bound
#' @param u Upper bound
#' @return Boolean indicating if value is within range
within <- function(x, l, u) {
  x >= l && x <= u
}

