#' Simulate Data
#'
#' This function simulates data for given dimensions and sample size, with specified sequencing depth.
#'
#' @param D An integer specifying the number of dimensions.
#' @param N An integer specifying the number of samples.
#' @param seq.depth A numeric value specifying the sequencing depth.
#' @return A list containing the simulated data matrix \code{Y}, the log-transformed matrix \code{logW}, and the covariance matrix \code{Sigma}.
#' @examples
#' # Example usage:
#' simulated_data <- simulate_data(3, 100, 1000)
#' @export
simulate_data <- function(D, N, seq.depth) {
  Sigma <- diag(D)
  Sigma[2, 1] <- Sigma[1, 2] <- -0.5
  Sigma[3, 1] <- Sigma[1, 3] <- +0.5

  logW <- rmvnorm(N, rep(0, D), Sigma)
  W <- exp(logW)
  Wpara <- t(miniclo(t(W)))
  Y <- Wpara * seq.depth

  return(list(Y = Y, logW = logW, Sigma = Sigma))
}
