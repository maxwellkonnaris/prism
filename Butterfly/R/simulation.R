#' Simulate Data
#'
#' @param D Number of dimensions
#' @param N Number of samples
#' @param seq.depth Sequencing depth
#' @return Simulated data matrix
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

