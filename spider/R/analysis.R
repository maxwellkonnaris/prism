#' Run Analysis on All Pairwise Taxa
#'
#' @param Y Data matrix
#' @param alpha Alpha parameter
#' @param rhobound Rho bound
#' @param S Number of simulations
#' @param num_cores Number of cores for parallel processing
#' @return A matrix of lists containing confidence intervals and true values for all pairs of taxa
#' @import pbapply
#' @import parallel
#' @import MCMCpack
run_analysis <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000, num_cores = parallel::detectCores() - 1) {
  N <- ncol(Y)
  D <- nrow(Y)
  
  rho1 <- seq(-rhobound, rhobound, by = 0.05)
  rho2 <- seq(-rhobound, rhobound, by = 0.05)
  x <- seq(0.05, 0.2, by = 0.01)
  pars <- base::expand.grid(rho1, rho2, x)
  colnames(pars) <- c("rho1", "rho2", "x")
  
  results <- base::matrix(list(), D, D)
  
  objective_function <- function(params, Yboot, d1, d2, alpha) {
    rho1 <- params[1]
    rho2 <- params[2]
    x <- params[3]
    
    rWpara <- base::matrix(NA, D, N)
    for (n in 1:N) {
      rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
    }
    rWpara <- base::log(rWpara)
    
    a <- stats::var(rWpara[d1, ])
    b <- stats::var(rWpara[d2, ])
    c <- stats::cor(rWpara[d1, ], rWpara[d2, ])
    
    sigma <- a * b * c + a * x * rho1 + b * x * rho2 + x^2
    return(base::min(sigma))
  }
  
  # Register the parallel backend
  cl <- parallel::makeCluster(num_cores)
  parallel::clusterExport(cl, c("Y", "alpha", "rhobound", "S", "objective_function", "D", "N"))
  parallel::clusterEvalQ(cl, library(MCMCpack))
  
  on.exit({
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  start_time <- base::Sys.time()
  
  total_pairs <- D * (D + 1) / 2
  pair_indices <- utils::combn(D, 2, simplify = FALSE)
  pair_indices <- base::c(pair_indices, lapply(1:D, function(x) c(x, x)))  # Add diagonal pairs
  
  results_list <- pbapply::pblapply(pair_indices, function(pair) {
    d1 <- pair[1]
    d2 <- pair[2]
    
    minmaxsigma <- base::matrix(NA, S, 2)
    for (s in 1:S) {
      Yboot <- Y[, base::sample(1:N, replace = TRUE)]
      
      res <- stats::optim(par = c(0, 0, 0.1), fn = objective_function, Yboot = Yboot, d1 = d1, d2 = d2, alpha = alpha, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, 0.2))
      
      minmaxsigma[s, ] <- c(base::min(res$value), base::max(res$value))
    }
    
    sortedmin <- base::sort(minmaxsigma[, 1])
    sortedmax <- base::sort(minmaxsigma[, 2])
    
    cilower <- stats::quantile(sortedmin, probs = 0.025)
    ciupper <- stats::quantile(sortedmax, probs = 0.975)
    
    finitesamplecovariance <- stats::cov(Y[d1,], Y[d2,])
    
    list(cilower = cilower, ciupper = ciupper, finitesamplecovariance = finitesamplecovariance)
  }, cl = cl)
  
  for (i in 1:base::length(pair_indices)) {
    pair <- pair_indices[[i]]
    d1 <- pair[1]
    d2 <- pair[2]
    results[[d1, d2]] <- results_list[[i]]
  }
  
  end_time <- base::Sys.time()
  elapsed_time <- end_time - start_time
  base::print(base::paste("Total time taken:", elapsed_time))
  
  return(results)
}
