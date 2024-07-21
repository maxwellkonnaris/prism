#' Run Analysis on All Pairwise Taxa
#'
#' @param Y Data matrix
#' @param alpha Alpha parameter
#' @param rhobound Rho bound
#' @param S Number of simulations
#' @return A matrix of lists containing confidence intervals and true values for all pairs of taxa
run_analysis <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000, x=NULL) {
  N <- ncol(Y)
  D <- nrow(Y)
  
  rho1 <- seq(-rhobound, rhobound, by = 0.05)
  rho2 <- seq(-rhobound, rhobound, by = 0.05)
  x <- seq(0.05, 0.2, by = 0.01)
  pars <- expand.grid(rho1, rho2, x)
  colnames(pars) <- c("rho1", "rho2", "x")
  
  results <- matrix(list(), D, D)
  
  objective_function <- function(params, Yboot, d1, d2, alpha) {
    rho1 <- params[1]
    rho2 <- params[2]
    x <- params[3]
    
    rWpara <- matrix(NA, D, N)
    for (n in 1:N) {
      rWpara[, n] <- rdirichlet(1, Yboot[, n] + alpha)
    }
    rWpara <- log(rWpara)
    
    a <- var(rWpara[d1,])
    b <- var(rWpara[d2,])
    c <- cor(rWpara[d1,], rWpara[d2,])
    
    sigma <- a * b * c + a * x * rho1 + b * x * rho2 + x^2
    return(min(sigma))
  }
  
  start_time <- Sys.time()
  
  counter <- 0
  pb <- pbapply::startpb(0, (D-1) * D / 2, style = 3)
  
  for (d1 in 1:(D-1)) {
    for (d2 in (d1+1):D) {
      minmaxsigma <- matrix(NA, S, 2)
      for (s in 1:S) {
        Yboot <- Y[, sample(1:N, replace = TRUE)]
        
        res <- optim(par = c(0, 0, 0.1), fn = objective_function, Yboot = Yboot, d1 = d1, d2 = d2, alpha = alpha, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, 0.2))
        
        minmaxsigma[s, ] <- c(res$value, res$value)
      }
      
      sortedmin <- sort(minmaxsigma[, 1])
      sortedmax <- sort(minmaxsigma[, 2])
      
      cilower <- quantile(sortedmin, probs = 0.025)
      ciupper <- quantile(sortedmax, probs = 0.975)
      
      results[[d1, d2]] <- list(c(cilower, ciupper), minmaxsigma)
      
      counter <- counter + 1
      pbapply::setpb(pb, counter)
    }
  }
  
  pbapply::closepb(pb)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste("Total time taken:", elapsed_time))
  
  return(results)
}
