run_analysis <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000) {
  library(progress)
  library(doSNOW)
  library(foreach)
  library(MCMCpack)
  library(stats)
  
  handlers(global = TRUE)
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
      rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
    }
    rWpara <- log(rWpara)
    
    a <- var(rWpara[d1, ])
    b <- var(rWpara[d2, ])
    c <- cor(rWpara[d1, ], rWpara[d2, ])
    
    sigma <- a * b * c + a * x * rho1 + b * x * rho2 + x^2
    return(min(sigma))
  }

  # Register the parallel backend
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  # Verify cluster registration
  if(!foreach::getDoParRegistered()) {
    stop("Parallel backend is not registered.")
  } else {
    cat("Parallel backend is registered.\n")
  }
  
  # Check the number of workers
  num_workers <- foreach::getDoParWorkers()
  cat("Number of workers: ", num_workers, "\n")
  
  start_time <- Sys.time()
  
  total_pairs <- D * (D + 1) / 2
  
  pair_indices <- combn(D, 2, simplify = FALSE)
  pair_indices <- c(pair_indices, lapply(1:D, function(x) c(x, x)))  # Add diagonal pairs

  pb <- progress::progress_bar$new(total = total_pairs, format = "  running [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 60)

  progress <- function(n) {
    pb$tick()
  } 
  
  opts <- list(progress = progress)
  
  results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
    d1 <- pair[1]
    d2 <- pair[2]
    
    minmaxsigma <- matrix(NA, S, 2)
    for (s in 1:S) {
      Yboot <- Y[, sample(1:N, replace = TRUE)]
      
      res <- optim(par = c(0, 0, 0.1), fn = objective_function, Yboot = Yboot, d1 = d1, d2 = d2, alpha = alpha, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, 0.2))
      
      minmaxsigma[s, ] <- c(min(res$value), max(res$value))
    }
    
    sortedmin <- sort(minmaxsigma[, 1])
    sortedmax <- sort(minmaxsigma[, 2])
    
    cilower <- quantile(sortedmin, probs = 0.025)
    ciupper <- quantile(sortedmax, probs = 0.975)

    finitesamplecovariance <- stats::cov(log(Y)[d1,], log(Y)[d2,])
    
    list(d1 = d1, d2 = d2, cilower = cilower, ciupper = ciupper, finitesamplecovariance = finitesamplecovariance)
  }
  
  formatted_results <- data.frame(
    comparison = character(),
    cilower = numeric(),
    ciupper = numeric(),
    finitesamplecovariance = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (res in results_list) {
    d1 <- res$d1
    d2 <- res$d2
    comparison <- paste(d1, ":", d2, sep = "")
    cilower <- res$cilower
    ciupper <- res$ciupper
    finitesamplecovariance <- res$finitesamplecovariance
    
    formatted_results <- rbind(formatted_results, data.frame(
      comparison = comparison,
      cilower = cilower,
      ciupper = ciupper,
      finitesamplecovariance = finitesamplecovariance,
      stringsAsFactors = FALSE
    ))
  }

  # Remove row names
  rownames(formatted_results) <- NULL
                                         
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste("Total time taken:", elapsed_time))
  
  return(formatted_results)
}
