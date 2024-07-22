run_analysis <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000) {
  # Load necessary libraries
  library(progress)
  library(doSNOW)
  library(foreach)
  library(MCMCpack)
  library(stats)
  library(profvis)  

  # Define global handlers for progress bars
  handlers(global = TRUE)

  # Get the number of columns (N) and rows (D) in the input matrix Y
  N <- ncol(Y)
  D <- nrow(Y)

  # Create a sequence for rho1, rho2, and x based on the given bounds
  rho1 <- seq(-rhobound, rhobound, by = 0.05)
  rho2 <- seq(-rhobound, rhobound, by = 0.05)
  x <- seq(0.05, 0.2, by = 0.01)
  
  # Generate all combinations of rho1, rho2, and x
  pars <- expand.grid(rho1, rho2, x)
  colnames(pars) <- c("rho1", "rho2", "x")
  
  # Print priors
  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Rho bounds:\n")
  print(rhobound)
  cat("Bootstrap sample size (S):\n")
  print(S)
  cat("Parameter grid (rho1, rho2, x):\n")
  print(pars)
  
  # Initialize a results matrix to store the results for each pair
  results <- matrix(list(), D, D)
  
  # Define the objective function used in the optimization
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
    return(sigma)
  }

  # Function to format elapsed time in a user-friendly format
  format_elapsed_time <- function(elapsed_time) {
    total_seconds <- as.numeric(elapsed_time, units = "secs")
    
    seconds <- total_seconds %% 60
    minutes <- (total_seconds %/% 60) %% 60
    hours <- (total_seconds %/% 3600) %% 24
    days <- total_seconds %/% 86400
    
    result <- c()
    
    if (days > 0) {
      result <- c(result, paste(days, "days"))
    }
    if (hours > 0) {
      result <- c(result, paste(hours, "hours"))
    }
    if (minutes > 0) {
      result <- c(result, paste(minutes, "minutes"))
    }
    if (seconds > 0 || length(result) == 0) {
      result <- c(result, paste(round(seconds, 2), "seconds"))
    }
    
    return(paste(result, collapse = ", "))
  }

  # Register the parallel backend
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  
  # Ensure the cluster is stopped after the function exits
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
  
  # Record the start time for profiling
  start_time <- Sys.time()
  
  # Total number of pairs to process
  total_pairs <- D * (D + 1) / 2
  
  # Generate all pairs of indices and add diagonal pairs
  pair_indices <- combn(D, 2, simplify = FALSE)
  pair_indices <- c(pair_indices, lapply(1:D, function(x) c(x, x)))  # Add diagonal pairs

  # Initialize progress bar
  pb <- progress::progress_bar$new(total = total_pairs, format = "  running [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 60)

  # Function to update progress bar
  progress <- function(n) {
    pb$tick()
  } 
  
  # Options for foreach to include progress updates
  opts <- list(progress = progress)

  # Precompute bootstrap samples
  bootstrap_samples <- replicate(S, sample(1:N, replace = TRUE), simplify = FALSE)
  
  # Run the analysis with profiling
  profvis({
    results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .combine = 'rbind', .options.snow = opts) %dopar% {
      d1 <- pair[1]
      d2 <- pair[2]
    
      minsigma_values <- numeric(S)
      maxsigma_values <- numeric(S)
      
      # Use parallel foreach for the inner loop
      results <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack')) %dopar% {
        Yboot <- Y[, bootstrap_samples[[s]]]
    
        # Find the minimum sigma
        res_min <- optim(par = c(0, 0, 0.1), fn = objective_function, Yboot = Yboot, d1 = d1, d2 = d2, alpha = alpha, 
                         method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, 0.2))
    
        # Find the maximum sigma by negating the objective function
        res_max <- optim(par = c(0, 0, 0.1), fn = function(params, Yboot, d1, d2, alpha) {
          -objective_function(params, Yboot, d1, d2, alpha)
        }, Yboot = Yboot, d1 = d1, d2 = d2, alpha = alpha, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, 0.2))
        
        c(minsigma = res_min$value, maxsigma = -res_max$value)
      }
    
      minsigma_values <- results[, "minsigma"]
      maxsigma_values <- results[, "maxsigma"]
    
      # Sort the results
      sortedmin <- sort(minsigma_values)
      sortedmax <- sort(maxsigma_values)
      
      # Compute the minimum, maximum, and confidence intervals
      minsigma <- min(sortedmin)
      maxsigma <- max(sortedmax)
      cilower <- quantile(sortedmin, probs = 0.025)
      ciupper <- quantile(sortedmax, probs = 0.975)
    
      # Compute finite sample covariance
      finitesamplecovariance <- stats::cov(log(Y)[d1,], log(Y)[d2,])
      
      list(d1 = d1, d2 = d2, cilower = cilower, ciupper = ciupper, minsigma = minsigma, maxsigma = maxsigma, finitesamplecovariance = finitesamplecovariance)
    }
  })
  
  # Combine results
  final_results <- do.call(rbind, results_list)
  
  # Create a data frame to store formatted results
  formatted_results <- data.frame(
    comparison = character(),
    cilower = numeric(),
    ciupper = numeric(),
    minsigma = numeric(),
    maxsigma = numeric(),
    finitesamplecovariance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Format results into a data frame
  for (res in final_results) {
    d1 <- res$d1
    d2 <- res$d2
    comparison <- paste(rownames(Y)[d1], ":", rownames(Y)[d2], sep = "")
    cilower <- res$cilower
    ciupper <- res$ciupper
    minsigma <- res$minsigma
    maxsigma <- res$maxsigma
    finitesamplecovariance <- res$finitesamplecovariance
    
    formatted_results <- rbind(formatted_results, data.frame(
      comparison = comparison,
      cilower = cilower,
      ciupper = ciupper,
      minsigma = minsigma,
      maxsigma = maxsigma,
      finitesamplecovariance = finitesamplecovariance,
      stringsAsFactors = FALSE
    ))
  }

  # Remove row names
  rownames(formatted_results) <- NULL
                                         
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(formatted_results)
}
