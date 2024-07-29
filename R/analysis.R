#' Run Analysis on All Pairwise Taxa
#'
#' This function runs a bootstrapped analysis on the input data matrix \code{Y}.
#'
#' @param Y A matrix of data with observations in columns and variables in rows.
#' @param alpha A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros.
#' @param rhobound A numeric value specifying the bound for the \code{rho1} and \code{rho2} parameters. Defaults to 0.8.
#' @param S An integer specifying the number of bootstrap samples. Increase to reduce Monte Carlo error. Defaults to 1000.
#' @param variance A boolean specifying whether to include diagonals of variance-covariance matrices (the variance). Defaults to FALSE.
#' @param upperx A numeric value specifying the upper bound of the variance of the scale. Defaults to 1.0.
#' @return A list of dataframes containing the results of the analysis including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances. \code{final_results} contains the data intended for forest_plot() and \code{all_inner_results} contains the data intended for sigmaplot().
#' @import progress
#' @import progressr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MCMCpack
#' @import stats
#' @import utils
#' @examples
#' # Example usage (You could also use the simulation function provided to generate sample data):
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10)
#' results <- estimate_covariance(Y)
#' @export

estimate_covariance <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 1.0, S = 1000, variance=FALSE, upperx = 2.0) {

  # Record the start time for profiling
  start_time <- Sys.time()

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

  # Define global handlers for progress bars
  handlers(global = TRUE)

  # Check if Y is a matrix and has appropriate dimensions
  if (!is.matrix(Y) || nrow(Y) < 2 || ncol(Y) < 2) {
    stop("Y must be a matrix with at least 2 rows and 2 columns.")
  }

  # Get the number of columns (N) and rows (D) in the input matrix Y
  N <- ncol(Y)
  D <- nrow(Y)

  # Create a sequence for rho1, rho2, and x based on the given bounds
  rho1 <- seq(-rhobound, rhobound, by = 0.05)
  rho2 <- seq(-rhobound, rhobound, by = 0.05)
  x <- seq(0.05, upperx, by = 0.025)
  
  # Generate all combinations of rho1, rho2, and x
  pars <- expand.grid(rho1, rho2, x)
  colnames(pars) <- c("rho1", "rho2", "x")
  
  # Print priors
  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Rho bounds:\n")
  print(paste0(-rhobound,":",rhobound))
  # Print the head of the data frame
  cat("Head of the parameter grid:\n")
  print(head(pars))
  # Print the tail of the data frame
  cat("Tail of the parameter grid:\n")
  print(tail(pars))
  cat("Dimensions of supplied matrix:\n")
  print(dim(Y))
  cat("Bootstrap sample size (S):\n")
  print(S)
  
  # Initialize a results matrix to store the results for each pair
  results <- matrix(list(), D, D)
  
  # Define the objective function used in the optimization
  objective_function <- function(params, a, b, c) {
    rho1 <- params[1]
    rho2 <- params[2]
    x <- params[3]
    
    sigma <- a * b * c + a * x * rho1 + b * x * rho2 + x^2
    return(sigma)
  }

  # Register the parallel backend
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  
  # Ensure the cluster is stopped after the function exits
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Verify cluster registration
  if (!foreach::getDoParRegistered()) {
    stop("Parallel backend is not registered.")
  } else {
    # Check the number of workers
    cat("Number of workers/cpus: ", foreach::getDoParWorkers(), "\n")
  }
  
  # Initialize progress bars
  total_pairs <- D * (D + 1) / 2
  pb_precomp <- progress::progress_bar$new(total = S, format = "  Precomputing bootstrap samples [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  pb <- progress::progress_bar$new(total = total_pairs, format = " Generating Sigmas [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  
  # Function to update progress bar
  progress_precomp <- function(n) {
    pb_precomp$tick()
  }
  
  progress <- function(n) {
    pb$tick()
  }
  
  # Options for foreach to include progress updates
  opts_precomp <- list(progress = progress_precomp)
  opts <- list(progress = progress)

  # Parallelize bootstrap precomputation
  bootstrap_samples <- foreach(s = 1:S, .combine = 'c', .options.snow = opts_precomp) %dopar% {
    sample(1:N, replace = TRUE)
  }
  
  # Reshape bootstrap_samples into a list of vectors
  bootstrap_samples <- split(bootstrap_samples, rep(1:S, each = N))
  
  # Generate all pairs of indices and add diagonal pairs
  pair_indices <- combn(D, 2, simplify = FALSE)
  if (variance) {
    pair_indices <- c(pair_indices, lapply(1:D, function(x) c(x, x)))  # Add diagonal pairs
  }

  # Check if pair_indices is populated correctly
  if (length(pair_indices) == 0) {
    stop("Error: pair_indices is not populated correctly. Aborting analysis.")
  }
                                          
  # Run the analysis
  cat("Running sigma estimation")
  
                                           
  results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
    tryCatch({
      d1 <- pair[1]
      d2 <- pair[2]
      comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
    
      minsigma_values <- numeric(S)
      maxsigma_values <- numeric(S)
      a_values <- numeric(S)
      b_values <- numeric(S)
      c_values <- numeric(S)
      
      # Use parallel foreach for the inner loop
      results_inner <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
          Yboot <- Y[, bootstrap_samples[[s]]]
          
          rWpara <- matrix(NA, D, N)
          for (n in 1:N) {
            rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
          }
          rWpara <- log(rWpara)
    
          a <- var(rWpara[d1, ])
          b <- var(rWpara[d2, ])
          c <- cor(rWpara[d1, ], rWpara[d2, ])
          
          # Find the minimum sigma
          res_min <- optim(par = c(0, 0, 0.5), fn = objective_function(params, a, b, c), a = a, b = b, c = c,
                           method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperx))
          
          # Find the maximum sigma by negating the objective function
          res_max <- optim(par = c(0, 0, 0.5), fn = -objective_function(params, a, b, c), a = a, b = b, c = c,
          		   method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperx))
          
          data.frame(
            d1 = d1,
            d2 = d2,
            s = s,
            minsigma = res_min$value,
            maxsigma = -res_max$value,
            min_rho1 = res_min$par[1],
            min_rho2 = res_min$par[2],
            min_x = res_min$par[3],
            max_rho1 = res_max$par[1],
            max_rho2 = res_max$par[2],
            max_x = res_max$par[3],
            a = a,
            b = b,
            c = c,
          )
      }
  
      # Gather the min and max optimized sigmas
      minsigma_values <- results_inner$minsigma
      maxsigma_values <- results_inner$maxsigma
      
      # Sort the results
      sortedmin <- sort(minsigma_values)
      sortedmax <- sort(maxsigma_values)
      
      # Compute the minimum, maximum, and confidence intervals
      minsigma <- min(sortedmin, na.rm = TRUE)
      maxsigma <- max(sortedmax, na.rm = TRUE)
      cilower <- quantile(sortedmin, probs = 0.025, na.rm = TRUE)
      ciupper <- quantile(sortedmax, probs = 0.975, na.rm = TRUE)
      
      # Obtain parameters for the minimum and maximum sigma values
      min_index <- which.min(minsigma_values)
      max_index <- which.max(maxsigma_values)
      
      min_rho1 <- results_inner$min_rho1[min_index]
      min_rho2 <- results_inner$min_rho2[min_index]
      min_x <- results_inner$min_x[min_index]
      
      max_rho1 <- results_inner$max_rho1[max_index]
      max_rho2 <- results_inner$max_rho2[max_index]
      max_x <- results_inner$max_x[max_index]
      
      a <- results_inner$a[min_index]
      b <- results_inner$b[min_index]
      c <- results_inner$c[min_index]
  
      list(
        resultsinner = results_inner,
        results = data.frame(
          comparison = comparison,
          taxa1 = rownames(Y)[d1],
          taxa2 = rownames(Y)[d2],
          ninetyfive_ci_lower = cilower,
          ninetyfive_ci_upper = ciupper,
          minsigma_absolute_minimum_covariance = minsigma,
          maxsigma_absolute_maximum_covariance = maxsigma,
          minsigma_correlation_taxa1_scale = min_rho1,
          minsigma_correlation_taxa2_scale = min_rho2,
          minsigma_scale_variance = min_x,
          maxsigma_correlation_taxa1_scale = max_rho1,
          maxsigma_correlation_taxa2_scale = max_rho2,
          maxsigma_scale_variance = max_x,
          relative_standard_dev_taxa1 = a,
          relative_standard_dev_taxa2 = b,
          relative_covariance = c,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
        message("Error in outer loop: ", e$message)
        return(NULL)
    })
  }
  
  # Combine the results into a data frame, transpose it, remove row names
  final_results <- do.call(rbind, lapply(results_list, function(x) x$results))
  final_results <- as.data.frame(final_results, stringsAsFactors = FALSE)
  rownames(final_results) <- NULL
  
  # Calculate p-values and adjust for multiple hypothesis testing
  final_results <- calculate_pval(final_results)
                                         
  # Combine all inner loop results
  all_inner_results <- do.call(rbind, lapply(results_list, function(x) x$resultsinner))
  all_inner_results <- as.data.frame(all_inner_results)
  all_inner_results$comparison <- paste(rownames(Y)[all_inner_results$d1],rownames(Y)[all_inner_results$d2],sep=":")
  
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(list(final_results = final_results, all_inner_results = all_inner_results))
}
