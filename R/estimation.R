#' Run Analysis on All Pairwise Taxa
#'
#' This function runs a bootstrapped analysis on the input data matrix \code{Y}.
#'
#' @param Y A matrix of data with observations in columns and variables in rows.
#' @param alpha A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros.
#' @param rhobound A numeric value specifying the bound for the \code{rho1} and \code{rho2} parameters. Defaults to 0.9.
#' @param S An integer specifying the number of bootstrap samples. Increase to reduce Monte Carlo error. Defaults to 1000.
#' @param upperscalevariance A numeric value specifying the upper bound of the variance of the scale. Defaults to 1.0.
#' @return A list of dataframes containing the results of the analysis including estimated 95% confidence intervals, minimum and maximum values for estimated absolute covariance, and parameters. \code{final_results} contains the data intended for forest_plot() and \code{all_inner_results} contains the data intended for sigmaplot().
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

estimate_covariance <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.9, S = 1000, upperscalevariance = 1.0) {

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
  taxa1scalecorrelation <- seq(-rhobound, rhobound, by = 0.05)
  taxa2scalecorrelation <- seq(-rhobound, rhobound, by = 0.05)
  scalevariance <- seq(0.05, upperscalevariance, by = 0.025)
  
  # Generate all combinations of rho1, rho2, and x
  pars <- expand.grid(taxa1scalecorrelation, taxa2scalecorrelation, scalevariance)
  colnames(pars) <- c("Taxa1-Scale Correlation", "Taxa2-Scale Correlation", "Scale Variance")
  
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
  objective_function <- function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {
    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalevariance <- params[3]
    
    sigma <- taxa1relativesd * taxa2relativesd * relativecorrelation + taxa1relativesd * scalevariance * taxa1scalecorrelation + taxa2relativesd * scalevariance * taxa2scalecorrelation + scalevariance^2
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
  total_pairs <- D * (D - 1) / 2
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

  # Check if pair_indices is populated correctly
  if (length(pair_indices) == 0) {
    stop("Error: pair_indices is not populated correctly. Aborting analysis.")
  }
                                          
  # Run the analysis
  cat("Running sigma estimation")
                                           
  results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
      d1 <- pair[1]
      d2 <- pair[2]
      comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
    
      minsigma_values <- numeric(S)
      maxsigma_values <- numeric(S)
      
      # Use parallel foreach for the inner loop
      results_inner <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
          Yboot <- Y[, bootstrap_samples[[s]]]
          
          rWpara <- matrix(NA, D, N)
          for (n in 1:N) {
            rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
          }
          rWpara <- log(rWpara)
    
          taxa1relativesd <- var(rWpara[d1, ])
          taxa2relativesd <- var(rWpara[d2, ])
          relativecorrelation <- cor(rWpara[d1, ], rWpara[d2, ])
          relativecovariance <- cov(rWpara[d1, ], rWpara[d2, ])
          
          # Find the minimum sigma
          res_min <- optim(par = c(0, 0, 0.1), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = objective_function, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperscalevariance))
          
          # Find the maximum sigma by negating the objective function
          res_max <- optim(par = c(0, 0, 0.1), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {-objective_function(params, taxa1relativesd, taxa2relativesd, relativecorrelation)}, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperscalevariance))
          
          data.frame(
            d1 = d1,
            d2 = d2,
            s = s,
            minsigma_absolute_minimum_covariance = res_min$value,
            maxsigma_absolute_maximum_covariance = -res_max$value,
            minsigma_correlation_relativetaxa1_scale = res_min$par[1],
            minsigma_correlation_relativetaxa2_scale = res_min$par[2],
            minsigma_scale_variance = res_min$par[3],
            maxsigma_correlation_relativetaxa1_scale = res_max$par[1],
            maxsigma_correlation_relativetaxa2_scale = res_max$par[2],
            maxsigma_scale_variance = res_max$par[3],
            taxa1relativesd = taxa1relativesd,
            taxa2relativesd = taxa2relativesd,
            relativecorrelation = relativecorrelation,
            relativecovariance = relativecovariance
          )
      }
      
      # Gather the min and max optimized sigmas
      minsigma_values <- results_inner$minsigma_absolute_minimum_covariance
      maxsigma_values <- results_inner$maxsigma_absolute_maximum_covariance
      
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
      
      min_rho1 <- results_inner$minsigma_correlation_relativetaxa1_scale[min_index]
      min_rho2 <- results_inner$minsigma_correlation_relativetaxa2_scale[min_index]
      min_x <- results_inner$minsigma_scale_variance[min_index]
      
      max_rho1 <- results_inner$maxsigma_correlation_relativetaxa1_scale[max_index]
      max_rho2 <- results_inner$maxsigma_correlation_relativetaxa2_scale[max_index]
      max_x <- results_inner$maxsigma_scale_variance[max_index]
      
      taxa1relativesd <- results_inner$taxa1relativesd[min_index]
      taxa2relativesd <- results_inner$taxa2relativesd[min_index]
      relativecorrelation <- results_inner$relativecorrelation[min_index]
      relativecovariance <- results_inner$relativecovariance[min_index]
  
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
          minsigma_correlation_relativetaxa1_scale = min_rho1,
          minsigma_correlation_relativetaxa2_scale = min_rho2,
          minsigma_scale_variance = min_x,
          maxsigma_correlation_relativetaxa1_scale = max_rho1,
          maxsigma_correlation_relativetaxa2_scale = max_rho2,
          maxsigma_scale_variance = max_x,
          relative_standard_dev_taxa1 = taxa1relativesd,
          relative_standard_dev_taxa2 = taxa2relativesd,
          relative_correlation = relativecorrelation,
          relative_covariance = relativecovariance,
          stringsAsFactors = FALSE
        )
      )
  }
  
  # Combine the results into a data frame, transpose it, remove row names
  final_results <- do.call(rbind, lapply(results_list, function(x) x$results))
  final_results <- as.data.frame(final_results, stringsAsFactors = FALSE)
  rownames(final_results) <- NULL

  # Calculate pvals
  final_results <- calculate_pval(final_results)
                                         
  # Combine all inner loop results
  all_inner_results <- do.call(rbind, lapply(results_list, function(x) x$resultsinner))
  all_inner_results <- as.data.frame(all_inner_results)
  all_inner_results$comparison <- paste(rownames(Y)[all_inner_results$d1],rownames(Y)[all_inner_results$d2],sep=":")
    # Ensure N and D are exported to the parallel environment
  clusterExport(cl, c("N", "D"))
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(list(final_results = final_results, all_inner_results = all_inner_results))
}

#' Run Analysis on All Pairwise Taxa with Bootstrap Convergence Diagnostics
#'
#' This function runs a bootstrapped analysis on the input data matrix \code{Y}, estimating covariance for all pairwise comparisons of taxa. It includes convergence diagnostics by incrementally increasing the number of bootstrap samples.
#'
#' @param Y A matrix of data with observations in columns and variables in rows.
#' @param alpha A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros.
#' @param rhobound A numeric value specifying the bound for the \code{rho1} and \code{rho2} parameters. Defaults to 0.9.
#' @param S An integer specifying the initial number of bootstrap samples. Defaults to 1000.
#' @param upperscalevariance A numeric value specifying the upper bound of the variance of the scale. Defaults to 1.0.
#' @return A list containing:
#' \item{convergence_results}{A list of results for each incrementally increased number of bootstrap samples (100, 500, 1000, 2000, 5000, 10000), including estimated 95% confidence intervals, minimum and maximum values for estimated absolute covariance, and parameters.}
#' \item{combined_results}{A dataframe combining the results from all the incremental bootstrap samples, used for plotting convergence diagnostics.}
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks the input matrix \code{Y} for appropriate dimensions.
#'   \item Defines the sequence for \code{rho1}, \code{rho2}, and scale variance based on the specified bounds.
#'   \item Registers a parallel backend to utilize multiple cores for computation.
#'   \item Runs bootstrap analysis for the specified number of samples (\code{S}), generating bootstrap samples and calculating pairwise covariances.
#'   \item Repeats the bootstrap analysis for incremental sample sizes (100, 500, 1000, 2000, 5000, 10000) to perform convergence diagnostics.
#'   \item Combines the results from all increments into a single dataframe for plotting.
#'   \item Plots the convergence diagnostics to visualize the stability of the estimates as the number of bootstrap samples increases.
#' }
#' @import progress
#' @import progressr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MCMCpack
#' @import stats
#' @import utils
#' @import ggplot2
#' @examples
#' # Example usage (You could also use the simulation function provided to generate sample data):
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10)
#' results <- estimate_covariance_convergence(Y)
#' @export
estimate_covariance_convergence <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.9, S = 1000, upperscalevariance = 1.0) {
  
  # Record the start time for profiling
  start_time <- Sys.time()

  # Record the times for the increments of bootstrap sample sizes
  bootstrap_times <- data.frame(S = numeric(), Time = numeric())

  # Function to format elapsed time in a user-friendly format
  format_elapsed_time <- function(elapsed_time) {
    total_seconds <- as.numeric(elapsed_time, units = "secs")
    seconds <- total_seconds %% 60
    minutes <- (total_seconds %/% 60) %% 60
    hours <- (total_seconds %/% 3600) %% 24
    days <- total_seconds %/% 86400
    result <- c()
    if (days > 0) result <- c(result, paste(days, "days"))
    if (hours > 0) result <- c(result, paste(hours, "hours"))
    if (minutes > 0) result <- c(result, paste(minutes, "minutes"))
    if (seconds > 0 || length(result) == 0) result <- c(result, paste(round(seconds, 2), "seconds"))
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
  taxa1scalecorrelation <- seq(-rhobound, rhobound, by = 0.05)
  taxa2scalecorrelation <- seq(-rhobound, rhobound, by = 0.05)
  scalevariance <- seq(0.05, upperscalevariance, by = 0.025)
  
  # Generate all combinations of rho1, rho2, and x
  pars <- expand.grid(taxa1scalecorrelation, taxa2scalecorrelation, scalevariance)
  colnames(pars) <- c("Taxa1-Scale Correlation", "Taxa2-Scale Correlation", "Scale Variance")
  
  # Print priors
  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Rho bounds:\n")
  print(paste0(-rhobound, ":", rhobound))
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
  objective_function <- function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {
    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalevariance <- params[3]
    sigma <- taxa1relativesd * taxa2relativesd * relativecorrelation + taxa1relativesd * scalevariance * taxa1scalecorrelation + taxa2relativesd * scalevariance * taxa2scalecorrelation + scalevariance^2
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
  total_pairs <- D * (D - 1) / 2
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

  # Function to run bootstrap analysis for a given number of samples
  run_bootstrap_analysis <- function(S, Y, N, D, alpha, rhobound, upperscalevariance) {

    # Run the analysis
    cat("Running sigma estimation for convergence bootstrap diagnostics")
    cat("Number of columns (N):", N, "\n")
    cat("Number of rows (D):", D, "\n")
    
    # Start timing
    bootstrap_start_time <- Sys.time()

    # Ensure N and D are exported to the parallel environment
    clusterExport(cl, c("N", "D"))
    
    # Parallelize bootstrap precomputation
    bootstrap_samples <- foreach(s = 1:S, .combine = 'c', .options.snow = opts_precomp) %dopar% {
      sample(1:N, replace = TRUE)
    }
    bootstrap_samples <- split(bootstrap_samples, rep(1:S, each = N))
    
    pair_indices <- combn(D, 2, simplify = FALSE)
    if (length(pair_indices) == 0) stop("Error: pair_indices is not populated correctly. Aborting analysis.")
    
    results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
      d1 <- pair[1]
      d2 <- pair[2]
      comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
    
      minsigma_values <- numeric(S)
      maxsigma_values <- numeric(S)
      
      results_inner <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
        Yboot <- Y[, bootstrap_samples[[s]]]
        rWpara <- matrix(NA, D, N)
        for (n in 1:N) {
          rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
        }
        rWpara <- log(rWpara)
    
        taxa1relativesd <- var(rWpara[d1, ])
        taxa2relativesd <- var(rWpara[d2, ])
        relativecorrelation <- cor(rWpara[d1, ], rWpara[d2, ])
        relativecovariance <- cov(rWpara[d1, ], rWpara[d2, ])
        
        res_min <- optim(par = c(0, 0, 0.1), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = objective_function, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperscalevariance))
        res_max <- optim(par = c(0, 0, 0.1), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {-objective_function(params, taxa1relativesd, taxa2relativesd, relativecorrelation)}, method = "L-BFGS-B", lower = c(-rhobound, -rhobound, 0.05), upper = c(rhobound, rhobound, upperscalevariance))
        
        data.frame(
          d1 = d1,
          d2 = d2,
          s = s,
          minsigma_absolute_minimum_covariance = res_min$value,
          maxsigma_absolute_maximum_covariance = -res_max$value,
          minsigma_correlation_relativetaxa1_scale = res_min$par[1],
          minsigma_correlation_relativetaxa2_scale = res_min$par[2],
          minsigma_scale_variance = res_min$par[3],
          maxsigma_correlation_relativetaxa1_scale = res_max$par[1],
          maxsigma_correlation_relativetaxa2_scale = res_max$par[2],
          maxsigma_scale_variance = res_max$par[3],
          taxa1relativesd = taxa1relativesd,
          taxa2relativesd = taxa2relativesd,
          relativecorrelation = relativecorrelation,
          relativecovariance = relativecovariance
        )
      }
      
      minsigma_values <- results_inner$minsigma_absolute_minimum_covariance
      maxsigma_values <- results_inner$maxsigma_absolute_maximum_covariance
      sortedmin <- sort(minsigma_values)
      sortedmax <- sort(maxsigma_values)
      minsigma <- min(sortedmin, na.rm = TRUE)
      maxsigma <- max(sortedmax, na.rm = TRUE)
      cilower <- quantile(sortedmin, probs = 0.025, na.rm = TRUE)
      ciupper <- quantile(sortedmax, probs = 0.975, na.rm = TRUE)
      
      min_index <- which.min(minsigma_values)
      max_index <- which.max(maxsigma_values)
      
      min_rho1 <- results_inner$minsigma_correlation_relativetaxa1_scale[min_index]
      min_rho2 <- results_inner$minsigma_correlation_relativetaxa2_scale[min_index]
      min_x <- results_inner$minsigma_scale_variance[min_index]
      
      max_rho1 <- results_inner$maxsigma_correlation_relativetaxa1_scale[max_index]
      max_rho2 <- results_inner$maxsigma_correlation_relativetaxa2_scale[max_index]
      max_x <- results_inner$maxsigma_scale_variance[max_index]
      
      taxa1relativesd <- results_inner$taxa1relativesd[min_index]
      taxa2relativesd <- results_inner$taxa2relativesd[min_index]
      relativecorrelation <- results_inner$relativecorrelation[min_index]
      relativecovariance <- results_inner$relativecovariance[min_index]
  
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
          minsigma_correlation_relativetaxa1_scale = min_rho1,
          minsigma_correlation_relativetaxa2_scale = min_rho2,
          minsigma_scale_variance = min_x,
          maxsigma_correlation_relativetaxa1_scale = max_rho1,
          maxsigma_correlation_relativetaxa2_scale = max_rho2,
          maxsigma_scale_variance = max_x,
          relative_standard_dev_taxa1 = taxa1relativesd,
          relative_standard_dev_taxa2 = taxa2relativesd,
          relative_correlation = relativecorrelation,
          relative_covariance = relativecovariance,
          stringsAsFactors = FALSE
        )
      )
    }
  
    final_results <- do.call(rbind, lapply(results_list, function(x) x$results))
    final_results <- as.data.frame(final_results, stringsAsFactors = FALSE)
    rownames(final_results) <- NULL
    
    all_inner_results <- do.call(rbind, lapply(results_list, function(x) x$resultsinner))
    all_inner_results <- as.data.frame(all_inner_results)
    all_inner_results$comparison <- paste(rownames(Y)[all_inner_results$d1], rownames(Y)[all_inner_results$d2], sep = ":")

    # End bootstrap timing
    bootstrap_end_time <- Sys.time()
    bootstrap_elapsed_time <- as.numeric(difftime(bootstrap_end_time, bootstrap_start_time, units = "secs"))
    
    # Record the time taken
    bootstrap_times <<- rbind(bootstrap_times, data.frame(S = S, Time = bootstrap_elapsed_time))
    
    return(list(final_results = final_results, all_inner_results = all_inner_results))
  }
  
  # Run bootstrap analysis for different sample sizes and store results
  sample_sizes <- c(100, 500, 1000, 2000, 5000, 10000)
  convergence_results <- lapply(sample_sizes, function(S) {
    run_bootstrap_analysis(S, Y, N, D, alpha, rhobound, upperscalevariance)
  })
  
  # Combine results for plotting
  combined_results <- do.call(rbind, lapply(1:length(sample_sizes), function(i) {
    data.frame(S = sample_sizes[i], convergence_results[[i]]$final_results)
  }))
  
# Modify the first plot to show ranges and add CI, grouped by comparison
  convergencediagnostics_plot <- ggplot(combined_results, aes(x = S)) +
    geom_ribbon(aes(ymin = minsigma_absolute_minimum_covariance, ymax = maxsigma_absolute_maximum_covariance), fill = "lightblue", alpha = 0.5) +
    geom_ribbon(aes(ymin = ninetyfive_ci_lower, ymax = ninetyfive_ci_upper), fill = "lightcoral", alpha = 0.5) +
    facet_wrap(~ comparison, scales = "free_y") +
    labs(title = "Bootstrap Convergence Diagnostics",
         x = "Number of Bootstrap Samples",
         y = "Estimate (with Confidence Interval)") +
    theme_bw() +
    theme(plot.background = element_rect(fill = "white"),
          strip.text = element_text(size = 6))  # Adjust facet label size for readability
  
  # Save the plot as JPG
  ggsave(convergencediagnostics_plot, filename = "convergence_diagnostics_facet.jpg", width = 12, height = 10)
  
  # Print the plot
  print(convergencediagnostics_plot)

  # Plot the time taken for each bootstrap sample size
  bootstrapcomputationtime_plot = ggplot(bootstrap_times, aes(x = S, y = Time)) +
                                    geom_line(color = "blue") +
                                    geom_point(color = "red") +
                                    labs(title = "Time Increase with Bootstrap Sample Size",
                                         x = "Number of Bootstrap Samples",
                                         y = "Time (seconds)") +
                                    theme_bw() +
                                    theme(plot.background = element_rect(fill = "white"))

  # Save the plot as JPG
  ggsave(bootstrapcomputationtime_plot, filename = "bootstrap_timeefficiency.jpg", width = 8, height = 6)
  print(bootstrapcomputationtime_plot)
  
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(list(convergence_results = convergence_results, combined_results = combined_results))
}
