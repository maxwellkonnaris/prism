#' Run Analysis on All Pairwise Taxa
#'
#' This function runs a bootstrapped analysis on the input data matrix \code{Y}.
#'
#' @param Y A matrix of data with observations in columns and variables in rows.
#' @param alpha A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros.
#' @param rhobound A numeric value specifying the bound for the \code{rho1} and \code{rho2} parameters. Defaults to 0.8.
#' @param S An integer specifying the number of bootstrap samples. Defaults to 1000.
#' @return A data frame containing the results of the analysis including estimated 95% confidence intervals, minimum and maximum values for estimated covariance, and finite sample covariances.
#' @import progress
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MCMCpack
#' @import stats
#' @import profvis
#' @examples
#' # Example usage (You could also use the simulation function provided to generate sample data):
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10)
#' results <- run_analysis(Y)
#' @export
run_gridanalysis <- function(Y, alpha = rep(0, nrow(Y)), rhobound = 0.8, S = 1000) {
  start_time <- Sys.time()

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

  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (!foreach::getDoParRegistered()) {
    stop("Parallel backend is not registered.")
  } else {
    cat("Number of workers/cpus: ", foreach::getDoParWorkers(), "\n")
  }
  
  total_pairs <- D * (D + 1) / 2
  pb_precomp <- progress::progress_bar$new(total = S, format = "  Precomputing bootstrap samples [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  pb <- progress::progress_bar$new(total = total_pairs, format = " Generating Sigmas [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  
  progress_precomp <- function(n) {
    pb_precomp$tick()
  }
  
  progress <- function(n) {
    pb$tick()
  }
  
  opts_precomp <- list(progress = progress_precomp)
  opts <- list(progress = progress)

  bootstrap_samples <- foreach(s = 1:S, .combine = 'c', .options.snow = opts_precomp) %dopar% {
    sample(1:N, replace = TRUE)
  }
  
  bootstrap_samples <- split(bootstrap_samples, rep(1:S, each = N))
  
  pair_indices <- combn(D, 2, simplify = FALSE)
  pair_indices <- c(pair_indices, lapply(1:D, function(x) c(x, x)))

  cat("Running sigma estimation")
  
  results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .combine = 'list', .options.snow = opts) %dopar% {
    d1 <- pair[1]
    d2 <- pair[2]
  
    sigma_values <- foreach(s = 1:S, .combine = 'c', .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
      Yboot <- Y[, bootstrap_samples[[s]]]
      sigma_vals <- apply(pars, 1, function(params) {
        objective_function(params, Yboot, d1, d2, alpha)
      })
      return(c(min(sigma_vals), max(sigma_vals)))
    }
  
    minsigma_values <- sigma_values[seq(1, length(sigma_values), 2)]
    maxsigma_values <- sigma_values[seq(2, length(sigma_values), 2)]
    
    sortedmin <- sort(minsigma_values)
    sortedmax <- sort(maxsigma_values)
    
    minsigma <- min(sortedmin)
    maxsigma <- max(sortedmax)
    cilower <- quantile(sortedmin, probs = 0.025)
    ciupper <- quantile(sortedmax, probs = 0.975)
  
    finitesamplecovariance <- stats::cov(log(Y)[d1,], log(Y)[d2,])
    
    list(d1 = d1, d2 = d2, cilower = cilower, ciupper = ciupper, minsigma = minsigma, maxsigma = maxsigma, finitesamplecovariance = finitesamplecovariance)
  }
  
  # Convert list of lists into a data frame
  final_results <- do.call(rbind, lapply(results_list, function(x) data.frame(matrix(unlist(x), ncol = 7, byrow = TRUE))))
  colnames(final_results) <- c("d1", "d2", "cilower", "ciupper", "minsigma", "maxsigma", "finitesamplecovariance")
  
  # Format results into a data frame
  formatted_results <- data.frame(
    comparison = apply(final_results, 1, function(row) paste(rownames(Y)[row["d1"]], ":", rownames(Y)[row["d2"]], sep = "")),
    cilower = final_results$cilower,
    ciupper = final_results$ciupper,
    minsigma = final_results$minsigma,
    maxsigma = final_results$maxsigma,
    finitesamplecovariance = final_results$finitesamplecovariance,
    stringsAsFactors = FALSE
  )
  
  # Remove row names
  rownames(formatted_results) <- NULL
  
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(formatted_results)
}
