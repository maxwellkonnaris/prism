#' Run Analysis on All Pairwise Taxa
#'
#' This function runs a bootstrapped analysis on the input data matrix \code{Y}.
#'
#' @param Y A matrix of data with observations in columns and variables in rows.
#' @param alpha A numeric vector of priors for the Dirichlet distribution. Defaults to a vector of zeros.
#' @param lowerrhobound A numeric value specifying the lower bound for the \code{lowerrhobound} parameter. Defaults to -1.0.
#' @param upperrhobound A numeric value specifying the upper bound for the \code{upperrhobound} parameter. Defaults to 1.0.
#' @param S An integer specifying the number of bootstrap samples. Increase to reduce Monte Carlo error. Defaults to 1000.
#' @param lowerscalestdev A numeric value specifying the lower bound of the variance of the scale. Defaults to 0.49.
#' @param upperscalestdev A numeric value specifying the upper bound of the variance of the scale. Defaults to 0.51.
#' @return A list of dataframes containing the results of the analysis including estimated 95% confidence intervals, minimum and maximum values for estimated absolute covariance, and parameters. \code{final_results} contains the data intended for forest_plot() and \code{all_inner_results} contains the data intended for sigmaplot().
#' @import progress
#' @import progressr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MCMCpack
#' @import stats
#' @import utils
#' @import nloptr
#' @examples
#' # Example usage (You could also use the simulation function provided to generate sample data):
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10)
#' results <- estimate_covariance(Y)
#' @export

estimate_covariance <- function(Y, alpha = 0.5, lowerrhobound = -1.0, upperrhobound = 1.0, S = 1000, lowerscalestdev = .49, upperscalestdev = .51, algorithm="COBYLA", scalestep=0.005) {

  ## COMPUTATIONAL TIME -------------------------------------------------------------------------------------------------------------------------
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
  
  ## END COMPUTATIONAL TIME SETUP ---------------------------------------------------------------------------------------------------------------
  
  ## SETUP --------------------------------------------------------------------------------------------------------------------------------------
  # Check if Y is a matrix, dataframe, or tibble, and has appropriate dimensions
  if (!(is.matrix(Y) || is.data.frame(Y) || inherits(Y, "tbl_df")) || nrow(Y) < 2 || ncol(Y) < 2) {
    stop("Y must be a matrix, dataframe, or tibble with at least 2 rows and 2 columns.")
  }

  # Get the number of columns (N) and rows (D) in the input matrix Y
  N <- ncol(Y)
  D <- nrow(Y)
  
  # Print priors
  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Rho bounds:\n")
  print(paste0(lowerrhobound,":",upperrhobound))
  cat("Scale standard deviation bounds:\n")
  print(paste0(lowerscalestdev,":",upperscalestdev))
  cat("Dimensions of supplied matrix:\n")
  print(dim(Y))
  cat("Bootstrap sample size (S):\n")
  print(S)

  # Definine the initial parameters for the Optimization                                      
  initialparameters = c(((upperrhobound+lowerrhobound)/2), ((upperrhobound+lowerrhobound)/2), ((upperscalestdev+lowerscalestdev)/2))

  if (algorithm == "GRID_SEARCH") {

    # Define parameter steps
    rho1 <- seq(lowerrhobound, upperrhobound, by=0.05)
    rho2 <- seq(lowerrhobound, upperrhobound, by=0.05)
    scalestdevstep <- seq(lowerscalestdev, upperscalestdev, scalestep)
    iterations <- length(rho1) * length(rho2) * length(scalestdevstep)
    
    # Create the grid of parameters
    pars <- expand.grid(rho1, rho2, scalestdevstep)

    print(pars)
    
  }
  
  # Initialize a results matrix to store the results for each pair
  results <- matrix(list(), D, D)
  ## END SETUP -----------------------------------------------------------------------------------------------------------------------------------

  ## OPTIMIZATION FUNCTIONS ----------------------------------------------------------------------------------------------------------------------
  # Objective function used to optimize the covariance whether minimum or maximum
  objective_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalestdev <- params[3]
    
    sigma <- relativecovariance + scalestdev * taxa1relativesd * taxa1scalecorrelation + scalestdev * taxa2relativesd * taxa2scalecorrelation + scalestdev^2

    return(sigma)
  }

  # Gradient function for covariance
  gradient_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalestdev <- params[3]
    
    # Calculate partial derivatives
    grad_rho1 <- taxa1relativesd * scalestdev
    grad_rho2 <- taxa2relativesd * scalestdev
    grad_x <- taxa1relativesd * taxa1scalecorrelation + taxa2relativesd * taxa2scalecorrelation + 2 * scalestdev
    
    # Return the gradient as a vector
    return(c(grad_rho1, grad_rho2, grad_x))
  }

  # Symmetric Positive Semi-Definite case constraint
  constraint_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {

    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalestdev <- params[3]  
    
    term1 <- taxa1relativesd*taxa1scalecorrelation + taxa2relativesd*taxa2scalecorrelation
    term2 <- sqrt(2) * sqrt(((taxa1relativesd^2 * taxa1scalecorrelation^2) + (taxa2relativesd^2 * taxa2scalecorrelation^2)))
                            
    term3 <- (1 / (2 * scalestdev)) * ((taxa1relativesd^2 + taxa2relativesd^2) 
            - sqrt((taxa1relativesd^2 - taxa2relativesd^2)^2 + 4 * relativecovariance^2) 
            + 4 * scalestdev^2)
    
    g_1 <- term1 - term2 + term3
    
    return(g_1)
  }

  constraint_gradient_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    taxa1scalecorrelation <- params[1]
    taxa2scalecorrelation <- params[2]
    scalestdev <- params[3]
    
    # Add a small epsilon to avoid division by zero
    epsilon <- 1e-8
    
    denominator <- sqrt(taxa1relativesd^2 * taxa1scalecorrelation^2 + taxa2relativesd^2 * taxa2scalecorrelation^2 + epsilon)
    
    grad_taxa1scalecorrelation <- taxa1relativesd - sqrt(2) * (taxa1relativesd^2 * taxa1scalecorrelation) / denominator
    grad_taxa2scalecorrelation <- taxa2relativesd - sqrt(2) * (taxa2relativesd^2 * taxa2scalecorrelation) / denominator
    
    term_to_simplify <- (taxa1relativesd^2 + taxa2relativesd^2) - sqrt((taxa1relativesd^2 - taxa2relativesd^2)^2 + 4 * relativecovariance^2)
    grad_scalestdev <- -1 / (2 * (scalestdev + epsilon)^2) * term_to_simplify + 2 * scalestdev
    
    return(c(grad_taxa1scalecorrelation, grad_taxa2scalecorrelation, grad_scalestdev))
  }

  # Wrapper function for the objective function
  objective_function_wrapper <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    objective_function(params, taxa1relativesd, taxa2relativesd, relativecovariance)
  }
  
  # Wrapper function for the constraint function
  constraint_function_wrapper <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    constraint_function(params, taxa1relativesd, taxa2relativesd, relativecovariance)
  }
  
  # Wrapper function for the gradient function
  gradient_function_wrapper <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    gradient_function(params, taxa1relativesd, taxa2relativesd, relativecovariance)
  }
  
  # Wrapper function for the constraint gradient function
  constraint_gradient_function_wrapper <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
    constraint_gradient_function(params, taxa1relativesd, taxa2relativesd, relativecovariance)
  }
  ## END OPTIMIZATION FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------

  ## CLUSTER RESOURCES --------------------------------------------------------------------------------------------------------------------------------------
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
  ## END CLUSTER RESOURCES ----------------------------------------------------------------------------------------------------------------------------------

  ## PROGRESS BARS ------------------------------------------------------------------------------------------------------------------------------------------
  # Define global handlers for progress bars
  handlers(global = TRUE)
  
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
  ## END PROGRESS BARS -------------------------------------------------------------------------------------------------------------------------------------
  
  ## BOOTSTRAP PRECOMPUTE ----------------------------------------------------------------------------------------------------------------------------------
  # Parallelize bootstrap precomputation
  bootstrap_samples <- foreach(s = 1:S, .combine = 'c', .options.snow = opts_precomp) %dopar% {
    sample(1:N, replace = TRUE)
  }
  
  # Reshape bootstrap_samples into a list of vectors
  bootstrap_samples <- split(bootstrap_samples, rep(1:S, each = N))
  ## END BOOTSTRAP PRECOMPUTE ------------------------------------------------------------------------------------------------------------------------------
  
  ## SIGMA ESTIMATION --------------------------------------------------------------------------------------------------------------------------------------
  # Generate all pairs of indices and add diagonal pairs
  pair_indices <- combn(D, 2, simplify = FALSE)

  # Check if pair_indices is populated correctly
  if (length(pair_indices) == 0) {
    stop("Error: pair_indices is not populated correctly. Aborting analysis.")
  }

  # Run the analysis
  cat("Running sigma estimation")
  results_list <- tryCatch({ foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
      d1 <- pair[1]
      d2 <- pair[2]
      comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
    
      minsigma_values <- numeric(S)
      maxsigma_values <- numeric(S)
      
      # Use parallel foreach for the inner loop
      results_inner <- tryCatch({ foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack'), .options.snow = opts) %dopar% {
          Yboot <- Y[, bootstrap_samples[[s]]]
          
          rWpara <- matrix(NA, D, N)
          for (n in 1:N) {
            rWpara[, n] <- MCMCpack::rdirichlet(1, Yboot[, n] + alpha)
          }

          # log transform relative abundances
          rWpara <- log(rWpara)
          rWpara = rWpara[c(d1,d2), ]
          
          taxa1relativesd <- sd(rWpara[1, ])
          taxa2relativesd <- sd(rWpara[2, ])
          relativecorrelation <- cor(rWpara[1, ], rWpara[2, ])
          relativecovariance <- cov(rWpara[1, ], rWpara[2, ])

          res_min <- NULL
          res_max <- NULL
        
          result <- switch(algorithm,
  
            "COBYLA" = {
              # Find the minimum sigma using nloptr with COBYLA
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
              )
              
              # Find the maximum sigma by negating the objective function using COBYLA
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
              )
              
              list(res_min = res_min, res_max = res_max)
            },
          
            "MMA" = {
              # Find the minimum sigma using nloptr with MMA
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list("algorithm"="NLOPT_LD_MMA", "maxeval" = 1000000, "ftol_rel" = 1e-5)
              )
              
              # Find the maximum sigma by negating the objective function using MMA
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list("algorithm"="NLOPT_LD_MMA", "maxeval" = 1000000, "ftol_rel" = 1e-5)
              )
              
              list(res_min = res_min, res_max = res_max)
            },
          
            "AUGLAG_COBYLA" = {
              # Find the minimum sigma using AUGLAG with COBYLA as the inner algorithm
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LN_COBYLA",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              # Find the maximum sigma using AUGLAG with COBYLA as the inner algorithm
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LN_COBYLA",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              list(res_min = res_min, res_max = res_max)
            },
          
            "AUGLAG_MMA" = {
              # Find the minimum sigma using AUGLAG with MMA as the inner algorithm
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              # Find the maximum sigma using AUGLAG with MMA as the inner algorithm
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              list(res_min = res_min, res_max = res_max)
            },
          
            "AUGLAG_SLSQP" = {
              # Find the minimum sigma using AUGLAG with SLSQP as the inner algorithm
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_SLSQP",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              # Find the maximum sigma using AUGLAG with SLSQP as the inner algorithm
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_SLSQP",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              list(res_min = res_min, res_max = res_max)
            },
          
            "AUGLAG_LBFGS" = {
              # Find the minimum sigma using AUGLAG with LBFGS as the inner algorithm
              res_min <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_LBFGS",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              # Find the maximum sigma using AUGLAG with LBFGS as the inner algorithm
              res_max <- nloptr(
                x0 = initialparameters,
                eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                lb = c(lowerrhobound, lowerrhobound, lowerscalestdev),
                ub = c(upperrhobound, upperrhobound, upperscalestdev),
                opts = list(
                  "algorithm" = "NLOPT_LD_AUGLAG",
                  "local_opts" = list(
                    "algorithm" = "NLOPT_LD_LBFGS",
                    "xtol_rel" = 1e-5,
                    "maxeval" = 1000000
                  ),
                  "maxeval" = 1000000,
                  "ftol_rel" = 1e-5
                )
              )
              
              list(res_min = res_min, res_max = res_max)
            },

            "GRID_SEARCH" = {

              # Step 1: Calculate sigma with tryCatch
              result_sigma <- tryCatch({
                pars %>%
                  mutate(sigma = relativecovariance + scalestdevstep * taxa1relativesd * rho1 +
                           scalestdevstep * taxa2relativesd * rho2 + scalestdevstep^2)
              }, error = function(e) {
                message("Error occurred during sigma calculation: ", e$message)
                stop("Stopping execution due to error in sigma calculation.")
              })
              
              # Step 2: Apply rowwise() and compute SPSD with tryCatch
              result_spsd <- tryCatch({
                result_sigma %>%
                  rowwise() %>%
                  mutate(SPSD = constraint_function(params = c(rho1, rho2, scalestdevstep), 
                                                    taxa1relativesd, 
                                                    taxa2relativesd, 
                                                    relativecovariance))
              }, error = function(e) {
                message("Error occurred during SPSD calculation: ", e$message)
                stop("Stopping execution due to error in SPSD calculation.")
              })
              
              # Step 3: Filter rows based on SPSD with tryCatch
              result_filtered <- tryCatch({
                result_spsd %>%
                  filter(SPSD >= 0) %>%
                  ungroup()
              }, error = function(e) {
                message("Error occurred during SPSD filtering: ", e$message)
                stop("Stopping execution due to error in SPSD filtering.")
              })
              
              # Second tryCatch: Find the row corresponding to minimum sigma
              res_min <- tryCatch({
                if (nrow(rpars) > 0) {
                  rpars %>%
                    slice(which.min(sigma)) %>%
                    mutate(objective = sigma,
                           solution = list(c(rho1, rho2, scalestdevstep)),  # Use list for multi-dimensional values
                           message = "GRIDSEARCH",
                           status = "GRIDSEARCH",
                           iterations = iterations)
                } else {
                  stop("No rows remaining after filtering based on SPSD constraints.")
                }
              }, error = function(e) {
                message("Error occurred while finding the row with minimum sigma: ", e$message)
                stop("Stopping execution due to error in finding minimum sigma.")
              })
              
              # Third tryCatch: Find the row corresponding to maximum sigma
              res_max <- tryCatch({
                if (nrow(rpars) > 0) {
                  rpars %>%
                    slice(which.max(sigma)) %>%
                    mutate(objective = -sigma,
                           solution = list(c(rho1, rho2, scalestdevstep)),  # Use list for multi-dimensional values
                           message = "GRIDSEARCH",
                           status = "GRIDSEARCH",
                           iterations = iterations)
                } else {
                  stop("No rows remaining after filtering based on SPSD constraints.")
                }
              }, error = function(e) {
                message("Error occurred while finding the row with maximum sigma: ", e$message)
                stop("Stopping execution due to error in finding maximum sigma.")
              })
              
              list(res_min = res_min, res_max = res_max)
            },
          
            stop("Invalid algorithm selected") # Default case if no match is found
          )

          res_min = result$res_min
          res_max = result$res_max

          if (is.null(res_min) || is.null(res_max)) {
          stop("Optimization failed for task.")
          }

          data.frame(
            d1 = d1,
            d2 = d2,
            s = s,
            minsigma_absolute_minimum_covariance = res_min$objective,
            minsigma_correlation_relativetaxa1_scale = res_min$solution[1],
            minsigma_correlation_relativetaxa2_scale = res_min$solution[2],
            minsigma_scale_variance = res_min$solution[3],
            minsigma_message = res_min$message,
            minsigma_status = res_min$status,
            minsigma_iterations = res_min$iterations,
            maxsigma_absolute_maximum_covariance = -res_max$objective,
            maxsigma_correlation_relativetaxa1_scale = res_max$solution[1],
            maxsigma_correlation_relativetaxa2_scale = res_max$solution[2],
            maxsigma_scale_variance = res_max$solution[3],
            maxsigma_message = res_max$message,
            maxsigma_status = res_max$status,
            maxsigma_iterations = res_max$iterations,
            taxa1relativesd = taxa1relativesd,
            taxa2relativesd = taxa2relativesd,
            relativecorrelation = relativecorrelation,
            relativecovariance = relativecovariance#,
            #variance_proportionality_taxa1 = variance_proportionality_d1,
            #variance_proportionality_taxa2 = variance_proportionality_d2
          )    
        }
      }, error = function(e) {
            # Stop the execution and show the error message
            stop("Error in optimization: ", e$message)
      })
      
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
  }, error = function(e) {
        # Handle the error and stop all parallel tasks
        stopCluster(cl)
        stop("Error in parallel computation: ", e$message)
  })
  
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
#' @param lowerrhobound A numeric value specifying the lower bound for the \code{lowerrhobound} parameter. Defaults to -1.0.
#' @param upperrhobound A numeric value specifying the upper bound for the \code{upperrhobound} parameter. Defaults to 0.6.
#' @param S A vector of numeric integers specifying the number of bootstrap samples to be iterated across. Increase to reduce Monte Carlo error. Defaults to \code{c(100, 500, 1000, 2000, 5000, 10000)}.
#' @param lowerscalestdev A numeric value specifying the lower bound of the variance of the scale. Defaults to 0.475.
#' @param upperscalestdev A numeric value specifying the upper bound of the variance of the scale. Defaults to 0.525.
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
estimate_covariance_convergence <- function(Y, S=c(100, 500, 1000, 2000, 5000, 10000), alpha = 0.5, lowerrhobound = -1.0, upperrhobound = 0.6, lowerscalestdev = .475, upperscalestdev = .525) {
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

  # Check if Y is a matrix, dataframe, or tibble, and has appropriate dimensions
  if (!(is.matrix(Y) || is.data.frame(Y) || inherits(Y, "tbl_df")) || nrow(Y) < 2 || ncol(Y) < 2) {
    stop("Y must be a matrix, dataframe, or tibble with at least 2 rows and 2 columns.")
  }

  # Get the number of columns (N) and rows (D) in the input matrix Y
  N <- ncol(Y)
  D <- nrow(Y)

  # Print priors
  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Rho bounds:\n")
  print(paste0(lowerrhobound,":",upperrhobound))
  cat("Scale standard deviation bounds:\n")
  print(paste0(lowerscalestdev,":",upperscalestdev))
  cat("Dimensions of supplied matrix:\n")
  print(dim(Y))
  cat("Bootstrap sample size (S):\n")
  print(S)

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

  # Function to run bootstrap analysis for a given number of samples
  run_bootstrap_analysis <- function(S, Y, N, D, alpha, lowerrhobound, upperrhobound, lowerscalestdev, upperscalestdev) {

    # Define the objective function used in the optimization
    objective_function <- function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {
      taxa1scalecorrelation <- params[1]
      taxa2scalecorrelation <- params[2]
      scalestdev <- params[3]
      sigma <- taxa1relativesd * taxa2relativesd * relativecorrelation + taxa1relativesd * scalestdev * taxa1scalecorrelation + taxa2relativesd * scalestdev * taxa2scalecorrelation + scalestdev^2
      return(sigma)
    }

    # Initialize a results matrix to store the results for each pair
    results <- matrix(list(), D, D)
    
    # Run the analysis
    cat("Running sigma estimation for bootstrap sample size (S):",S,"\n")
    
    # Start timing
    bootstrap_start_time <- Sys.time()

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
    
        taxa1relativesd <- sd(rWpara[d1, ])
        taxa2relativesd <- sd(rWpara[d2, ])
        relativecorrelation <- cor(rWpara[d1, ], rWpara[d2, ])
        relativecovariance <- cov(rWpara[d1, ], rWpara[d2, ])
        
        res_min <- optim(par = c(0, 0, 0.15), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = objective_function, method = "L-BFGS-B", lower = c(lowerrhobound, lowerrhobound, lowerscalestdev), upper = c(upperrhobound, upperrhobound, upperscalestdev))
        res_max <- optim(par = c(0, 0, 0.15), taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecorrelation = relativecorrelation, fn = function(params, taxa1relativesd, taxa2relativesd, relativecorrelation) {-objective_function(params, taxa1relativesd, taxa2relativesd, relativecorrelation)}, method = "L-BFGS-B", lower = c(lowerrhobound, lowerrhobound, lowerscalestdev), upper = c(upperrhobound, upperrhobound, upperscalestdev))
        
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
  convergence_results <- lapply(S, function(SS) {
    run_bootstrap_analysis(SS, Y, N, D, alpha, lowerrhobound, upperrhobound, lowerscalestdev, upperscalestdev)
  })
  
  # Combine results for plotting
  combined_results <- do.call(rbind, lapply(1:length(S), function(i) {
    data.frame(S = S[i], convergence_results[[i]]$final_results)
  }))

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
  ggsave(bootstrapcomputationtime_plot, filename = "bootstrap_timeefficiency.jpg", width = 5, height = 5)
  
  # Modify the first plot to show ranges and add CI, grouped by comparison with flipped axes and categorical bootstrap samples
  convergencediagnostics_plot <- ggplot(combined_results, aes(y = factor(S))) +
    geom_errorbarh(aes(xmin = minsigma_absolute_minimum_covariance, xmax = maxsigma_absolute_maximum_covariance, color = "Covariance Range"), height = 0.2) +
    geom_errorbarh(aes(xmin = ninetyfive_ci_lower, xmax = ninetyfive_ci_upper, color = ifelse(ninetyfive_ci_lower <= 0 & ninetyfive_ci_upper >= 0, "95% CI Covers Zero", "95% CI Does Not Covers Zero")), height = 0.2) +
    facet_wrap(~ comparison, scales = "free_x") +
    labs(title = "Bootstrap Convergence Diagnostics",
         y = "Number of Bootstrap Samples",
         x = "Covariance Set Estimate (with Confidence Interval)") +
    scale_color_manual(values = c("Covariance Range" = "black", 
                                  "95% CI Covers Zero" = "red",
                                  "95% CI Does Not Cover Zero" = "green"),
                       name = "Legend") +
    theme_bw() +
    theme(plot.background = element_rect(fill = "white"),
          strip.text = element_text(size = 6))  # Adjust facet label size for readability
  
  # Save the plot as JPG
  ggsave(convergencediagnostics_plot, filename = "convergence_diagnostics_facet.jpg", width = 30, height = 20)
  
  # Print the plot
  print(convergencediagnostics_plot)
  
  # Calculate and print the total elapsed time
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  formatted_time <- format_elapsed_time(elapsed_time)
  print(paste("Total time taken:", formatted_time))
  
  return(list(convergence_results = convergence_results, combined_results = combined_results))
}
