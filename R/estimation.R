#' Run Bootstrapped Analysis on Pairwise Taxa
#'
#' This function performs a bootstrapped analysis on the input matrix \code{Y}, estimating 
#' covariance, correlation, and scale for all pairwise combinations of taxa. It supports 
#' various optimization algorithms to estimate the minimum and maximum covariance.
#'
#' @param Y A matrix of observed counts, where rows represent variables (e.g., taxa) and columns represent observations (samples). The matrix should have at least two rows and two columns.
#' 
#' @param alpha A numeric vector of Dirichlet priors with the same length as the number of rows in \code{Y}. If \code{alpha} is a scalar, it will be replicated for each row. Defaults to \code{0.5}.
#' 
#' @param lowerrhobound A numeric vector of lower bounds for correlation parameters (\eqn{\rho}). Each element specifies the lower bound for a row of \code{Y}. If a scalar is provided, it will be replicated for each row. Defaults to \code{-1.0}.
#' 
#' @param upperrhobound A numeric vector of upper bounds for correlation parameters (\eqn{\rho}). Each element specifies the upper bound for a row of \code{Y}. If a scalar is provided, it will be replicated for each row. Defaults to \code{1.0}.
#' 
#' @param S An integer specifying the number of bootstrap samples. Larger values reduce Monte Carlo error but increase computation time. Defaults to \code{1000}.
#' 
#' @param lowerscalestdev A numeric value specifying the lower bound of the standard deviation of the scale. Defaults to \code{0.49}.
#' 
#' @param upperscalestdev A numeric value specifying the upper bound of the standard deviation of the scale. Defaults to \code{0.51}.
#' 
#' @param algorithm A character string specifying the optimization algorithm to be used. 
#' Can be one of \code{"COBYLA"}, \code{"MMA"}, \code{"SLSQP"}, \code{"AUGLAG_COBYLA"}, \code{"AUGLAG_MMA"}, \code{"AUGLAG_SLSQP"}, or \code{"GRID_SEARCH"}. Defaults to \code{"COBYLA"}.
#' 
#' @param outputdirectory A character string specifying the directory to save results for grid search. If \code{NULL}, the current working directory is used. Defaults to \code{NULL}.
#'
#' @return A list with two elements:
#' \item{\code{final_results}}{A dataframe containing the minimum and maximum covariance estimates for each pair of taxa along with confidence intervals and other relevant metrics.}
#' \item{\code{all_inner_results}}{A dataframe containing the full set of results from the bootstrapped analysis, including correlation, covariance, and scale estimates for each bootstrap sample.}
#'
#' @details
#' The function \code{estimate_covariance} performs a comprehensive analysis of pairwise relationships between taxa in the input matrix \code{Y}. For each pair of taxa, it:
#' \enumerate{
#'   \item Generates \code{S} bootstrap samples.
#'   \item Estimates covariance and correlation parameters using the specified optimization algorithm.
#'   \item Computes 95\% confidence intervals for the estimated covariances.
#'   \item Aggregates results across all pairs and bootstrap samples.
#' }
#'
#' Supported optimization algorithms include:
#' \describe{
#'   \item{\code{"COBYLA"}}{Constrained Optimization BY Linear Approximations.}
#'   \item{\code{"MMA"}}{Method of Moving Asymptotes.}
#'   \item{\code{"AUGLAG_COBYLA"}}{Augmented Lagrangian with COBYLA as the inner optimizer.}
#'   \item{\code{"AUGLAG_MMA"}}{Augmented Lagrangian with MMA as the inner optimizer.}
#'   \item{\code{"GRID_SEARCH"}}{Exhaustive grid search over specified parameter ranges.}
#' }
#'
#' The \code{final_results} dataframe is intended for visualization functions such as \code{forest_plot()}, while \code{all_inner_results} is suitable for detailed analysis and plotting with functions like \code{sigmaplot()}.
#'
#' @section Input and Output Structures:
#'
#' \strong{Input Matrix (\code{Y}):}
#'
#' The input matrix \code{Y} should be structured with taxa as rows and samples as columns. Each entry represents the observed count for a specific taxon in a given sample.
#'
#' \preformatted{
#'      Sample1 Sample2 Sample3 ... SampleN
#' Taxon1    x11     x12     x13        x1N
#' Taxon2    x21     x22     x23        x2N
#' ...       ...     ...     ...         ...
#' TaxonD    xD1     xD2     xD3        xDN
#' }
#'
#' \strong{Output Dataframes:}
#'
#' \emph{final_results}:
#'
#' Contains summary statistics for each pair of taxa, including covariance estimates and confidence intervals.
#'
#' \preformatted{
#'   comparison    taxa1    taxa2 proportion_intervals_dontcoverzero ninetyfive_ci_lower ninetyfive_ci_upper minsigma_absolute_minimum_covariance maxsigma_absolute_maximum_covariance ...
#'   Taxon1:Taxon2 Taxon1    Taxon2                            0.95                -0.8                0.7                            -0.85                             0.75 ...
#'   Taxon1:Taxon3 Taxon1    Taxon3                            0.90                -0.6                0.6                            -0.65                             0.55 ...
#'   ...            ...       ...                              ...                  ...                  ...                              ...                               ...
#' }
#'
#' \emph{all_inner_results}:
#'
#' Contains detailed results for each bootstrap sample and pair of taxa.
#'
#' \preformatted{
#'   d1 d2 s minsigma_absolute_minimum_covariance minsigma_correlation_relativetaxa1_scale minsigma_correlation_relativetaxa2_scale minsigma_scale_sd maxsigma_absolute_maximum_covariance ...
#'    1  2 1                            -0.80                               0.10                               -0.05                      0.50                             0.75 ...
#'    1  2 2                            -0.85                               0.12                               -0.04                      0.51                             0.73 ...
#'    ... ... ...                             ...                                 ...                                 ...                       ...                               ...
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic Usage with Default Parameters
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10, ncol = 100)
#' rownames(Y) <- paste0("Taxon", 1:10)
#' colnames(Y) <- paste0("Sample", 1:100)
#' 
#' results <- estimate_covariance(Y)
#' 
#' # View the final results
#' head(results$final_results)
#' 
#' # Example Output:
#' \dontrun{
#'   comparison   taxa1  taxa2 proportion_intervals_dontcoverzero ninetyfive_ci_lower ninetyfive_ci_upper minsigma_absolute_minimum_covariance maxsigma_absolute_maximum_covariance
#'   Taxon1:Taxon2 Taxon1 Taxon2                           0.95                 -0.8                  0.7                             -0.85                              0.75
#'   Taxon1:Taxon3 Taxon1 Taxon3                           0.90                 -0.6                  0.6                             -0.65                              0.55
#'   ...            ...     ...                             ...                   ...                    ...                               ...                                ...
#' }
#'
#' # Example 2: Using a Different Optimization Algorithm (MMA)
#' results_mma <- estimate_covariance(Y, algorithm = "MMA")
#' 
#' # Example 3: Custom Priors and Bounds
#' alpha_custom <- rep(0.5, 10)
#' lowerrhobound_custom <- rep(-0.9, 10)
#' upperrhobound_custom <- rep(0.9, 10)
#' results_custom <- estimate_covariance(Y, alpha = alpha_custom, 
#'                                       lowerrhobound = lowerrhobound_custom, 
#'                                       upperrhobound = upperrhobound_custom)
#' 
#' # Example 4: Specifying Output Directory for Grid Search
#' results_grid <- estimate_covariance(Y, algorithm = "GRID_SEARCH", 
#'                                     outputdirectory = "grid_search_results/")
#' }
#'
#' @import progress
#' @import progressr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MCMCpack
#' @import stats
#' @import utils
#' @import nloptr
#' @import filelock
#' @import ggridges
#' @export
estimate_covariance <- function(Y, alpha = 0.5, uncertaintydistribution = "Multinomial Dirichlet", externalscalemeasurements = NULL, lowerrhobound = rep(-1.0, nrow(Y)), upperrhobound = rep(1.0, nrow(Y)), S = 1000, lowerscalestdev = 0.450, upperscalestdev = 0.650, algorithm = "GRID_SEARCH", outputdirectory = NULL) {
  
  ## COMPUTATIONAL TIME -------------------------------------------------------------------------------------------------------------------------------------
  start_time <- Sys.time()
  ## END COMPUTATIONAL TIME SETUP ---------------------------------------------------------------------------------------------------------------------------
  ## SETUP --------------------------------------------------------------------------------------------------------------------------------------------------
  # Check if Y is a matrix, dataframe, or tibble, and has appropriate dimensions
  if (!(is.matrix(Y) || is.data.frame(Y) || inherits(Y, "tbl_df")) || nrow(Y) < 2 || ncol(Y) < 2) {
    stop("Y must be a matrix, dataframe, or tibble with at least 2 rows and 2 columns.")
  }

  # Get the number of columns (N) and rows (D) in the input matrix Y
  N <- ncol(Y)
  D <- nrow(Y)

  cat("Priors used for the analysis:\n")
  cat("Alpha:\n")
  print(alpha)
  cat("Dimensions of supplied Y matrix:\n")
  print(paste0("Number of Taxa: ", D))
  print(paste0("Number of Samples: ", N))
  cat("Approximating relative counts with the:\n") 
  print(paste0(uncertaintydistribution))  
  cat("Bootstrap sample size (S):\n")
  print(S)
  cat("Algorithm selected:\n")
  print(algorithm)

  if (!is.null(externalscalemeasurements) && is.matrix(externalscalemeasurements) && ncol(Y) == nrow(externalscalemeasurements)) {
      cat("External scale measurements were provided:\n")
      cat("Dimensions of supplied external scale measurements matrix:\n")
      replicates <- ncol(externalscalemeasurements)
      sampletotals <- nrow(externalscalemeasurements)
      print(paste0("Number of Sample-scale Measurement Pairs: ", sampletotals))
      print(paste0("Number of Replicates: ", replicates))
      cat("Estimating Rho bounds and scale SD from the external scale measurements:\n")
  } else {
      cat("Using default Rho bounds:\n")
      print(paste0(lowerrhobound, ":", upperrhobound))
      cat("Using default Scale standard deviation bounds:\n")
      print(paste0(lowerscalestdev, ":", upperscalestdev))
  }

  ## END SETUP ----------------------------------------------------------------------------------------------------------------------------------------------
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
  progressr::handlers(global = TRUE)

  # Define total comparisons
  total_pairs <- D * (D - 1) / 2

  # Create progress bars
  pb <- progress::progress_bar$new(total = total_pairs, format = "Generating Sigmas [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  pb_rhosd <- progress::progress_bar$new(total = D*S, format = "Generating Rho and SD [:bar] :percent in :elapsed | eta: :eta", clear = FALSE, width = 100)
  
  # Function to update progress bar
  progress <- function(n) {
    pb$tick()
  }
  progress_rhosd <- function(n) {
    pb_rhosd$tick()
  }
  # Options for foreach to include progress updates
  opts <- list(progress = progress)
  opts_rhosd <- list(progress = progress_rhosd)
  ## END PROGRESS BARS -------------------------------------------------------------------------------------------------------------------------------------
  ## ACCOUNTING FOR UNCERTAINTY IN FINITE SAMPLING ---------------------------------------------------------------------------------------------------------
  ## calculate bootstrap resampling -- accounting for finite sampling
  boostrap_samples <- matrix(NA, N, S)
  for (s in 1:S) {
    boostrap_samples[,s] <- sample(1:N, replace=TRUE)
  }
  ## END ACCOUNTING FOR UNCERTAINTY IN FINITE SAMPLING -----------------------------------------------------------------------------------------------------
  ## ACCOUNTING FOR UNCERTAINTY IN OBSERVED RELATIVE ABUNDANCES --------------------------------------------------------------------------------------------
  # Dimensions: (n_taxa, n_samples, n_iter) -- populate matrix of NAs
  rWparaoriginal <- array(NA, dim = c(D, N, S))
  
  # calculate posterior samples -- accounting for uncertainty in the observed relative abundances
  if (uncertaintydistribution == "Multinomial Dirichlet") {
    # generate S Dirichlet samples for each sample (column)
    for (n in 1:N) {
        rWparaoriginal[,n,] <- rdirichlet(S, Y[,n] + alpha) 
    } 
  } else if (uncertaintydistribution == "Multinomial Logistic Normal") {
    # generate S Multinomial logistic Normal posterior samples for each sample (column) using fido
    otu_table = otu_table(Y, taxa_are_rows = TRUE)
    otu_table = otu_table + alpha
    X <- matrix(1, ncol=N, nrow=1)
    upsilon <- D+3 
    Omega <- diag(D)
    G <- cbind(diag(D-1), -1)
    Xi <- (upsilon-D)*G%*%Omega%*%t(G)
    Theta <- matrix(0, D-1, nrow(X))
    Gamma <- diag(nrow(X))
    
    priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
    priors <- to_clr(priors)  
    names_covariates(priors) <- rownames(X)
    priors$Y <- otu_table 
    posterior <- refit(priors, optim_method="lbfgs")
    rWparaoriginal <- to_proportions(posterior)$Eta
  }
  # Log transform relative abundances
  rWparaoriginal <- log(rWparaoriginal)
  ## END ACCOUNTING FOR UNCERTAINTY IN OBSERVED RELATIVE ABUNDANCES ---------------------------------------------------------------------------------------
  ## ESTIMATING RHO AND SD --------------------------------------------------------------------------------------------------------------------------------
  if (!is.null(externalscalemeasurements) && is.matrix(externalscalemeasurements) && ncol(Y) == nrow(externalscalemeasurements)) {
    rhoandsd_list <- foreach(s = 1:S, .packages = c('stats')) %dopar% {
        n <- length(externalscalemeasurements)
        sample_indices <- boostrap_samples[, s]
        S2 <- var(log(externalscalemeasurements[sample_indices]))
        chi2_lower <- qchisq(alpha / 2, df = n - 1)
        chi2_upper <- qchisq(1 - alpha / 2, df = n - 1)
        var_lower <- (n - 1) * S2 / chi2_upper
        var_upper <- (n - 1) * S2 / chi2_lower
        scalestdev_s <- c(sqrt(var_lower), sqrt(var_upper))
      
        # Initialize matrix for rhobounds for each taxa
        rhobounds_s <- matrix(NA, nrow = D, ncol = 2)  # [D x 2]
        z_critical <- qnorm(1 - alpha / 2)
        for (taxa in 1:D) {
            # Compute correlation
            r <- cor(rWparaoriginal[taxa, sample_indices, s], log(externalscalemeasurements[sample_indices]))
            # Fisher Z-transformation
            z <- 0.5 * log((1 + r) / (1 - r))
            # Standard error of z
            se_z <- 1 / sqrt(n - 3)
            
            # Confidence intervals in Fisher Z-space
            z_lower <- z - z_critical * se_z
            z_upper <- z + z_critical * se_z
            
            # Inverse Fisher Z-transformation to get rho bounds
            rho_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)
            rho_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)
            
            rhobounds_s[taxa, 1] <- rho_lower
            rhobounds_s[taxa, 2] <- rho_upper
        }
        
        list(scalestdev_s = scalestdev_s, rhobounds_s = rhobounds_s)
    }
    # scalestdev_matrix --------------------------------------------------- [S x 2] 
    scalestdev <- do.call(rbind, lapply(rhoandsd_list, function(x) x$scalestdev_s))
    
    # rhobounds ---------------------------------------------------------------------- [D x 2 x S]
    rhobounds <- array(unlist(lapply(rhoandsd_list, function(x) x$rhobounds_s)), dim = c(D, 2, S))
  
    # remove any unneeded space
    rm(rhoandsd_list)

    # plot the scale SD
    plot_scalestdev_histograms(scalestdev, S)
                                     
    # plot the rho
    plot_rho_ridges(rhoubounds, D, S)

  } 
  ## END ESTIMATING RHO AND SD ----------------------------------------------------------------------------------------------------------------------------
  ## ESTIMATING COVARIANCE --------------------------------------------------------------------------------------------------------------------------------
  cat("Running sigma estimation\n")

  # Generate all pairs of indices
  pair_indices <- combn(D, 2, simplify = FALSE)
  
  # Check if pair_indices is populated correctly
  if (length(pair_indices) == 0) {
    stop("Error: pair_indices is not populated correctly. Aborting analysis.")
  }

  # Outer loop is not parallelized
  results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack', 'dplyr'), .options.snow = opts) %do% {
      
      d1 <- pair[1]
      d2 <- pair[2]
      comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
      
      # Use sequential foreach for the inner loop
      results_inner <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack', 'nloptr', 'dplyr')) %dopar% {
        
        rWpara <- rWparaoriginal[c(d1, d2), boostrap_samples[, s], s]
        taxa1relativesd <- sd(rWpara[1, ])
        taxa2relativesd <- sd(rWpara[2, ])
        relativecorrelation <- cor(rWpara[1, ], rWpara[2, ])
        relativecovariance <- cov(rWpara[1, ], rWpara[2, ])
        if (!is.null(rhobounds)) {
          rho_upper_bound_1 = rhobounds[d1, 2, s]
          rho_upper_bound_2 = rhobounds[d2, 2, s]
          rho_lower_bound_1 = rhobounds[d1, 1, s]
          rho_lower_bound_2 = rhobounds[d2, 1, s]
          upperscalestdev = scalestdev[s,2]
          lowerscalestdev = scalestdev[s,1]
        } else {
          rho_lower_bound_1 = lowerrhobound[d1]
          rho_lower_bound_2 = lowerrhobound[d2]
          rho_upper_bound_1 = upperrhobound[d1]
          rho_upper_bound_2 = lowerrhobound[d2]
        }
        initialparameters <- c(((rho_upper_bound_1+rho_lower_bound_1) / 2), ((rho_upper_bound_2+rho_lower_bound_2) / 2), ((upperscalestdev + lowerscalestdev) / 2))
        
        res_min <- NULL
        res_max <- NULL
        
        result <- switch(algorithm,
                         
                         "COBYLA" = {
                           # Find the minimum sigma using nloptr with COBYLA
                           res_min <- nloptr(
                             x0 = initialparameters,
                             eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                             eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                             lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                             ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                             opts = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
                           )
                           
                           # Find the maximum sigma by negating the objective function using COBYLA
                           res_max <- nloptr(
                             x0 = initialparameters,
                             eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                             eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                             lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                             ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                             opts = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
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
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list("algorithm"="NLOPT_LD_MMA", "maxeval" = 10000, "ftol_rel" = 1e-4)
                          )
                          
                          # Find the maximum sigma by negating the objective function using MMA
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list("algorithm"="NLOPT_LD_MMA", "maxeval" = 10000, "ftol_rel" = 1e-4)
                          )
                          
                          list(res_min = res_min, res_max = res_max)
                        },
  
                        "SLSQP" = {
                          # Set optimization options
                          opts <- list(
                            "algorithm" = "NLOPT_LD_SLSQP",
                            "xtol_rel" = 1e-4,
                            "ftol_rel" = 1e-4,
                            "maxeval" = 10000
                          )
                          
                          # Perform the optimization to find the minimum sigma
                          res_min <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = opts
                          )
                          
                          # Perform the optimization to find the maximum sigma (negate the objective function)
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = opts
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
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LN_COBYLA",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
                            )
                          )
                          
                          # Find the maximum sigma using AUGLAG with COBYLA as the inner algorithm
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LN_COBYLA",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
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
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_MMA",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
                            )
                          )
                          
                          # Find the maximum sigma using AUGLAG with MMA as the inner algorithm
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_MMA",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
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
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_SLSQP",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
                            )
                          )
                          
                          # Find the maximum sigma using AUGLAG with SLSQP as the inner algorithm
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_SLSQP",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
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
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_LBFGS",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
                            )
                          )
                          
                          # Find the maximum sigma using AUGLAG with LBFGS as the inner algorithm
                          res_max <- nloptr(
                            x0 = initialparameters,
                            eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_grad_f = function(params) -gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            eval_jac_g_ineq = function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                            lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                            ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                            opts = list(
                              "algorithm" = "NLOPT_LD_AUGLAG",
                              "local_opts" = list(
                                "algorithm" = "NLOPT_LD_LBFGS",
                                "xtol_rel" = 1e-4,
                                "maxeval" = 10000
                              ),
                              "maxeval" = 10000,
                              "ftol_rel" = 1e-4
                            )
                          )
                          
                          list(res_min = res_min, res_max = res_max)
                        },
                         
                         "GRID_SEARCH" = {
                           
                           # Define parameter steps
                           scale <- (upperscalestdev - lowerscalestdev) / 5
                           bound1 <- (rho_upper_bound_1 - rho_lower_bound_1) / 5
                           bound2 <- (rho_upper_bound_2 - rho_lower_bound_2) / 5
                           
                           rho1 <- seq(rho_lower_bound_1, rho_upper_bound_1, by = bound1)
                           rho2 <- seq(rho_lower_bound_2, rho_upper_bound_2, by = bound2)
                           scalestdevstep <- seq(lowerscalestdev, upperscalestdev, by = scale)
                           iterations <- length(rho1) * length(rho2) * length(scalestdevstep)
                           
                           # Create the grid of parameters
                           pars <- expand.grid(rho1 = rho1, rho2 = rho2, scalestdevstep = scalestdevstep, iterations = iterations)
                           
                           # Calculate the constraint
                           constraint_values <- vectorized_constraint_function(
                             pars$rho1,
                             pars$rho2,
                             pars$scalestdevstep,
                             taxa1relativesd,
                             taxa2relativesd,
                             relativecovariance
                           )
                           
                           # Add the calculated constraint as a new column to 'pars'
                           pars$SPSD <- constraint_values
                           
                           # Calculate rpars without filtering based on SPSD
                           rpars <- pars %>%
                             dplyr::mutate(sigma = relativecovariance + scalestdevstep * taxa1relativesd * rho1 + scalestdevstep * taxa2relativesd * rho2 + scalestdevstep^2,
                                           s = s,
                                           comparison = paste0(d1, ":", d2),
                                           taxa1relativesd = taxa1relativesd,
                                           taxa2relativesd = taxa2relativesd,
                                           relativecovariance = relativecovariance
                             )
                           
                           # Append the rpars to the file for this pair
                           pair_file_name <- paste0("gridresults_taxa_", d1, "_", d2, ".txt")
                           append_to_pair_file(rpars, pair_file_name, outputdirectory)
                           
                           # Filter rows where SPSD is >= 0
                           rpars <- rpars %>%
                             dplyr::filter(SPSD >= 0)
    
                           if (nrow(rpars) == 0) {
                             
                             # NULL
                             res_min <- list(
                               objective = NULL,
                               solution = c(NULL, NULL, NULL),
                               message = "GRIDSEARCH",
                               status = "GRIDSEARCH",
                               iterations = NULL
                             )
                             
                             # NULL
                             res_max <- list(
                               objective = NULL,
                               solution = c(NULL, NULL, NULL),
                               message = "GRIDSEARCH",
                               status = "GRIDSEARCH",
                               iterations = NULL
                             )
                           } else {
                             # Find min and max
                             min_sigma_row <- rpars[which.min(rpars$sigma), , drop = FALSE]
                             max_sigma_row <- rpars[which.max(rpars$sigma), , drop = FALSE]
                             
                             # Extract values from min_sigma_row
                             res_min <- list(
                               objective = min_sigma_row$sigma,
                               solution = c(min_sigma_row$rho1, min_sigma_row$rho2, min_sigma_row$scalestdevstep),
                               message = "GRIDSEARCH",
                               status = "GRIDSEARCH",
                               iterations = min_sigma_row$iterations
                             )
                             
                             # Extract values from max_sigma_row
                             res_max <- list(
                               objective = -max_sigma_row$sigma,
                               solution = c(max_sigma_row$rho1, max_sigma_row$rho2, max_sigma_row$scalestdevstep),
                               message = "GRIDSEARCH",
                               status = "GRIDSEARCH",
                               iterations = max_sigma_row$iterations
                             )
                           }
                           
                           list(res_min = res_min, res_max = res_max)
                          
                         },
                         stop("Invalid algorithm selected") # Default case if no match is found
        )
        
        res_min <- result$res_min
        res_max <- result$res_max
        
        data.frame(
          d1 = d1,
          d2 = d2,
          s = s,
          minsigma_absolute_minimum_covariance = ifelse(is.null(res_min$objective), NA, res_min$objective),
          minsigma_correlation_relativetaxa1_scale = ifelse(is.null(res_min$solution[1]), NA, res_min$solution[1]),
          minsigma_correlation_relativetaxa2_scale = ifelse(is.null(res_min$solution[2]), NA, res_min$solution[2]),
          minsigma_scale_sd = ifelse(is.null(res_min$solution[3]), NA, res_min$solution[3]),
          minsigma_message = res_min$message,
          minsigma_status = res_min$status,
          minsigma_iterations = ifelse(is.null(res_min$iterations), NA, res_min$iterations),
          maxsigma_absolute_maximum_covariance = ifelse(is.null(res_max$objective), NA, -res_max$objective),
          maxsigma_correlation_relativetaxa1_scale = ifelse(is.null(res_max$solution[1]), NA, res_max$solution[1]),
          maxsigma_correlation_relativetaxa2_scale = ifelse(is.null(res_max$solution[2]), NA, res_max$solution[2]),
          maxsigma_scale_sd = ifelse(is.null(res_max$solution[3]), NA, res_max$solution[3]),
          maxsigma_message = res_max$message,
          maxsigma_status = res_max$status,
          maxsigma_iterations = ifelse(is.null(res_max$iterations), NA, res_max$iterations),
          taxa1relativesd = taxa1relativesd,
          taxa2relativesd = taxa2relativesd,
          relativecorrelation = relativecorrelation,
          relativecovariance = relativecovariance,
          d1lowerrhobound = rho_lower_bound_1,
          d1upperrhobound = rho_upper_bound_1,
          d2lowerrhobound = rho_lower_bound_2,
          d2upperrhobound = rho_upper_bound_2,
          scalesdlowerbound = lowerscalestdev,
          scalesdupperbound = upperscalestdev
        )
      }

    if (nrow(results_inner) == 0) {
      stop("Error: No valid inner results generated")
    }

    # Check how many bootstraps failed the constraint entirely
    numbootstrapsfailedspsd <- sum(is.na(results_inner$minsigma_absolute_minimum_covariance) | is.na(results_inner$maxsigma_absolute_maximum_covariance))

    # Filter results_inner to keep rows without NA in the specified columns
    results_inner <- results_inner %>%
      filter(!is.na(minsigma_absolute_minimum_covariance) & !is.na(maxsigma_absolute_maximum_covariance))

    # Calculate the number of intervals where both minsigma and maxsigma do not cover zero
    positive_non_zero_intervals <- sum((results_inner$minsigma_absolute_minimum_covariance > 0 & results_inner$maxsigma_absolute_maximum_covariance > 0))
    negative_non_zero_intervals <- sum((results_inner$minsigma_absolute_minimum_covariance < 0 & results_inner$maxsigma_absolute_maximum_covariance < 0))
    proportion_positiveintervals_dontcoverzero <- positive_non_zero_intervals / nrow(results_inner)
    proportion_negativeintervals_dontcoverzero <- negative_non_zero_intervals / nrow(results_inner)
    proportion_intervals_dontcoverzero <- (positive_non_zero_intervals + negative_non_zero_intervals) / nrow(results_inner)
    
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
    min_x <- results_inner$minsigma_scale_sd[min_index]
    max_rho1 <- results_inner$maxsigma_correlation_relativetaxa1_scale[max_index]
    max_rho2 <- results_inner$maxsigma_correlation_relativetaxa2_scale[max_index]
    max_x <- results_inner$maxsigma_scale_sd[max_index]
    
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
        proportion_intervals_dontcoverzero = proportion_intervals_dontcoverzero,
        proportion_positiveintervals_dontcoverzero = proportion_positiveintervals_dontcoverzero,
        proportion_negativeintervals_dontcoverzero = proportion_negativeintervals_dontcoverzero,
        numbootstrapsfailedspsd = numbootstrapsfailedspsd,
        ninetyfive_ci_lower = cilower,
        ninetyfive_ci_upper = ciupper,
        minsigma_absolute_minimum_covariance = minsigma,
        maxsigma_absolute_maximum_covariance = maxsigma,
        minsigma_correlation_relativetaxa1_scale = min_rho1,
        minsigma_correlation_relativetaxa2_scale = min_rho2,
        minsigma_scale_sd = min_x,
        maxsigma_correlation_relativetaxa1_scale = max_rho1,
        maxsigma_correlation_relativetaxa2_scale = max_rho2,
        maxsigma_scale_sd = max_x,
        relative_standard_dev_taxa1 = taxa1relativesd,
        relative_standard_dev_taxa2 = taxa2relativesd,
        relative_correlation = relativecorrelation,
        relative_covariance = relativecovariance,
        d1lowerrhobound = rho_lower_bound_1,
        d1upperrhobound = rho_upper_bound_1,
        d2lowerrhobound = rho_lower_bound_2,
        d2upperrhobound = rho_upper_bound_2,
        scalesdlowerbound = lowerscalestdev,
        scalesdupperbound = upperscalestdev,
        stringsAsFactors = FALSE
      )
    )
  }
  
  ## END SIGMA ESTIMATION --------------------------------------------------------------------------------------------------------------------------------------
  
  # Remove all lock files
  if (is.null(outputdirectory)) {
    outputdirectory <- getwd()
  }
  lock_files <- list.files(path = outputdirectory, pattern = "\\.txt\\.lock$", full.names = TRUE)
  file.remove(lock_files)
  
  # Combine the results into a data frame, transpose it, remove row names
  final_results <- do.call(rbind, lapply(results_list, function(x) x$results))
  final_results <- as.data.frame(final_results, stringsAsFactors = FALSE)
  rownames(final_results) <- NULL

  final_results <- calculate_pval(final_results)
  
  # Combine all inner loop results
  all_inner_results <- do.call(rbind, lapply(results_list, function(x) x$resultsinner))
  all_inner_results <- as.data.frame(all_inner_results)
  all_inner_results$comparison <- paste(rownames(Y)[all_inner_results$d1], rownames(Y)[all_inner_results$d2], sep = ":")
  
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
          minsigma_scale_sd = res_min$par[3],
          maxsigma_correlation_relativetaxa1_scale = res_max$par[1],
          maxsigma_correlation_relativetaxa2_scale = res_max$par[2],
          maxsigma_scale_sd = res_max$par[3],
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
      min_x <- results_inner$minsigma_scale_sd[min_index]
      
      max_rho1 <- results_inner$maxsigma_correlation_relativetaxa1_scale[max_index]
      max_rho2 <- results_inner$maxsigma_correlation_relativetaxa2_scale[max_index]
      max_x <- results_inner$maxsigma_scale_sd[max_index]
      
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
          minsigma_scale_sd = min_x,
          maxsigma_correlation_relativetaxa1_scale = max_rho1,
          maxsigma_correlation_relativetaxa2_scale = max_rho2,
          maxsigma_scale_sd = max_x,
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



      
