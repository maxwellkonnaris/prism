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
#' @param seed An optional integer value to set the seed for random number generation, ensuring reproducible results. Defaults to \code{NULL}, which means the seed is not set within the function.
#'
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
#' results_custom <- prism.covariance(Y, alpha = alpha_custom, 
#'                                       lowerrhobound = lowerrhobound_custom, 
#'                                       upperrhobound = upperrhobound_custom)
#' 
#' # Example 4: Specifying Output Directory for Grid Search
#' results_grid <- prism.covariance(Y, algorithm = "GRID_SEARCH", 
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
#' @import fido
#' @export
prism.covariance <- function(Y, alpha = 0.5, uncertaintydistribution = "multinomialdirichlet", bootstrap = TRUE, externalscalemeasurements = NULL, 
                             lowerrhobound = rep(-1.0, nrow(Y)), upperrhobound = rep(1.0, nrow(Y)), S = 1000, lowerscalestdev = 0.450, upperscalestdev = 0.650, 
                             algorithm = "GRID_SEARCH", prefix="setprefix", pvalue=TRUE, outputdirectory = NULL, logfile = "log_prismcovariance.txt", seed = NULL) {
      
      ## START LOGGING -----------------------------------------------------------------------------------------------------------------------------------------
      sink(logfile, append = TRUE, type = "message")
      ## END LOGGING --------------------------------------------------------------------------------------------------------------------------------------------
      
      ## COMPUTATIONAL TIME -------------------------------------------------------------------------------------------------------------------------------------
      start_time_total <- Sys.time()  # Start time for the entire function
      message(sprintf("Start time: %s", format_elapsed_time(start_time_total)))
      ## END COMPUTATIONAL TIME SETUP ---------------------------------------------------------------------------------------------------------------------------
      
      ## SETUP --------------------------------------------------------------------------------------------------------------------------------------------------
      # Check if Y is a matrix, dataframe, or tibble, and has appropriate dimensions
      if (!(is.matrix(Y) || is.data.frame(Y) || inherits(Y, "tbl_df")) || nrow(Y) < 2 || ncol(Y) < 2) {
        stop("Y must be a matrix, dataframe, or tibble with at least 2 rows and 2 columns.")
      }
    
      # Set the seed if provided
      if (!is.null(seed)) {
        set.seed(seed)
      }
    
      # Get the number of columns (N) and rows (D) in the input matrix Y
      N <- ncol(Y)
      D <- nrow(Y)
    
      # Check if rownames are NULL, and if so, assign default rownames
      if (is.null(rownames(Y))) {
        message("Taxa labels not specified, default taxa labels have been assigned")
        rownames(Y) <- paste0("Taxa", 1:D)
      }
    
      message("Priors used for the analysis:")
      message(sprintf("Alpha: %s", alpha))
      message("Dimensions of supplied Y matrix:")
      message(sprintf("Number of Taxa: %s", D))
      message(sprintf("Number of Samples: %s", N))
      message(sprintf("Approximating relative counts with the: %s", uncertaintydistribution))
      message(sprintf("Bootstrap sample size (S): %s", S))
      message(sprintf("Algorithm selected: %s", algorithm))
      
      if (!is.null(externalscalemeasurements) && 
          (is.matrix(externalscalemeasurements) || is.vector(externalscalemeasurements))) {
          
          if (is.vector(externalscalemeasurements)) {
              externalscalemeasurements <- matrix(externalscalemeasurements, ncol = 1)
          }
      
          if (ncol(Y) == nrow(externalscalemeasurements)) {
              message("External scale measurements were provided.")
              replicates <- ncol(externalscalemeasurements)
              sampletotals <- nrow(externalscalemeasurements)
              message(sprintf("Dimensions of supplied external scale measurements matrix: Number of Sample-scale Measurement Pairs %s and Number of Replicates %s", sampletotals, replicates))
              message("Estimating Rho bounds and scale SD from the external scale measurements.")
              externalscalemeasurements <- log(externalscalemeasurements)
          } else {
              stop("Error: Mismatch in dimensions between external scale measurements and Y.")
          }
      } else {
          pbounds = paste0(lowerrhobound, ":", upperrhobound)
          sdbounds = paste0(lowerscalestdev, ":", upperscalestdev)
          message(sprintf("Using default Rho bounds: %s ", pbounds))
          message(sprintf("Using default Scale standard deviation bounds %s ", sdbounds))
          rm(pbounds)
          rm(sdbounds)
      }
    
      message(sprintf("Input Y %s", head(Y)))
      ## END SETUP ----------------------------------------------------------------------------------------------------------------------------------------------
      
      ## CLUSTER RESOURCES --------------------------------------------------------------------------------------------------------------------------------------
      num_cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(num_cores)
      doSNOW::registerDoSNOW(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      # Verify cluster registration
      if (!foreach::getDoParRegistered()) {
          stop("Parallel backend is not registered.")
      } else {
          # Check the number of workers
          message(sprintf("Number of workers/cpus: %s ", foreach::getDoParWorkers()))
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
      # Preallocate the matrix of bootstrap samples
      bootstrap_samples <- matrix(NA, nrow = N, ncol = S)
      
      # Populate the bootstrap_samples matrix based on whether bootstrapping is needed
      if (bootstrap) {
        # Generate the bootstrap resampling indices directly
        bootstrap_samples <- replicate(S, sample(1:N, replace = TRUE))
      } else {
        # If no bootstrapping, use the original indices repeated S times
        bootstrap_samples <- matrix(rep(1:N, S), nrow = N, ncol = S, byrow = FALSE)
      }
      ## END ACCOUNTING FOR UNCERTAINTY IN FINITE SAMPLING -----------------------------------------------------------------------------------------------------
      
      ## ACCOUNTING FOR UNCERTAINTY IN OBSERVED RELATIVE ABUNDANCES --------------------------------------------------------------------------------------------
      # Dimensions: (n_taxa, n_samples, n_iter) -- populate matrix of NAs
      rWparaoriginal <- array(NA, dim = c(D, N, S))
      
      # calculate posterior samples -- accounting for uncertainty in the observed relative abundances
      if (uncertaintydistribution == "multinomialdirichlet") {
        # generate S Dirichlet samples for each sample (column)
        for (n in 1:N) {
            # Generate S Dirichlet samples for sample n
            dirichlet_samples <- rdirichlet(S, Y[,n] + alpha)  # Output: S x D
            # Transpose to D x S to match array dimensions
            transposed_samples <- t(dirichlet_samples)         # Output: D x S
            # Assign to the nth slice of the array
            rWparaoriginal[,n,] <- transposed_samples          # Assign D x S
        } 
      } else if (uncertaintydistribution == "multinomiallognormal") {
        # generate S Multinomial logistic Normal posterior samples for each sample (column) using fido
        otu_table = phyloseq::otu_table(Y, taxa_are_rows = TRUE)
        otu_table = otu_table + alpha
        X <- matrix(1, ncol=N, nrow=1)
        upsilon <- D+3 
        Omega <- diag(D)
        G <- cbind(diag(D-1), -1)
        Xi <- (upsilon-D)*G%*%Omega%*%t(G)
        Theta <- matrix(0, D-1, nrow(X))
        Gamma <- diag(nrow(X))
        
        priors <- fido::pibble(otu_table, X, upsilon, Theta, Gamma, Xi, n_samples = S)  
        names_covariates(priors) <- rownames(X)
        posterior <- refit(priors, optim_method="lbfgs")
        rWparaoriginal <- fido::to_proportions(posterior)$Eta
      } else {
        stop("Please select an uncertainty distribution to model the proportions. uncertaintydistribution = {'multinomialdirichlet','multinomiallognormal'}")
      }
        
      # Log transform relative abundances
      rWparaoriginal <- log(rWparaoriginal)
    
      prism.posteriorsamples(rWparaoriginal, filename = prefix)
      ## END ACCOUNTING FOR UNCERTAINTY IN OBSERVED RELATIVE ABUNDANCES ---------------------------------------------------------------------------------------
      
      ## ESTIMATING RHO AND SD --------------------------------------------------------------------------------------------------------------------------------
      message("Start SD and Rho estimation")  
    
      result_rho_sd <- estimate_rho_and_sd(externalscalemeasurements = externalscalemeasurements, Y = Y, S = S, D = D, rWparaoriginal = rWparaoriginal,
                                           bootstrap_samples = bootstrap_samples, alpha = alpha, prefix = prefix)
    
      # Access the results if the function did not return NULL
      if (!is.null(result_rho_sd)) {
        scalestdev <- result_rho_sd$scalestdev
        rhobounds <- result_rho_sd$rhobounds
        message("End SD and Rho estimation") 
      } else {
        print("SD and Rho estimation was incorrectly performed, returned NULL.")
      }           
      ## END ESTIMATING RHO AND SD ----------------------------------------------------------------------------------------------------------------------------
      
      ## ESTIMATING COVARIANCE --------------------------------------------------------------------------------------------------------------------------------
    
      # Share large data objects with the cluster
      parallel::clusterExport(cl, varlist = c("Y", "rWparaoriginal", "bootstrap_samples", "rhobounds", "scalestdev"))
                                         
      message("Running Sigma estimation")
      sigmastart = Sys.time()
      # Generate all pairs of indices
      pair_indices <- combn(D, 2, simplify = FALSE)
      
      # Check if pair_indices is populated correctly
      if (length(pair_indices) == 0) {
        stop("Error: pair_indices is not populated correctly. Aborting analysis.")
      }
    
      # Outer loop is not parallelized
      results_list <- foreach(pair = pair_indices, .packages = c('stats', 'MCMCpack', 'dplyr'), .options.snow = opts) %dopar% {
          
          d1 <- pair[1]
          d2 <- pair[2]
          comparison <- paste(rownames(Y)[d1], rownames(Y)[d2], sep = ":")
          message(sprintf("Start Current Comparison: %s", comparison))
          comparisonstart = Sys.time()
        
          # Use sequential foreach for the inner loop
          results_inner <- foreach(s = 1:S, .combine = 'rbind', .packages = c('stats', 'MCMCpack', 'nloptr', 'dplyr')) %do% {
            
            rWpara <- rWparaoriginal[c(d1, d2), bootstrap_samples[, s], s]
            taxa1relativesd <- sd(rWpara[1, ])
            taxa2relativesd <- sd(rWpara[2, ])
            relativecorrelation <- cor(rWpara[1, ], rWpara[2, ])
            relativecovariance <- cov(rWpara[1, ], rWpara[2, ])
            if (is.null(rhobounds)) {
              rho_lower_bound_1 = lowerrhobound[d1]
              rho_lower_bound_2 = lowerrhobound[d2]
              rho_upper_bound_1 = upperrhobound[d1]
              rho_upper_bound_2 = upperrhobound[d2]
            } else {
              rho_upper_bound_1 = rhobounds[d1, 2, s]
              rho_upper_bound_2 = rhobounds[d2, 2, s]
              rho_lower_bound_1 = rhobounds[d1, 1, s]
              rho_lower_bound_2 = rhobounds[d2, 1, s]
              upperscalestdev = scalestdev[s,2]
              lowerscalestdev = scalestdev[s,1]
            }
            # Define initial parameters for optimization
            initialparameters <- c(((rho_upper_bound_1+rho_lower_bound_1) / 2), ((rho_upper_bound_2+rho_lower_bound_2) / 2), ((upperscalestdev + lowerscalestdev) / 2))
    
            # Function for optimizations (See helper_functions)
            result_sigma = optimize_sigma(d1 = d1, d2 = d2, s = s, taxa1relativesd = taxa1relativesd, taxa2relativesd = taxa2relativesd, relativecovariance = relativecovariance, rho_lower_bound_1 = rho_lower_bound_1, 
                                         rho_lower_bound_2 = rho_lower_bound_2, rho_upper_bound_1 = rho_upper_bound_1, rho_upper_bound_2 = rho_upper_bound_2, lowerscalestdev = lowerscalestdev, upperscalestdev = upperscalestdev, 
                                         initialparameters = initialparameters, algorithm = algorithm, outputdirectory = outputdirectory)
      
            res_min <- result_sigma$res_min
            res_max <- result_sigma$res_max
            
            data.table(
              d1 = d1, 
              d2 = d2, 
              s = s,
              minsigma_absolute_minimum_covariance = ifelse(is.null(res_min$objective), {message(sprintf("res_min$objective is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$objective),
              minsigma_correlation_relativetaxa1_scale = ifelse(is.null(res_min$solution[1]), {message(sprintf("res_min$solution[1] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$solution[1]),
              minsigma_correlation_relativetaxa2_scale = ifelse(is.null(res_min$solution[2]), {message(sprintf("res_min$solution[2] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$solution[2]),
              minsigma_scale_sd = ifelse(is.null(res_min$solution[3]), {message(sprintf("res_min$solution[3] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$solution[3]),
              minsigma_message = ifelse(is.null(res_min$message), {message(sprintf("res_min$message is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$message),
              minsigma_status = ifelse(is.null(res_min$status), {message(sprintf("res_min$status is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$status),
              minsigma_iterations = ifelse(is.null(res_min$iterations), {message(sprintf("res_min$iterations is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_min$iterations),
              maxsigma_absolute_maximum_covariance = ifelse(is.null(res_max$objective), {message(sprintf("res_max$objective is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, -res_max$objective),
              maxsigma_correlation_relativetaxa1_scale = ifelse(is.null(res_max$solution[1]), {message(sprintf("res_max$solution[1] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$solution[1]),
              maxsigma_correlation_relativetaxa2_scale = ifelse(is.null(res_max$solution[2]), {message(sprintf("res_max$solution[2] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$solution[2]),
              maxsigma_scale_sd = ifelse(is.null(res_max$solution[3]), {message(sprintf("res_max$solution[3] is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$solution[3]),
              maxsigma_message = ifelse(is.null(res_max$message), {message(sprintf("res_max$message is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$message),
              maxsigma_status = ifelse(is.null(res_max$status), {message(sprintf("res_max$status is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$status),
              maxsigma_iterations = ifelse(is.null(res_max$iterations), {message(sprintf("res_max$iterations is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, res_max$iterations),
              taxa1relativesd = ifelse(is.null(taxa1relativesd), {message(sprintf("taxa1relativesd is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, taxa1relativesd),
              taxa2relativesd = ifelse(is.null(taxa2relativesd), {message(sprintf("taxa2relativesd is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, taxa2relativesd),
              relativecorrelation = ifelse(is.null(relativecorrelation), {message(sprintf("relativecorrelation is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, relativecorrelation),
              relativecovariance = ifelse(is.null(relativecovariance), {message(sprintf("relativecovariance is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, relativecovariance),
              d1lowerrhobound = ifelse(is.null(rho_lower_bound_1), {message(sprintf("rho_lower_bound_1 is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, rho_lower_bound_1),
              d1upperrhobound = ifelse(is.null(rho_upper_bound_1), {message(sprintf("rho_upper_bound_1 is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, rho_upper_bound_1),
              d2lowerrhobound = ifelse(is.null(rho_lower_bound_2), {message(sprintf("rho_lower_bound_2 is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, rho_lower_bound_2),
              d2upperrhobound = ifelse(is.null(rho_upper_bound_2), {message(sprintf("rho_upper_bound_2 is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, rho_upper_bound_2),
              scalesdlowerbound = ifelse(is.null(lowerscalestdev), {message(sprintf("lowerscalestdev is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, lowerscalestdev),
              scalesdupperbound = ifelse(is.null(upperscalestdev), {message(sprintf("upperscalestdev is NULL for d1 %s and d2 %s in bootstrap %s", d1, d2, s)); NA}, upperscalestdev)
            )
        }
        comparisonend = Sys.time()
        elapsed_time = comparisonend - comparisonstart
        formatted_time <- format_elapsed_time(elapsed_time)
        message(sprintf("Comparison %s --- Total time taken %s", comparison, formatted_time))
        message(sprintf("End Comparison: %s ", comparison))
        
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
        sigmarange <- maxsigma - minsigma 
        ciintervalrange <- ciupper - cilower
        
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
        rho_lower_bound_1 <- results_inner$d1lowerrhobound[min_index]
        rho_upper_bound_1 <- results_inner$d1upperrhobound[min_index]
        rho_lower_bound_2 <- results_inner$d2lowerrhobound[min_index]
        rho_upper_bound_2 <- results_inner$d2upperrhobound[min_index]
        lowerscalestdev <- results_inner$scalesdupperbound[min_index]
        upperscalestdev <- results_inner$scalesdlowerbound[min_index]
        
        # Construct list of dataframes considering the inner results of each optimization and of the final estimates for each pair indices.
        results_inner <- ifelse(is.null(results_inner), {message(sprintf("results_inner is NULL for d1 %s and d2 %s", d1, d2)); ""}, results_inner)
        
        results_df <- data.table(
          comparison = ifelse(is.null(comparison), {message(sprintf("comparison is NULL for d1 %s and d2 %s", d1, d2)); NA}, comparison),
          taxa1 = ifelse(is.null(rownames(Y)[d1]), {message(sprintf("taxa1 is NULL for d1 %s and d2 %s", d1, d2)); NA}, rownames(Y)[d1]),
          taxa2 = ifelse(is.null(rownames(Y)[d2]), {message(sprintf("taxa2 is NULL for d1 %s and d2 %s", d1, d2)); NA}, rownames(Y)[d2]),
          proportion_intervals_dontcoverzero = ifelse(is.null(proportion_intervals_dontcoverzero), {message(sprintf("proportion_intervals_dontcoverzero is NULL for d1 %s and d2 %s", d1, d2)); NA}, proportion_intervals_dontcoverzero),
          proportion_positiveintervals_dontcoverzero = ifelse(is.null(proportion_positiveintervals_dontcoverzero), {message(sprintf("proportion_positiveintervals_dontcoverzero is NULL for d1 %s and d2 %s", d1, d2)); NA}, proportion_positiveintervals_dontcoverzero),
          proportion_negativeintervals_dontcoverzero = ifelse(is.null(proportion_negativeintervals_dontcoverzero), {message(sprintf("proportion_negativeintervals_dontcoverzero is NULL for d1 %s and d2 %s", d1, d2)); NA}, proportion_negativeintervals_dontcoverzero),
          numbootstrapsfailedspsd = ifelse(is.null(numbootstrapsfailedspsd), {message(sprintf("numbootstrapsfailedspsd is NULL for d1 %s and d2 %s", d1, d2)); NA}, numbootstrapsfailedspsd),
          ninetyfive_ci_lower = ifelse(is.null(cilower), {message(sprintf("ninetyfive_ci_lower is NULL for d1 %s and d2 %s", d1, d2)); NA}, cilower),
          ninetyfive_ci_upper = ifelse(is.null(ciupper), {message(sprintf("ninetyfive_ci_upper is NULL for d1 %s and d2 %s", d1, d2)); NA}, ciupper),
          minsigma_absolute_minimum_covariance = ifelse(is.null(minsigma), {message(sprintf("minsigma_absolute_minimum_covariance is NULL for d1 %s and d2 %s", d1, d2)); NA}, minsigma),
          maxsigma_absolute_maximum_covariance = ifelse(is.null(maxsigma), {message(sprintf("maxsigma_absolute_maximum_covariance is NULL for d1 %s and d2 %s", d1, d2)); NA}, maxsigma),
          sigmarange = ifelse(is.null(sigmarange), {message(sprintf("range is NULL for d1 %s and d2 %s", d1, d2)); NA}, sigmarange),
          cirange = ifelse(is.null(ciintervalrange), {message(sprintf("cirange is NULL for d1 %s and d2 %s", d1, d2)); NA}, ciintervalrange),
          minsigma_correlation_relativetaxa1_scale = ifelse(is.null(min_rho1), {message(sprintf("minsigma_correlation_relativetaxa1_scale is NULL for d1 %s and d2 %s", d1, d2)); NA}, min_rho1),
          minsigma_correlation_relativetaxa2_scale = ifelse(is.null(min_rho2), {message(sprintf("minsigma_correlation_relativetaxa2_scale is NULL for d1 %s and d2 %s", d1, d2)); NA}, min_rho2),
          minsigma_scale_sd = ifelse(is.null(min_x), {message(sprintf("minsigma_scale_sd is NULL for d1: %s d2: %s", d1, d2)); NA}, min_x),
          maxsigma_correlation_relativetaxa1_scale = ifelse(is.null(max_rho1), {message(sprintf("maxsigma_correlation_relativetaxa1_scale is NULL for d1 %s and d2 %s", d1, d2)); NA}, max_rho1),
          maxsigma_correlation_relativetaxa2_scale = ifelse(is.null(max_rho2), {message(sprintf("maxsigma_correlation_relativetaxa2_scale is NULL for d1 %s and d2 %s", d1, d2)); NA}, max_rho2),
          maxsigma_scale_sd = ifelse(is.null(max_x), {message(sprintf("maxsigma_scale_sd is NULL for d1: %s d2: %s", d1, d2)); NA}, max_x),
          relative_standard_dev_taxa1 = ifelse(is.null(taxa1relativesd), {message(sprintf("relative_standard_dev_taxa1 is NULL for d1 %s and d2 %s", d1, d2)); NA}, taxa1relativesd),
          relative_standard_dev_taxa2 = ifelse(is.null(taxa2relativesd), {message(sprintf("relative_standard_dev_taxa2 is NULL for d1 %s and d2 %s", d1, d2)); NA}, taxa2relativesd),
          relative_correlation = ifelse(is.null(relativecorrelation), {message(sprintf("relative_correlation is NULL for d1 %s and d2 %s", d1, d2)); NA}, relativecorrelation),
          relative_covariance = ifelse(is.null(relativecovariance), {message(sprintf("relative_covariance is NULL for d1 %s and d2 %s", d1, d2)); NA}, relativecovariance),
          d1lowerrhobound = ifelse(is.null(rho_lower_bound_1), {message(sprintf("d1lowerrhobound is NULL for d1 %s and d2 %s", d1, d2)); NA}, rho_lower_bound_1),
          d1upperrhobound = ifelse(is.null(rho_upper_bound_1), {message(sprintf("d1upperrhobound is NULL for d1 %s and d2 %s", d1, d2)); NA}, rho_upper_bound_1),
          d2lowerrhobound = ifelse(is.null(rho_lower_bound_2), {message(sprintf("d2lowerrhobound is NULL for d1 %s and d2 %s", d1, d2)); NA}, rho_lower_bound_2),
          d2upperrhobound = ifelse(is.null(rho_upper_bound_2), {message(sprintf("d2upperrhobound is NULL for d1 %s and d2 %s", d1, d2)); NA}, rho_upper_bound_2),
          scalesdlowerbound = ifelse(is.null(lowerscalestdev), {message(sprintf("scalesdlowerbound is NULL for d1 %s and d2 %s", d1, d2)); NA}, lowerscalestdev),
          scalesdupperbound = ifelse(is.null(upperscalestdev), {message(sprintf("scalesdupperbound is NULL for d1 %s and d2 %s", d1, d2)); NA}, upperscalestdev),
          stringsAsFactors = FALSE
        )
    
        # Combine the results into a list of two data frames
        results_list <- list(resultsinner = results_inner, results = results_df)
      }
      sigmaend = Sys.time()
      elapsed_time = sigmaend - sigmastart
      message(sprintf("Sigma Total time taken %s", format_elapsed_time(elapsed_time)))
      message("End Sigma estimation")
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
    
      if (pvalue) {
        final_results <- calculate_pval(final_results)
      }
                                             
      # Combine all inner loop results
      all_inner_results <- do.call(rbind, lapply(results_list, function(x) x$resultsinner))
      all_inner_results <- as.data.frame(all_inner_results)
      all_inner_results$comparison <- paste(rownames(Y)[all_inner_results$d1], rownames(Y)[all_inner_results$d2], sep = ":")
    
      # Calculate and print the total elapsed time
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time_total
    
      message(sprintf("Total time taken: %s", format_elapsed_time(elapsed_time)))
      sink(type = "message")
                               
      ## END LOGGING --------------------------------------------------------------------------------------------------------------------------------------------
      return(list(final_results = final_results, all_inner_results = all_inner_results))
    }
