#' Run Covariance Estimation and Comparison Pipeline
#'
#' This function runs a covariance estimation pipeline, comparing several methods (PRISM, SpiecEasi, SparCC, CClasso, Proportionality, and Banocc) on input data.
#' It allows the user to save the resulting plots and outputs with a user-specified filename prefix.
#'
#' @param Y A matrix of abundance data to apply the covariance estimation methods.
#' @param prefix A character string to prefix the filenames of saved plots and results.
#' @param S The number of samples to draw in the PRISM covariance method (default: 4000).
#' @param output_directory The directory to save output files (default: current working directory).
#' @param algorithm The algorithm used for PRISM (default: "GRID_SEARCH").
#' @param uncertaintydistribution The uncertainty distribution used for PRISM (default: "multinomiallognormal").
#' @return A list of covariance results for each method and saved plots.
#' @export
prism.method_comparison <- function(Y, 
                                    externalscalemeasurements = NULL,
                                    trueabundances = NULL,
                                    filename = "simulation_", 
                                    S = 4000, 
                                    output_directory = getwd(),
                                    algorithm = "GRID_SEARCH",
                                    uncertaintydistribution = "multinomiallognormal",
                                    logfile = "log_prismcovariance.txt",
                                    pvalue = FALSE) {

        # Check if externalscalemeasurements is NULL
    if (is.null(externalscalemeasurements)) {
      stop("External scale measurements cannot be NULL")
    }

    num_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    rdat3 <- Y + 0.5
    normalize_counts <- function(count_matrix) {
      return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
    }
    
    normalized_data <- normalize_counts(as.matrix(rdat3))
  
    ps <- phyloseq::phyloseq(otu_table(normalized_data, taxa_are_rows = TRUE))
    ps <- t(ps)

    cat("Start Banocc\n")
    ### BANOCC Model ###
  
    ### Guidelines fo choose gamma parameters for the Glasso ###
    # The shape parameter a determines how quickly the prior density decreasys, while the rate parameter b determines how much prior weight is placed on small λ values rather than large λ values.
    # The shape of the prior on ojk and ojj for several values of λ. Smaller λ results in greater shrinkage towards zero. B-C The prior probability in the interval (−0.001,0.001) for each off-diagonal element ojk| λ∼Laplace(λ) across small (B) or large (C) values of λ. Small values (<0.1) of λ show the greatest shrinkage, while beyond λ = 1 the shrinkage becomes negligible, as shown by the maximal shrinkage for λ > 0.2 being <0.005.


    # Log-Basis Mean
    # We recommend using an uninformative prior for the log-basis mean: centered at zero and with large variance.
    
    # GLASSO Shrinkage Parameter
    # We recommend using a prior with large probability mass close to zero; because λ has a gamma prior, this means that the shape parameter a should be less than one. The rate parameter b determines the variability; in cases with either small (order of 10) or very large (p > n) numbers of features b should be large so that the variance of the gamma distribution, a/b^2, is small. Otherwise, a small value of b will make the prior more uninformative.
        # Define priors for L
    L_val <- list(
      L_prior_1 = 100 * diag(ncol(ps)),
      L_prior_2 = 500 * diag(ncol(ps)),
      L_prior_3 = 100 * diag(ncol(ps)),
      L_prior_4 = 100 * diag(ncol(ps))
    )
    
    # Define priors for alpha and beta (Gamma distribution parameters)
    alpha_val <- list(alpha_prior_1 = 0.5, 
                       alpha_prior_2 = 0.5, 
                       alpha_prior_3 = 1, 
                       alpha_prior_4 = 2)
  
    beta_val <- list(beta_prior_1 = 5, 
                   beta_prior_2 = 5, 
                   beta_prior_3 = 0.5, 
                   beta_prior_4 = 0.5)
    
    # Define fixed initial values
    init_list <- list(
                  list(m = rep(0, ncol(ps)), O = 10*diag(ncol(ps)), lambda = 0.5),
                  list(m = rep(0, ncol(ps)), O = 100*diag(ncol(ps)), lambda = 0.5),
                  list(m = rep(0, ncol(ps)), O = 10*diag(ncol(ps)), lambda = 0.1),
                  list(m = rep(0, ncol(ps)), O = 100*diag(ncol(ps)), lambda = 0.1),
                  list(m = rep(0, ncol(ps)), O = 10*diag(ncol(ps)), lambda = 1),
                  list(m = rep(0, ncol(ps)), O = 100*diag(ncol(ps)), lambda = 1)
                  )
      
  
    # Define the function to run sensitivity analysis with different priors and initial values
    run_banocc_sensitivity <- function(C, compiled_model, n_prior, L_val, alpha_val, beta_val, init_list, chains = 6, iter = 6000, warmup = 3000, cores = num_cores) {
      
      # Store results in a list
      results_list <- list()
      output_list <- list()
      
      # Loop over each set of priors 
      for (i in seq_len(length(L_val))) {
        
        # Get the current L, alpha, beta, and initial values
        L_current <- L_val[[i]]
        a_current <- alpha_val[[i]]
        b_current <- beta_val[[i]]
        
        # Print progress message
        cat("Running model", i, "of", length(L_val), "\n",
            "L =", names(L_val)[i], "\n",
            "alpha =", names(alpha_val)[i], "\n",
            "beta =", names(beta_val)[i], "\n")
        
        # Run the BAnOCC model with the specified priors and initial values
        fit <- banocc::run_banocc(C = C,
                                  compiled_banocc_model = compiled_model,
                                  n = n_prior,
                                  L = L_current,     # L for this model run
                                  a = a_current,     # alpha (shape) for this model run
                                  b = b_current,     # beta (rate) for this model run
                                  init = init_list, # Initial values for this model run
                                  chains = chains,   # Number of chains
                                  iter = iter,       # Total iterations
                                  warmup = warmup,   # Warmup iterations
                                  cores = cores)     # Number of cores for parallel processing
        
        # Store the fit result in the results list
        results_list[[paste0("model_", i)]] <- fit
        
        # Extract the covariance matrix and credible intervals
        cov_matrix <- banocc::get_banocc_output(banoccfit = fit, conf_alpha = 0)
        cov_matrix_CI <- banocc::get_banocc_output(banoccfit = fit, conf_alpha = 0.05)
        
        # Store them in output list
        output_list[[paste0("cov_matrix_", i)]] <- cov_matrix
        output_list[[paste0("cov_matrix_CI_", i)]] <- cov_matrix_CI
      }
      
      # Return the list of results and the extracted covariance matrices
      return(list(fit_results = results_list, cov_matrices = output_list))
    }
    
    # Compile BAnOCC model
    compiled_banocc_model <- rstan::stan_model(model_code = banocc::banocc_model)
    
    # Run the sensitivity analysis function
    banoccresults <- run_banocc_sensitivity(
      C = ps,
      compiled_model = compiled_banocc_model,
      n_prior = rep(0, ncol(ps)),
      L_val = L_val,
      alpha_val = alpha_val,
      beta_val = beta_val,
      param_grid = param_grid,
      init_list = init_list,
      chains = 6,
      iter = 6000,
      warmup = 3000,
      cores = num_cores
    )

  
    # Access the results for each model
    fit_results <- banoccresults$fit_results
    cov_matrices <- banoccresults$cov_matrices
    
    # Extract the covariance matrices and credible intervals
    cov_matrix_model_one <- cov_matrices[[paste0("cov_matrix_", 1)]]
    cov_matrix_CI_model_one <- cov_matrices[[paste0("cov_matrix_CI_", 1)]]
    
    cov_matrix_model_two <- cov_matrices[[paste0("cov_matrix_", 2)]]
    cov_matrix_CI_model_two <- cov_matrices[[paste0("cov_matrix_CI_", 2)]]
    
    cov_matrix_model_three <- cov_matrices[[paste0("cov_matrix_", 3)]]
    cov_matrix_CI_model_three <- cov_matrices[[paste0("cov_matrix_CI_", 3)]]

    cov_matrix_model_four <- cov_matrices[[paste0("cov_matrix_", 4)]]
    cov_matrix_CI_model_four <- cov_matrices[[paste0("cov_matrix_CI_", 4)]]
  
    
    # Function to extract results into a data frame
    extract_cov_results <- function(cov_matrix, cov_matrix_CI, label) {
      results_df <- data.frame(
        comparison = character(),
        ninetyfive_ci_lower = numeric(),
        ninetyfive_ci_upper = numeric(),
        minsigma_absolute_minimum_covariance = numeric(),
        maxsigma_absolute_maximum_covariance = numeric(),
        p_value = numeric()
      )
      
      # Loop through the upper triangular part of the covariance matrix to extract comparisons
      for (i in 1:(ncol(cov_matrix$Estimates.median) - 1)) {
        for (j in (i + 1):ncol(cov_matrix$Estimates.median)) {
          
          # Extract the 95% CI bounds for the same pair
          ninetyfive_ci_lower <- cov_matrix_CI$CI.hpd$lower[i, j]
          ninetyfive_ci_upper <- cov_matrix_CI$CI.hpd$upper[i, j]
          
          # Extract the min and max bounds for the same pair
          lower <- cov_matrix$CI.hpd$lower[i, j]
          upper <- cov_matrix$CI.hpd$upper[i, j]

          taxa_i <- paste0("Taxa", i)
          taxa_j <- paste0("Taxa", j)
          comparison_name <- paste0(taxa_i, ":", taxa_j)
          
          # Append this comparison to the data frame
          results_df <- rbind(results_df, data.frame(
            comparison = comparison_name,
            ninetyfive_ci_lower = ninetyfive_ci_lower,
            ninetyfive_ci_upper = ninetyfive_ci_upper,
            minsigma_absolute_minimum_covariance = lower,
            maxsigma_absolute_maximum_covariance = upper,
            p_value = NA
          ))
        }
      }
      
      # Split the comparison into taxa1 and taxa2
      results_df <- results_df %>%
        separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
      
      return(results_df)
    }
    
    # Extract results for each sensitivity model
    banoccresults_one <- extract_cov_results(cov_matrix_model_one, cov_matrix_CI_model_one, "one")
    banoccresults_two <- extract_cov_results(cov_matrix_model_two, cov_matrix_CI_model_two, "two")
    banoccresults_three <- extract_cov_results(cov_matrix_model_three, cov_matrix_CI_model_three, "three")
    banoccresults_four <- extract_cov_results(cov_matrix_model_four, cov_matrix_CI_model_four, "four")
  
    cat("END Banocc\n")
  
    cat("Start PRISM\n")
    
    # Initialize a list to store results if we have multiple external scale measurements
    result_list <- list()
    
    # If externalscalemeasurements is a single vector (or data frame column)
    if (!is.list(externalscalemeasurements)) {
      cat("Single external scale measurement provided\n")
      
      # Perform PRISM analysis with a single vector of external scale measurements
      results <- PRISM::prism.covariance(Y = Y, S = S, uncertaintydistribution = uncertaintydistribution,  
                                         externalscalemeasurements = externalscalemeasurements, algorithm = algorithm, 
                                         outputdirectory = output_directory, prefix=filename, pvalue=pvalue, logfile = logfile)
      
      # Generate filenames for final and inner results
      final_results_filename <- paste0(filename, "finalresults_", uncertaintydistribution, "_", algorithm, ".csv")
      all_inner_results_filename <- paste0(filename, "allinnerresults_", uncertaintydistribution, "_", algorithm, ".csv")
      
      # Save results to CSV files
      write.csv(results$final_results, file = file.path(output_directory, final_results_filename))
      write.csv(results$all_inner_results, file = file.path(output_directory, all_inner_results_filename))
      
      # Optionally generate and save forest plot filename (if needed)
      forest_plot_filename <- paste0(filename, "forestplot_", uncertaintydistribution, "_", algorithm)
      
      # Store the final results
      prismresults <- results$final_results
      
    } else {
      # If externalscalemeasurements is a list of vectors/data frame columns
      cat("Multiple external scale measurements provided\n")
      
      # Iterate over the list of vectors and run PRISM on each
      for (i in seq_along(externalscalemeasurements)) {
        cat("Processing external scale measurement", i, "\n")
        
        # Extract the current external scale measurement vector (it could be a data frame column or vector)
        current_measurement <- as.matrix(externalscalemeasurements[[i]])
        
        # Run PRISM for each vector of external scale measurements
        results <- PRISM::prism.covariance(Y = Y, S = S, uncertaintydistribution = uncertaintydistribution,  
                                           externalscalemeasurements = current_measurement, algorithm = algorithm, 
                                           outputdirectory = output_directory, prefix=filename, pvalue=pvalue, logfile = logfile)
        
        # Append iteration number to filenames to avoid overwriting
        final_results_filename <- paste0(filename, "finalresults_", uncertaintydistribution, "_", algorithm, "_", i, ".csv")
        all_inner_results_filename <- paste0(filename, "allinnerresults_", uncertaintydistribution, "_", algorithm, "_", i, ".csv")
        
        # Save the results to CSV files
        write.csv(results$final_results, file = file.path(output_directory, final_results_filename))
        write.csv(results$all_inner_results, file = file.path(output_directory, all_inner_results_filename))
        
        # Optionally generate and save forest plot filename for each iteration (if needed)
        forest_plot_filename <- paste0(filename, "forestplot_", uncertaintydistribution, "_", algorithm, "_", i)
        
        # Store the results for this iteration in the result list
        result_list[[i]] <- results$final_results
      }
      
      # Assign the final results list to prismresults
      prismresults <- result_list
  }
    

  if (is.list(externalscalemeasurements)) {
    for (i in seq_along(prismresults)) {
      # Dynamically create a variable name like prismresults_1, prismresults_2, etc.
      variable_name <- paste0("prismresults_", i)
      
      # Assign each element of prismresults to the dynamically created variable name
      assign(variable_name, prismresults[[i]], envir = .GlobalEnv)
    }
  }
  
  cat("END PRISM\n")

  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  cat("Start SpiecEasi\n")

  lambda_min_ratios <- c(1e-1, 1e-2, 1e-3, 1e-4) ## parameter defines the minimum ratio of the regularization path. smaller lambda min ratio will explore more sparse models (i.e., fewer edges in the network). larger lambda min ratio will allow denser networks. Smaller Networks or Less Complexity: If you expect the microbial network to be sparse (e.g., few interactions), a lower value like 1e-3 or 1e-4 might be useful. Larger Networks or Complex Interactions: If you believe there is more interaction complexity, keeping the ratio at 1e-2 or even increasing it slightly might help find more connections without overfitting.
  nlambda_values <- 100 ## We are not converned about computational efficiency, but if you are then this is a parameter you would decrease to explore less lambda iterations. 100 is chosen for increased precision.
  stars_thresholds <- c(0.01, 0.05, 0.1) ## Low Threshold (e.g., 0.01): This means you are more strict about accepting network variability, resulting in a sparser and more stable network. High Threshold (e.g., 0.2): This is more lenient, allowing for more variability and resulting in a denser network with potentially more false positives. Default Value: A value of 0.05 is often used, which balances stability with network density. If You Want Higher Stability: A lower value like 0.01 or 0.02 ensures that the connections in the network are stable across subsamples. This is useful if you want to be conservative and prioritize stable, reliable edges. If You Are Tolerant to Unstable Edges: You might increase the threshold to 0.1 or 0.2 if you are willing to accept a denser network and more potential false positives.
  rep_num_values <- 100 ## Again we are not concerned with computational efficiency here, so 100 is chosen for better stability from subsampling stars. specifies the minimum fraction of subsamples where an edge must appear to be considered stable.

  results_mb <- list()
  results_gl <- list()

 # Define the grid of parameters to test
  lambda_min_ratios <- c(1e-1, 1e-2, 1e-3, 1e-4) # parameter defines the minimum ratio of the regularization path. smaller lambda min ratio will explore more sparse models (i.e., fewer edges in the network). larger lambda min ratio will allow denser networks. Smaller Networks or Less Complexity: If you expect the microbial network to be sparse (e.g., few interactions), a lower value like 1e-3 or 1e-4 might be useful. Larger Networks or Complex Interactions: If you believe there is more interaction complexity, keeping the ratio at 1e-2 or even increasing it slightly might help find more connections without overfitting.
  nlambda_values <- 100 # We are not converned about computational efficiency, but if you are then this is a parameter you would decrease to explore less lambda iterations. 100 is chosen for increased precision.
  stars_thresholds <- c(0.01, 0.05, 0.1) #Low Threshold (e.g., 0.01): This means you are more strict about accepting network variability, resulting in a sparser and more stable network. High Threshold (e.g., 0.2): This is more lenient, allowing for more variability and resulting in a denser network with potentially more false positives. Default Value: A value of 0.05 is often used, which balances stability with network density. If You Want Higher Stability: A lower value like 0.01 or 0.02 ensures that the connections in the network are stable across subsamples. This is useful if you want to be conservative and prioritize stable, reliable edges. If You Are Tolerant to Unstable Edges: You might increase the threshold to 0.1 or 0.2 if you are willing to accept a denser network and more potential false positives.
  rep_num_values <- 100 # Again we are not concerned with computational efficiency here, so 100 is chosen for better stability from subsampling stars. specifies the minimum fraction of subsamples where an edge must appear to be considered stable.
  
  # Create a results list to store outputs
  results_mb <- list()
  results_gl <- list()
  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Loop over each combination of hyperparameters
  counter <- 1
  for (lambda_min in lambda_min_ratios) {
      for (thresh in stars_thresholds) {

          # Run spiec.easi with the current combination of parameters
          tryCatch({
            #regression coefficients (in neighborhood selection)
            mbspiec_easi_result <- spiec.easi(t(as.matrix(Y)), 
                                            method = "mb", 
                                            lambda.min.ratio = lambda_min, 
                                            nlambda = nlambda_values,  
                                            sel.criterion = "stars",  
                                            pulsar.params = list(rep.num = rep_num_values, thresh = thresh, seed=10241994, ncores=num_cores),
                                            pulsar.select = TRUE,
                                            verbose = TRUE)
            # precision matrix (in graphical lasso)                           
            glspiec_easi_result <- spiec.easi(t(as.matrix(Y)), 
                                            method = "glasso", 
                                            lambda.min.ratio = lambda_min, 
                                            nlambda = nlambda_values,  
                                            sel.criterion = "stars",  
                                            pulsar.params = list(rep.num = rep_num_values, thresh = thresh, seed=10241994, ncores=num_cores),
                                            pulsar.select = TRUE)

            # Store the spiec_easi_result and additional metrics in the results list
            results_mb[[counter]] <- list(method = "mb",
                                      lambda_min = lambda_min,
                                      nlambda = nlambda_values,
                                      thresh = thresh,
                                      rep_num = rep_num_values,
                                      num_edges = sum(getRefit(mbspiec_easi_result)),
                                      sparsity = sum(getRefit(mbspiec_easi_result)) / (ncol(t(Y)) * (ncol(t(Y)) - 1) / 2),
                                      beta_matrix_spiec_easi = getOptBeta(mbspiec_easi_result),
                                      symmetricbeta_matrix_spiec_easi = symBeta(getOptBeta(mbspiec_easi_result), mode='maxabs'),
                                      stability = getStability(mbspiec_easi_result),
                                      num_edges_pair = sum(getRefit(mbspiec_easi_result))/2,
                                      spiec_easi_result = mbspiec_easi_result)  # Store full spiec_easi result
            
            # Store the spiec_easi_result and additional metrics in the results list
            results_gl[[counter]] <- list(method = "glasso",
                                      lambda_min = lambda_min,
                                      nlambda = nlambda_values,
                                      thresh = thresh,
                                      rep_num = rep_num_values,
                                      num_edges = sum(getRefit(glspiec_easi_result)),
                                      sparsity = sum(getRefit(glspiec_easi_result)) / (ncol(t(Y)) * (ncol(t(Y)) - 1) / 2),
                                      precision_matrix_spiec_easi = getOptiCov(glspiec_easi_result),
                                      cov_matrix_spiec_easi = getOptCov(glspiec_easi_result),
                                      stability = getStability(glspiec_easi_result),
                                      num_edges_pair = sum(getRefit(glspiec_easi_result))/2,
                                      spiec_easi_result = glspiec_easi_result)  # Store full spiec_easi result
            
            counter <- counter + 1
            
          }, error = function(e) {
            message(paste("Error with lambda_min:", lambda_min, "nlambda:", nlambda_values, "thresh:", thresh, "rep_num:", rep_num_values))
            message(e)
          })
        }
      }

  # Define a function to select the best model
  select_best_model <- function(results_list) {
    best_model <- NULL
    best_stability <- 0
    for (res in results_list) {
      # You may want to prioritize models with highest stability and reasonable sparsity
      if (res$stability > best_stability && res$num_edges > 0) {
        best_stability <- res$stability
        best_model <- res
      }
    }
    return(best_model)
  }
  
  # Select the best model for both methods (MB and GL)
  best_mb_model <- select_best_model(results_mb)
  best_gl_model <- select_best_model(results_gl)

  # Extract the covariance matrix from the best GL model
  cov_matrix_spiec_easi <- best_gl_model$cov_matrix_spiec_easi
  
  spieceasiresults <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()
  )
  
  # Loop through the upper triangle of the covariance matrix to extract comparisons
  for (i in 1:(ncol(cov_matrix_spiec_easi) - 1)) {
    for (j in (i + 1):ncol(cov_matrix_spiec_easi)) {
      cov_value <- cov_matrix_spiec_easi[i, j]
      comparison_name <- paste0("Taxa", i, ":Taxa", j)
      spieceasiresults <- rbind(spieceasiresults, data.frame(
        comparison = comparison_name,
        ninetyfive_ci_lower = cov_value,
        ninetyfive_ci_upper = cov_value,
        minsigma_absolute_minimum_covariance = cov_value,
        maxsigma_absolute_maximum_covariance = cov_value,
        p_value = NA
      ))
    }
  }

  spieceasiresults <- spieceasiresults %>%
          separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
  
  cat("END SpiecEasi\n")
  cat("Start SparCC\n")
  
  sparcc_result <- tryCatch(PRISM::SparCC_count(x = t(Y)), 
                            error = function(e) stop("SparCC failed: ", e$message))
  
  cov_matrix_sparcc <- sparcc_result$cov_w
  cov_matrix_sparcc_95lower <- sparcc_result$cov_ci_lower
  cov_matrix_sparcc_95upper <- sparcc_result$cov_ci_upper
  
  sparccresults <- data.frame(
      comparison = character(),
      ninetyfive_ci_lower = numeric(),
      ninetyfive_ci_upper = numeric(),
      minsigma_absolute_minimum_covariance = numeric(),
      maxsigma_absolute_maximum_covariance = numeric(),
      p_value = numeric()
    )
    
  for (i in 1:(ncol(cov_matrix_sparcc) - 1)) {
      for (j in (i + 1):ncol(cov_matrix_sparcc)) {
        cov_value <- cov_matrix_sparcc[i, j]
        comparison_name <- paste0("Taxa", i, ":Taxa", j)
        sparccresults <- rbind(sparccresults, data.frame(
          comparison = comparison_name,
          ninetyfive_ci_lower = cov_matrix_sparcc_95lower[i,j],
          ninetyfive_ci_upper = cov_matrix_sparcc_95upper[i,j],
          minsigma_absolute_minimum_covariance = cov_value,
          maxsigma_absolute_maximum_covariance = cov_value,
          p_value = NA
        ))
      }
  }
  
  sparccresults <- sparccresults %>%
            separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
  
  cat("END SparCC\n")
  cat("Start CClasso\n")
  
  cclasso_result <- tryCatch(PRISM::cclasso(t(Y), counts = TRUE, pseudo = 0.5, k_cv = 10, lam_int = c(1e-4, 10), k_max = 100, n_boot = 100), 
                             error = function(e) stop("CCLasso failed: ", e$message))
  
  variances <- cclasso_result$var_w
  correlations <- cclasso_result$cor_w
  cov_matrix_cclasso <- correlations * sqrt(outer(variances, variances))
  pvals = cclasso_result$p_vals
  
  cclassoresults <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()
  )
  
  for (i in 1:(ncol(cov_matrix_cclasso) - 1)) {
    for (j in (i + 1):ncol(cov_matrix_cclasso)) {
      cov_value <- cov_matrix_cclasso[i, j]
      comparison_name <- paste0("Taxa", i, ":Taxa", j)
      cclassoresults <- rbind(cclassoresults, data.frame(
        comparison = comparison_name,
        ninetyfive_ci_lower = cov_value,
        ninetyfive_ci_upper = cov_value,
        minsigma_absolute_minimum_covariance = cov_value,
        maxsigma_absolute_maximum_covariance = cov_value,
        p_value = pvals[i,j]
      ))
    }
  }
                             
  cclassoresults <- cclassoresults %>%
          separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
                             
  cat("END CClasso\n")
  cat("Start Propr\n")

  # Function to calculate proportionality metrics and store results in a list
  calculate_proportionality <- function(normalized_data) {
    
    metrics <- c("rho", "phi", "phs")  # Define the metrics to be used
    tails_options <- c("both", "right", "right")  # Different tails options for rho, phi, and phs
    
    results_list <- list()  # Empty list to store results
    
    for (k in 1:length(metrics)) {
      
      metric <- metrics[k]
      tails <- tails_options[k]
      
      # Try-Catch block to handle proportionality calculation errors
      prop_obj <- tryCatch(propr(t(normalized_data), metric = metric, ivar = "clr", alpha=0, p = 100),
                           error = function(e) stop("Proportionality failed: ", e$message))
      
      # Apply cutoff criteria
      pr_obj <- updateCutoffs(
      prop_obj,
      tails = tails,  # 'both' for rho, 'right' for phi and phs
      ncores = num_cores  # Parallel processing
      )
  
      # Assuming pr_obj is the propr object
      corr_matrix_proportionality <- getMatrix(pr_obj)
  
      # Get the data from getResults() for the corresponding metric (rho, phs, etc.)
      results_data <- getSignificantResultsFDR(pr_obj)
  
      # Create a comparison column in results_data to facilitate merging
      results_data$comparison <- paste0(results_data$Pair, ":", results_data$Partner)
  
      # Initialize an empty data frame to store the results
      propr_results <- data.frame(
        comparison = character(),
        ninetyfive_ci_lower = numeric(),
        ninetyfive_ci_upper = numeric(),
        minsigma_absolute_minimum_covariance = numeric(),
        maxsigma_absolute_maximum_covariance = numeric(),
        p_value = numeric()  # Optional, placeholder for p-values if needed
      )
  
      # Loop through the upper triangle of the covariance matrix to extract comparisons
      for (i in 1:(ncol(corr_matrix_proportionality) - 1)) {
        for (j in (i + 1):ncol(corr_matrix_proportionality)) {
          
          # Create the comparison name (e.g., "Taxa1:Taxa2")
          comparison_name <- paste0("Taxa", i, ":", "Taxa", j)
          
          # Check if this comparison exists in the results_data table
          if (comparison_name %in% results_data$comparison) {
            # If present, extract the propr value from results_data
            propr_value <- results_data$propr[results_data$comparison == comparison_name]
          } else {
            # If not present, set propr_value as NA
            propr_value <- NaN
          }
          
          # Append the result to the data frame
          propr_results <- rbind(propr_results, data.frame(
            comparison = comparison_name,
            ninetyfive_ci_lower = propr_value,  # Fill this with propr_value
            ninetyfive_ci_upper = propr_value,  # Fill this with propr_value
            minsigma_absolute_minimum_covariance = propr_value,
            maxsigma_absolute_maximum_covariance = propr_value,
            p_value = NA  # Placeholder, can be removed if p-values are not needed
          ))
        }
      }
  
      # Separate the comparison column into taxa1 and taxa2 columns
  
      propr_results <- propr_results %>%
        separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
  
      
      # Store the results in the list
      results_list[[metric]] <- list(
        dataframe = propr_results,
        pr_object = pr_obj
      )
    }
    
    # Return the list containing results for rho, phi, and phs
    return(results_list)
  }
  
  # Example usage with your normalized data
  proportionality_results <- calculate_proportionality(normalized_data)
  
  # Accessing the results for rho, phi, and phs
  proprresults_rho <- proportionality_results[["rho"]][["dataframe"]]
  proprresults_phi <- proportionality_results[["phi"]][["dataframe"]]
  proprresults_phs <- proportionality_results[["phs"]][["dataframe"]]       

  cat("END Propr\n")
  
  #Save Banocc results
  saveRDS(banoccresults, file = "banocc_sensitivity_results.rds")
  saveRDS(banoccresults_weak, paste0(filename,"Banocc_results_weak.rds"))
  saveRDS(banoccresults_md, paste0(filename,"Banocc_results_md.rds"))
  saveRDS(banoccresults_strong, paste0(filename,"Banocc_results_strong.rds"))
                          
  # Save PRISM results
  save(prismresults, paste0(filename,"PRISM_results.Rdata"))
  
  # Save SpiecEasi results
  saveRDS(spieceasiresults, paste0(filename,"SpiecEasi_results.rds"))
  save(results_mb, results_gl, file = "spiec_easi_results.RData")
                           
  # Save SparCC results
  saveRDS(sparccresults, paste0(filename,"SparCC_results.rds"))
  
  # Save CClasso results
  saveRDS(cclassoresults, paste0(filename,"CClasso_results.rds"))

  # Save the entire proportionality_results list to an RData file
  save(proportionality_results, file = "proportionality_results.RData")

  # Visualize and save the plots
  forest_plot_filename <- paste0("simulated_forestplot_all_",uncertaintydistribution,"_",algorithm)
                        
  # Initialize covarianceresults with the fixed components first
  covarianceresults <- list(
    BanoCC_one = banoccresults_one, 
    BanoCC_two = banoccresults_two, 
    BanoCC_three = banoccresults_three, 
    BanoCC_four = banoccresults_four, 
    SpiecEasi = spieceasiresults, 
    SparCC = sparccresults, 
    CClasso = cclassoresults, 
    Propr_rho = proprresults_rho,
    Propr_phi = proprresults_phi,
    Propr_phs = proprresults_phs
  )
  
  # Add the PRISM results dynamically based on the length of prismresults
  if (is.list(prismresults)) {
    for (i in seq_along(prismresults)) {
      # Create a name for each PRISM result, first is 'PRISM', the rest 'PRISM_2', 'PRISM_3', etc.
      if (i == 1) {
        result_name <- "PRISM"
      } else {
        result_name <- paste0("PRISM_", i)
      }
      
      # Add each PRISM result to the covarianceresults list
      covarianceresults[[result_name]] <- prismresults[[i]]
    }
  }

  # Extract only the PRISM results from covarianceresults
  prism_results <- covarianceresults[grepl("^PRISM", names(covarianceresults))]
  
  # Extract the remaining (non-PRISM) results from covarianceresults
  other_results <- covarianceresults[!grepl("^PRISM", names(covarianceresults))]
  
  # Combine them, putting PRISM results first
  covarianceresults <- c(prism_results, other_results)

  
  if (!is.null(trueabundances)) { 
    covarianceresults$TrueAbundances = trueabundances
  }
  
  # Save the results as RDS files
  saveRDS(covarianceresults, file = file.path(paste0(filename, "covariance_comparisons.rds")))
  
  # Create forest plot and save
  PRISM::prism.forestplot(covarianceresults, 
                          save = "png", 
                          filename = file.path(output_directory, filename, forest_plot_filename), 
                          color_y_axis_by_ci = TRUE)

  PRISM::prism.circlenetwork(covarianceresults, metric = "covariance", combine_plots = TRUE, filename=filename)
  
  cat("Pipeline complete.\n")
}
