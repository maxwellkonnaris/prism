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
  cat("Start PRISM\n")
  
  results <- PRISM::prism.covariance(Y = Y, S = S, uncertaintydistribution = uncertaintydistribution,  
                                     externalscalemeasurements = externalscalemeasurements, algorithm = algorithm, 
                                     outputdirectory = output_directory, prefix=filename, pvalue=pvalue, logfile = logfile)
  
  final_results_filename <- paste0(filename, "finalresults_", uncertaintydistribution, "_", algorithm, ".csv")
  all_inner_results_filename <- paste0(filename, "allinnerresults_", uncertaintydistribution, "_", algorithm, ".csv")
  
  write.csv(results$final_results, file = file.path(output_directory, final_results_filename))
  write.csv(results$all_inner_results, file = file.path(output_directory, all_inner_results_filename))
  
  forest_plot_filename <- paste0(filename, "forestplot_", uncertaintydistribution, "_", algorithm)
  prismresults = results$final_results
  cat("END PRISM\n")
  
  rdat3 <- Y + 0.5
  normalize_counts <- function(count_matrix) {
    return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
  }
  
  normalized_data <- normalize_counts(as.matrix(rdat3))
  
  cat("Start SpiecEasi\n")
  spiec_easi_result <- tryCatch(spiec.easi(t(as.matrix(Y)), 
                                           method = "glasso", 
                                           lambda.min.ratio = 1e-4, 
                                           nlambda = 100,  
                                           sel.criterion = "stars",  
                                           pulsar.params = list(rep.num = 100, thresh = 0.01)),
                               error = function(e) stop("SPIEC-EASI failed: ", e$message))
  
  cov_matrix_spiec_easi <- getOptCov(spiec_easi_result)
  
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
  
  cov_matrix_sparcc <- sparcc_result$cor_w
  
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
        ninetyfive_ci_lower = cov_value,
        ninetyfive_ci_upper = cov_value,
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
  
  cclasso_result <- tryCatch(PRISM::cclasso(t(normalized_data)), 
                             error = function(e) stop("CCLasso failed: ", e$message))
  
  variances <- cclasso_result$var_w
  correlations <- cclasso_result$cor_w
  cov_matrix_cclasso <- correlations * sqrt(outer(variances, variances))
  
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
        p_value = NA
      ))
    }
  }
                             
  cclassoresults <- cclassoresults %>%
          separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
                             
  cat("END CClasso\n")
  cat("Start Propr\n")

  prop_result <- tryCatch(propr(t(Y), metric = "phi", ivar = "clr", p = 100),
                              error = function(e) stop("Proportionality failed: ", e$message))
  pr <- updateCutoffs(
          prop_result,
          custom_cutoffs = c(0.1),  # number of cutoffs to estimate FDR
          tails = 'right',  # consider only the positive values ('right') or both sides ('both')
          ncores = 4  # parallelize here
        ) 
  corr_matrix_proportionality <- pr@matrix
  
  
  proprresults <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Optional, can leave out if not needed
  )
  
  # Loop through the upper triangle of the covariance matrix to extract comparisons
  for (i in 1:(ncol(corr_matrix_proportionality) - 1)) {
    for (j in (i + 1):ncol(corr_matrix_proportionality)) {
      # Extract the covariance value
      cov_value <- corr_matrix_proportionality[i, j]
  
      # Create the comparison name (e.g., "Taxa1:Taxa2")
      comparison_name <- paste0(colnames(corr_matrix_proportionality)[i], ":", colnames(corr_matrix_proportionality)[j])
      comparison_name <- paste0("Taxa",i, ":", "Taxa",j) # fix and remove
      
      # Append this comparison to the data frame
      proprresults <- rbind(proprresults, data.frame(
        comparison = comparison_name,
        ninetyfive_ci_lower = cov_value,
        ninetyfive_ci_upper = cov_value,
        minsigma_absolute_minimum_covariance = cov_value,
        maxsigma_absolute_maximum_covariance = cov_value,
        p_value = NA  # Placeholder, can remove this if not needed
      ))
    }
  }

  proprresults <- proprresults %>%
          separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)                        
  
  cat("Start Banocc\n")
  max_cores <- parallel::detectCores() - 1  # Reserve 1 core
  ### BANOCC Model ###

  # Define the function to run sensitivity analysis with different priors and initial values
  run_banocc_sensitivity <- function(C, compiled_model, n_prior, L_list, a_list, b_list, init_list, chains = 4, iter = 6000, warmup = 3000, cores = 4) {
    
    # Ensure the lists of priors are the same length
    if (length(L_list) != length(a_list) || length(L_list) != length(b_list) || length(L_list) != length(init_list)) {
      stop("L_list, a_list, b_list, and init_list must have the same length!")
    }
    
    # Store results in a list
    results_list <- list()
    output_list <- list()
    
    # Loop over each set of priors and initial values
    for (i in seq_along(L_list)) {
      
      # Print progress message
      cat("Running model", i, "with L =", L_list[[i]], "a =", a_list[[i]], "b =", b_list[[i]], "\n")
      
      # Run the BAnOCC model with the specified priors and initial values
      fit <- banocc::run_banocc(C = C,
                                compiled_banocc_model = compiled_model,
                                n = n_prior,
                                L = L_list[[i]],   # L for this model run
                                a = a_list[[i]],   # a (shape) for this model run
                                b = b_list[[i]],   # b (rate) for this model run
                                init = init_list[[i]], # Initial values for this model run
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
  
  # Define priors and initial values for sensitivity analysis
  ### Guidelines fo choose gamma parameters for the Glasso ###
  # Prior covariance matrix for m ----- Scaling by 10 is a weakly informative choice, allowing moderate spread for the values of m. A higher value (e.g., 100) would imply less confidence in the prior, giving m more freedom to take extreme values.
  # Conversely, a smaller value (e.g., 1) would shrink m closer to the prior mean.
  #More confident priors: Use a smaller scaling factor (e.g., 𝐿=1×diag(𝑝) 
  #Less confident priors: Increase the scaling factor (e.g., 𝐿=100×diag(𝑝) to allow for more flexibility in the model estimates.
  
  ### Guidelines for Choosing α\alphaα and β\betaβ:
  #1. **Weakly Informative Priors**:
  #   - If you do not have strong prior information about λ\lambdaλ, you can use weakly informative priors. A common choice is to set both α\alphaα and β\betaβ to values that encourage shrinkage but allow the model to explore other values.- **Suggested values**: α=0.5\alpha = 0.5α=0.5, β=0.01\beta = 0.01β=0.01
  #        - This choice places significant probability mass near zero, encouraging shrinkage while not completely ruling out larger values of λ\lambdaλ.2. **Based on Prior Knowledge**:
  #   - If you have domain knowledge or previous studies that suggest what values of λ\lambdaλ are reasonable, you can adjust α\alphaα and β\betaβ to reflect that.
  #        - For example, if you know λ\lambdaλ is likely to be small but nonzero, you might use α=1\alpha = 1α=1 and β=0.1\beta = 0.1β=0.1, leading to a prior that concentrates around small positive values.3. **Controlling Shrinkage**:
  #    - If you want more aggressive shrinkage (forcing λ\lambdaλ closer to zero), set a small value for α\alphaα (e.g., α=0.1\alpha = 0.1α=0.1) and a larger β\betaβ (e.g., β=1\beta = 1β=1).- If you want to reduce the amount of shrinkage and give the model more freedom, you can increase α\alphaα (e.g., α=2\alpha = 2α=2) and set β\betaβ to a smaller value (e.g., β=0.5\beta = 0.5β=0.5).

  # Priors for L (covariance matrix)
  L_prior_1 <- 1 * diag(ncol(ps))
  L_prior_2 <- 10 * diag(ncol(ps))
  L_prior_3 <- 100 * diag(ncol(ps))
  
  # Priors for a (shape) and b (rate)
  a_prior_1 <- 0.5; b_prior_1 <- 0.01
  a_prior_2 <- 1; b_prior_2 <- 0.1
  a_prior_3 <- 0.1; b_prior_3 <- 1
  
  # Initial values for each model
  init_1 <- list(
    list(m = rep(0, ncol(ps)), O = diag(ncol(ps)), lambda = 0.02),
    list(m = runif(ncol(ps), -0.1, 0.1), O = diag(ncol(ps)), lambda = runif(1, 0.01, 0.1)),
    list(m = rnorm(ncol(ps), 0, 0.5), O = diag(ncol(ps)), lambda = 0.05),
    list(m = rep(-1, ncol(ps)), O = diag(ncol(ps)), lambda = 0.01)
  )
  
  init_2 <- list(
    list(m = rep(0, ncol(ps)), O = diag(ncol(ps)), lambda = 0.5),
    list(m = runif(ncol(ps), -1, 1), O = diag(ncol(ps)), lambda = runif(1, 0.1, 2)),
    list(m = rnorm(ncol(ps)), O = diag(ncol(ps)), lambda = runif(1, 0.05, 1)),
    list(m = rep(-1, ncol(ps)), O = diag(ncol(ps)), lambda = 0.01)
  )
  
  init_3 <- list(
    list(m = rep(0, ncol(ps)), O = diag(ncol(ps)), lambda = 1),
    list(m = runif(ncol(ps), -5, 5), O = diag(ncol(ps)), lambda = runif(1, 0.5, 3)),
    list(m = rnorm(ncol(ps), 0, 2), O = diag(ncol(ps)), lambda = runif(1, 0.5, 2)),
    list(m = rep(-2, ncol(ps)), O = diag(ncol(ps)), lambda = 0.1)
  )
  
  # Combine priors into lists for sensitivity analysis
  L_list <- list(L_prior_1, L_prior_2, L_prior_3)
  a_list <- list(a_prior_1, a_prior_2, a_prior_3)
  b_list <- list(b_prior_1, b_prior_2, b_prior_3)
  init_list <- list(init_1, init_2, init_3)
  
  # Run the sensitivity analysis function
  banoccresults <- run_banocc_sensitivity(C = ps,
                                    compiled_model = compiled_banocc_model,
                                    n_prior = rep(0, ncol(ps)),
                                    L_list = L_list,
                                    a_list = a_list,
                                    b_list = b_list,
                                    init_list = init_list,
                                    chains = 4,
                                    iter = 6000,
                                    warmup = 3000,
                                    cores = max_cores)
  
  # Access the results for each model
  fit_results <- banoccresults$fit_results
  cov_matrices <- banoccresults$cov_matrices
  
  cov_matrix_CI_model_lo <- cov_matrices$cov_matrix_CI_1
  cov_matrix_CI_model_md <- cov_matrices$cov_matrix_CI_2
  cov_matrix_CI_model_lg <- cov_matrices$cov_matrix_CI_3

    # Initialize empty data frames for each sensitivity run
  banoccresults_lo <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Optional, can leave out if not needed
  )
  
  banoccresults_md <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Optional, can leave out if not needed
  )
  
  banoccresults_lg <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Optional, can leave out if not needed
  )
  
  # Function to extract results into a data frame
  extract_cov_results <- function(cov_matrix, cov_matrix_CI, label) {
    results_df <- data.frame(
      comparison = character(),
      ninetyfive_ci_lower = numeric(),
      ninetyfive_ci_upper = numeric(),
      minsigma_absolute_minimum_covariance = numeric(),
      maxsigma_absolute_maximum_covariance = numeric(),
      p_value = numeric()  # Optional, can leave out if not needed
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
        
        # Create the comparison name (e.g., "Taxa1:Taxa2")
        comparison_name <- paste0("Taxa", i, ":Taxa", j)
        
        # Append this comparison to the data frame
        results_df <- rbind(results_df, data.frame(
          comparison = comparison_name,
          ninetyfive_ci_lower = ninetyfive_ci_lower,
          ninetyfive_ci_upper = ninetyfive_ci_upper,
          minsigma_absolute_minimum_covariance = lower,
          maxsigma_absolute_maximum_covariance = upper,
          p_value = NA  # Placeholder, can remove or compute if needed
        ))
      }
    }

    results_df <- results_df %>%
          separate(comparison, into = c("taxa1", "taxa2"), sep = ":", remove = FALSE)
    
    return(results_df)
  }
  
  # Extract results for each sensitivity model
  banoccresults_lo <- extract_cov_results(cov_matrix_banocc_lo, cov_matrix_banocc_CI_lo, "lo")
  banoccresults_md <- extract_cov_results(cov_matrix_banocc_md, cov_matrix_banocc_CI_md, "md")
  banoccresults_lg <- extract_cov_results(cov_matrix_banocc_lg, cov_matrix_banocc_CI_lg, "lg")

  cat("END Banocc\n")
  
  #Save Banocc results
  saveRDS(banoccresults, file = "banocc_sensitivity_results.rds")
  saveRDS(banoccresults_lo, paste0(filename,"Banocc_results_lo.rds"))
  saveRDS(banoccresults_md, paste0(filename,"Banocc_results_md.rds"))
  saveRDS(banoccresults_lg, paste0(filename,"Banocc_results_lg.rds"))
                          
  # Save PRISM results
  saveRDS(prismresults, paste0(filename,"PRISM_results.rds"))
  
  # Save SpiecEasi results
  saveRDS(spieceasiresults, paste0(filename,"SpiecEasi_results.rds"))
  
  # Save SparCC results
  saveRDS(sparccresults, paste0(filename,"SparCC_results.rds"))
  
  # Save CClasso results
  saveRDS(cclassoresults, paste0(filename,"CClasso_results.rds"))
  
  # Save Proportionality results
  saveRDS(proprresults, paste0(filename,"Propr_results.rds"))

  # Visualize and save the plots
  forest_plot_filename <- paste0("simulated_forestplot_all_",uncertaintydistribution,"_",algorithm)
                        
  covarianceresults = list(PRISM = prismresults, 
                           BanoCC_lo = banoccresults_lo, 
                           BanoCC_md = banoccresults_md, 
                           BanoCC_lg = banoccresults_lg, 
                           SpiecEasi = spieceasiresults, 
                           SparCC = sparccresults, 
                           CClasso = cclassoresults, 
                           Propr = proprresults)
  
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
