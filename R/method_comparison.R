#' Run Covariance Estimation and Comparison Pipeline
#'
#' This function runs a covariance estimation pipeline, comparing several methods (PRISM, SpiecEasi, SparCC, CClasso, Proportionality, and Banocc) on input data.
#' It allows the user to save the resulting plots and outputs with a user-specified filename prefix.
#'
#' @param Y A matrix of abundance data to apply the covariance estimation methods.
#' @param filename_prefix A character string to prefix the filenames of saved plots and results.
#' @param S The number of samples to draw in the PRISM covariance method (default: 4000).
#' @param output_directory The directory to save output files (default: current working directory).
#' @param algorithm The algorithm used for PRISM (default: "GRID_SEARCH").
#' @param uncertaintydistribution The uncertainty distribution used for PRISM (default: "multinomiallognormal").
#' @return A list of covariance results for each method and saved plots.
#' @export
prism.method_comparison <- function(Y, 
                                    flow,
                                    filename_prefix = "simulation_", 
                                    S = 4000, 
                                    output_directory = getwd(),
                                    algorithm = "GRID_SEARCH",
                                    uncertaintydistribution = "multinomiallognormal") {
  cat("Start PRISM\n")
  
  results <- PRISM::prism.covariance(Y = Y, S = S, uncertaintydistribution = uncertaintydistribution,  
                                     externalscalemeasurements = flow, algorithm = algorithm, 
                                     outputdirectory = output_directory)
  
  final_results_filename <- paste0(filename_prefix, "finalresults_", uncertaintydistribution, "_", algorithm, ".csv")
  all_inner_results_filename <- paste0(filename_prefix, "allinnerresults_", uncertaintydistribution, "_", algorithm, ".csv")
  
  write.csv(results$final_results, file = file.path(output_directory, final_results_filename))
  write.csv(results$all_inner_results, file = file.path(output_directory, all_inner_results_filename))
  
  forest_plot_filename <- paste0(filename_prefix, "forestplot_", uncertaintydistribution, "_", algorithm)
  
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
  
  cat("END CClasso\n")
  cat("Start Propr\n")

  prop_result <- tryCatch(propr(t(Y), metric = "rho", ivar = "clr", p = 100),
                              error = function(e) stop("Proportionality failed: ", e$message))
  pr <- updateCutoffs(
          prop_result,
          number_of_cutoffs = 100,  # number of cutoffs to estimate FDR
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
  
  cat("Start Banocc\n")
  max_cores <- parallel::detectCores() - 1  # Reserve 1 core
  ### BANOCC Model ###
  # Pre data
  ps <- phyloseq::phyloseq(otu_table(normalized_data, taxa_are_rows = TRUE))
  ps <- t(ps)
  compiled_banocc_model <- rstan::stan_model(model_code = banocc::banocc_model)
  banocc_results <- banocc::run_banocc(ps, compiled_banocc_model = compiled_banocc_model,
                                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                                           iter = 6000, chains = 4, cores = max_cores, init=list( list(r = 0.1),  # For chain 1
    												list(r = 0.4),  # For chain 2
    												list(r = 0.9),  # For chain 3
    												list(r = 0.7)   # For chain 4
                                                                                        ))
  
  cov_matrix_banocc <- banocc::get_banocc_output(banoccfit = banocc_results, conf_alpha = 0)
  cov_matrix_banocc_CI <- banocc::get_banocc_output(banoccfit = banocc_results, conf_alpha = 0.05)
  
  # Save PRISM results
  saveRDS(prismresults, "PRISM_results.rds")
  
  # Save SpiecEasi results
  saveRDS(spieceasiresults, "SpiecEasi_results.rds")
  
  # Save SparCC results
  saveRDS(sparccresults, "SparCC_results.rds")
  
  # Save CClasso results
  saveRDS(cclassoresults, "CClasso_results.rds")
  
  # Save Proportionality results
  saveRDS(proprresults, "Propr_results.rds")
  
  banoccresults <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Optional, can leave out if not needed
  )
  
  # Loop through the upper tri$CI.hpd$upperangle of the covariance matrix to extract comparisons
  for (i in 1:(ncol(cov_matrix_banocc$Estimates.median) - 1)) {
    for (j in (i + 1):ncol(cov_matrix_banocc$Estimates.median)) {
  
      # Extract the 95% CI bounds for the same pair
      ninetyfive_ci_lower <- cov_matrix_banocc_CI$CI.hpd$lower[i, j]
      ninetyfive_ci_upper <- cov_matrix_banocc_CI$CI.hpd$upper[i, j]
      
      # Extract the min and max bounds for the same pair
      lower <- cov_matrix_banocc$CI.hpd$lower[i, j]
      upper <- cov_matrix_banocc$CI.hpd$upper[i, j]
      
      # Create the comparison name (e.g., "Taxa1:Taxa2")
      comparison_name <- paste0("Taxa", i, ":Taxa", j)
  
      # Append this comparison to the data frame
      banoccresults <- rbind(banoccresults, data.frame(
        comparison = comparison_name,
        ninetyfive_ci_lower = ninetyfive_ci_lower,
        ninetyfive_ci_upper = ninetyfive_ci_upper,
        minsigma_absolute_minimum_covariance = lower,
        maxsigma_absolute_maximum_covariance = upper,
        p_value = NA  # Placeholder, can remove this if not needed
      ))
    }
  }
  
  #Save Banocc results
  saveRDS(banocc_results, "Banocc_results.rds")
  
  cat("END Banocc\n")
  # Visualize and save the plots
  forest_plot_filename <- paste0("simulated_forestplot_all_",uncertaintydistribution,"_",algorithm)
                        
  covarianceresults = list(PRISM = prismresults, 
                           BanoCC = banoccresults, 
                           SpiecEasi = spieceasiresults, 
                           SparCC = sparccresults, 
                           CClasso = cclassoresults, 
                           Propr = proprresults)
  
  if (!is.null(forest_plot_data)) { 
    covarianceresults$TrueAbundances = forest_plot_data
  }
  
  # Save the results as RDS files
  saveRDS(covarianceresults, file = file.path(output_directory, paste0(filename_prefix, "covariance_comparisons.rds")))
  
  # Create forest plot and save
  PRISM::prism.forestplot(covarianceresults, 
                          save = "png", 
                          filename = file.path(output_directory, forest_plot_filename), 
                          color_y_axis_by_ci = TRUE)

  PRISM::prism.circlenetwork(covarianceresults, metric = "covariance", combine_plots = TRUE)
  
  cat("Pipeline complete.\n")
}
