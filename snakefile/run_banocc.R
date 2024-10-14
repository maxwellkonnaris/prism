# scripts/run_banocc.R

library(banocc)
library(phyloseq)
library(rstan)
library(parallel)

# Load data
load('data/abundances.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Prepare data for BAnOCC
normalize_counts <- function(count_matrix) {
  return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
}
normalized_data <- normalize_counts(as.matrix(Y))
ps <- phyloseq::phyloseq(otu_table(normalized_data, taxa_are_rows = TRUE))
ps <- t(ps)

  
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

# Save results
saveRDS(banoccresults, file = "results/banocc_sensitivity_results.rds")
saveRDS(banoccresults_one, file="results/bannocc_one.rds")
saveRDS(banoccresults_two, file="results/bannocc_two.rds")
saveRDS(banoccresults_three, file="results/bannocc_three.rds")
saveRDS(banoccresults_four, file="results/bannocc_four.rds")
