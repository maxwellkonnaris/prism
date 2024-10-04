#' Multivariate Normal Random Variable Generator
#'
#' This function generates random samples from a multivariate normal distribution.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param mu A numeric vector specifying the mean of the distribution.
#' @param Sigma A numeric matrix specifying the covariance matrix of the distribution.
#' @return A numeric matrix of random samples, where each column represents a sample.
#' @examples
#' # Example usage:
#' mu <- c(0, 0)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' samples <- rmvnorm(100, mu, Sigma)
rmvnorm <- function(n, mu, Sigma) {
  p <- length(mu)
  r <- matrix(rnorm(n * p), p, n)
  L <- t(chol(Sigma))
  r <- L %*% r
  sweep(r, 1, mu, FUN = `+`)
}

#' Check if Value is Within Range
#'
#' This function checks if a given value is within a specified range.
#'
#' @param x A numeric value to check.
#' @param l A numeric value specifying the lower bound of the range.
#' @param u A numeric value specifying the upper bound of the range.
#' @return A logical value indicating whether \code{x} is within the range [\code{l}, \code{u}].
#' @examples
#' # Example usage:
#' within(5, 1, 10) # Returns TRUE
#' within(0, 1, 10) # Returns FALSE
within <- function(x, l, u) {
  x >= l && x <= u
}

#' Calculate and Adjust P-Values from Confidence Intervals
#'
#' This function calculates p-values for the estimated confidence intervals from the \code{estimate_covariance} function and adjusts them for multiple hypothesis testing using the Bonferroni correction and the Benjamini-Hochberg FDR method.
#'
#' @param results A dataframe containing the results from the \code{estimate_covariance} function, including the confidence intervals.
#' @return A dataframe with the original results plus columns for the calculated p-values, Bonferroni-adjusted p-values, and FDR-adjusted p-values.
#' @import stats
#' @examples
#' # Example usage:
#' # Assuming `results` is the output dataframe from the `estimate_covariance` function
#' adjusted_results <- calculate_pval(results)
calculate_pval <- function(results) {
  # Check if the necessary columns are present in the results dataframe
  required_columns <- c("ninetyfive_ci_lower", "ninetyfive_ci_upper", "comparison")
  if (!all(required_columns %in% colnames(results))) {
    stop("The input results dataframe must contain the following columns: ninetyfive_ci_lower, ninetyfive_ci_upper, comparison")
  }
  
  # Function to calculate p-value from confidence interval
  get_p_value <- function(lower, upper) {
    if (lower > 0 || upper < 0) {
      return(2 * min(pnorm(lower, lower.tail = FALSE), pnorm(upper)))
    } else {
      return(2 * min(pnorm(0, lower, upper - lower), pnorm(0, upper, upper - lower)))
    }
  }
  
  # Initialize p_value, bonferroni_p_value, and bh_p_value columns with NA
  results$p_value <- NA
  results$bonferroni_p_value <- NA
  results$bh_p_value <- NA
  
  # Filter rows without NA in confidence interval columns
  valid_rows <- !is.na(results$ninetyfive_ci_lower) & !is.na(results$ninetyfive_ci_upper)
  
  if (sum(valid_rows) > 0) {
    # Calculate p-values for valid rows
    results$p_value[valid_rows] <- apply(results[valid_rows, ], 1, function(row) {
      lower <- as.numeric(row["ninetyfive_ci_lower"])
      upper <- as.numeric(row["ninetyfive_ci_upper"])
      get_p_value(lower, upper)
    })
    
    # Number of hypotheses/tests
    m <- sum(valid_rows)
    
    # Bonferroni correction
    results$bonferroni_p_value[valid_rows] <- pmin(results$p_value[valid_rows] * m, 1)
    
    # Benjamini-Hochberg FDR method
    results_valid <- results[valid_rows, ]
    results_valid <- results_valid[order(results_valid$p_value), ]
    results_valid$bh_p_value <- p.adjust(results_valid$p_value, method = "BH")
    
    # Assign the BH p-values back to the original results
    results$bh_p_value[valid_rows] <- results_valid$bh_p_value
  }
  
  return(results)
}


#' Calculate Skewness
#'
#' This function calculates the skewness of a numeric vector.
#' 
#' Skewness is a measure of the asymmetry of the probability distribution of a real-valued random variable about its mean. 
#' The skewness value can be positive or negative, or undefined. 
#' A negative skewness indicates that the tail on the left side of the probability density function is longer or fatter than the right side.
#' A positive skewness indicates that the tail on the right side is longer or fatter than the left side.
#' 
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be removed before the computation. Defaults to FALSE.
#' @return The skewness of the numeric vector.
#' @examples
#' data <- c(1, 2, 2, 3, 4, 6, 7, 8, 9, NA)
#' skewness(data, na.rm = TRUE)
skewness <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  skew <- (n / ((n - 1) * (n - 2))) * sum(((x - mean_x) / sd_x)^3)
  return(skew)
}

#' Calculate Kurtosis
#'
#' This function calculates the kurtosis of a numeric vector.
#'
#' Kurtosis is a measure of the "tailedness" of the probability distribution of a real-valued random variable. 
#' The kurtosis of the normal distribution is 3. Excess kurtosis is often defined as kurtosis minus 3, 
#' such that the standard normal distribution has a kurtosis of zero.
#' Positive kurtosis indicates a distribution with heavier tails and a sharper peak than the normal distribution.
#' Negative kurtosis indicates a distribution with lighter tails and a flatter peak than the normal distribution.
#' 
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be removed before the computation. Defaults to FALSE.
#' @return The kurtosis of the numeric vector.
#' @examples
#' data <- c(1, 2, 2, 3, 4, 6, 7, 8, 9, NA)
#' kurtosis(data, na.rm = TRUE)
kurtosis <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  kurt <- (n * (n + 1) / ((n - 1) * (n - 2) * (n - 3))) * sum(((x - mean_x) / sd_x)^4) -
              (3 * (n - 1)^2 / ((n - 2) * (n - 3)))
  return(kurt)
}

#' Calculate Comprehensive Summary Statistics
#'
#' This function calculates a comprehensive set of summary statistics (mean, standard deviation, variance, median, minimum, 1st quartile, 3rd quartile, and maximum) for each numeric column in a dataframe, grouped by a specified column. Optionally, the summary statistics can be saved as a CSV file.
#' This is meant to provide summary statistics for the bootstrap results of the \code{estimate_covariance()} function in PRISM.
#'
#' Bootstrap estimation is a powerful statistical technique, but its accuracy can be affected by several factors:
#' - **Sample Size**: Larger samples tend to provide more accurate bootstrap estimates.
#' - **Number of Resamples**: More resamples generally lead to more stable and accurate estimates, but also increase computational cost.
#' - **Distribution of Data**: The underlying distribution of the data can affect the bootstrap results. Non-normal data can lead to biased estimates.
#' - **Presence of Outliers**: Outliers can significantly impact the estimates of summary statistics.
#' - **Dependence Structure**: Dependencies between data points can violate the assumptions of the bootstrap method, affecting the accuracy of the estimates.
#'
#' @param results A list containing the dataframes. The function expects a dataframe named \code{`all_inner_results`} within the results list from estimate_covariance().
#' @param group_col A string specifying the name of the column to group by. Default is \code{"comparison"}.
#' @param exclude_cols A character vector of column names to exclude from the summary statistics calculation. Default is \code{c(d1, d2, s)}.
#' @param save_as_csv A logical value indicating whether to save the summary statistics as a CSV file. Default is FALSE.
#' @param csv_path A string specifying the file path to save the CSV file if `save_as_csv` is TRUE. Default is "bootstrap_summary.csv".
#' 
#' @return A dataframe with the calculated summary statistics for each numeric column, grouped by the specified column. If `save_as_csv` is TRUE, the summary statistics are also saved as a CSV file.
#' @export
#'
#' @examples
#' # Load necessary libraries
#' library(dplyr)
#'
#' # Example results list with all_inner_results dataframe
#' results <- list(
#'   all_inner_results = data.frame(
#'     comparison = rep(c("Group1", "Group2"), each = 5),
#'     value1 = rnorm(10),
#'     value2 = runif(10),
#'     d1 = rnorm(10),
#'     d2 = rnorm(10)
#'   )
#' )
#'
#' # Calculate summary statistics excluding columns 'd1' and 'd2'
#' summary_stats <- prism.summary_stats(results$all_inner_results, group_col = "comparison", exclude_cols = c("d1", "d2"), save_as_csv = TRUE, csv_path = "bootstrap_summary_stats.csv")
#' print(summary_stats)
prism.bootstrap_summary <- function(data, group_col = "comparison", exclude_cols = c("d1", "d2", "s"), save_as_csv = FALSE, csv_path = "bootstrap_summary.csv") {
  
  # Ensure the group column and exclude columns are character vectors
  group_col <- as.character(group_col)
  exclude_cols <- as.character(exclude_cols)
  
  # Select columns to include in the summary statistics calculation
  data <- data %>%
    select(-all_of(exclude_cols))
  
  # Calculate summary statistics grouped by the specified column
  summary_df <- data %>%
    group_by_at(group_col) %>%
    summarise_if(is.numeric, list(
      mean = ~mean(., na.rm = TRUE),
      sd = ~sd(., na.rm = TRUE),
      variance = ~var(., na.rm = TRUE),
      median = ~median(., na.rm = TRUE),
      min = ~min(., na.rm = TRUE),
      q1 = ~quantile(., probs = 0.25, na.rm = TRUE),
      q3 = ~quantile(., probs = 0.75, na.rm = TRUE),
      max = ~max(., na.rm = TRUE),
      kurtosis = ~PRISM::kurtosis(., na.rm = TRUE),
      skewness = ~PRISM::skewness(., na.rm = TRUE)
    ), .groups = 'drop')
  
  # Save the summary statistics as a CSV file if required
  if (save_as_csv) {
    write.csv(summary_df, csv_path, row.names = FALSE)
  }

  return(summary_stats)
}

#' Calculate Proportionality Metrics (Rho) for All Pairs of Taxa
#'
#' This function calculates the proportionality metric (\eqn{\rho}) for all pairs of taxa
#' in a given dataset of relative abundances. The \eqn{\rho} metric quantifies the consistency
#' of the log-ratios of relative abundances between pairs of taxa across multiple samples.
#'
#' @param data A data frame of relative abundances where rows are samples and columns are taxa.
#' @return A named vector of \eqn{\rho} values for each pair of taxa.
#' @examples
#' # Example relative abundance data
#' relative_abundances <- data.frame(
#'   Taxon1 = c(0.5, 0.6, 0.55, 0.7),
#'   Taxon2 = c(0.3, 0.25, 0.35, 0.2),
#'   Taxon3 = c(0.2, 0.15, 0.1, 0.1),
#'   Taxon4 = c(0.1, 0.1, 0.05, 0.05)
#' )
#' # Calculate proportionality metrics
#' rho_values <- prism.proportionality(relative_abundances)
#' print(rho_values)
#' @export
prism.proportionality <- function(data) {
  # Check for zeros in the data
  if (any(data == 0)) {
    stop("Data contains zero values. Please remove or replace zeros before calculating proportionality metrics.")
  }
  
  # Define a function to calculate the rho proportionality metric
  proportionality <- function(x, y) {
    log_ratio <- log(x / y)
    var_log_ratio <- var(log_ratio, na.rm = TRUE)

    return(var_log_ratio)
  }
  
  # Calculate rho for each pair of taxa
  taxa_pairs <- combn(names(data), 2, simplify = FALSE)
  rho_values <- sapply(taxa_pairs, function(pair) {
    proportionality(data[[pair[1]]], data[[pair[2]]])
  })
  
  # Convert rho values to a named vector
  names(rho_values) <- sapply(taxa_pairs, function(pair) paste(pair, collapse = "-"))
  
  return(rho_values)
}
#' Format elapsed time in a user-friendly format
#'
#' This function takes an elapsed time and formats it in a user-friendly string 
#' with days, hours, minutes, and seconds.
#'
#' @param elapsed_time A difftime object representing the elapsed time.
#' @return A string representing the formatted elapsed time.
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

#' Append results to a file with error handling
#'
#' This function appends results to a specified file, creating a lock file to 
#' avoid concurrency issues.
#'
#' @param results A data frame or matrix to be written to the file.
#' @param pair_file_name The name of the file where the results will be appended.
#' @param outputdirectory The directory where the file will be stored. If NULL, 
#' the current working directory is used.
#' @return NULL
append_to_pair_file <- function(results, pair_file_name, outputdirectory) {
  tryCatch({
    if (is.null(outputdirectory)) {
      outputdirectory <- getwd()
    }
    if (substr(outputdirectory, nchar(outputdirectory), nchar(outputdirectory)) != "/") {
      outputdirectory <- paste0(outputdirectory, "/")
    }
    pair_file_name <- paste0(outputdirectory, pair_file_name)
    lock_file <- paste0(pair_file_name, ".lock")
    
    # Acquire file lock before writing
    lock <- filelock::lock(lock_file)
    
    # Check if the file exists to decide whether to add the header or not
    if (file.exists(pair_file_name)) {
      # Append data without header if the file already exists
      write.table(results, file = pair_file_name, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    } else {
      # Write data with header if the file does not exist
      write.table(results, file = pair_file_name, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    
  }, error = function(e) {
    message("Error while writing to file: ", pair_file_name, "\n", e)
  }, finally = {
    # Ensure the lock is released in the finally block
    filelock::unlock(lock)
  })
}

#' Objective function for optimizing covariance
#'
#' This function computes the covariance based on relative standard deviations 
#' and correlation parameters.
#'
#' @param params A numeric vector with correlation and standard deviation parameters.
#' @param taxa1relativesd Relative standard deviation for Taxa 1.
#' @param taxa2relativesd Relative standard deviation for Taxa 2.
#' @param relativecovariance The relative covariance.
#' @return A numeric value representing the covariance.
objective_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
  taxa1scalecorrelation <- params[1]
  taxa2scalecorrelation <- params[2]
  scalestdev <- params[3]
  
  sigma <- relativecovariance + scalestdev * taxa1relativesd * taxa1scalecorrelation + scalestdev * taxa2relativesd * taxa2scalecorrelation + scalestdev^2
  
  return(sigma)
}

#' Gradient function for covariance optimization
#'
#' This function calculates the gradient of the covariance function with respect 
#' to the correlation and standard deviation parameters.
#'
#' @param params A numeric vector with correlation and standard deviation parameters.
#' @param taxa1relativesd Relative standard deviation for Taxa 1.
#' @param taxa2relativesd Relative standard deviation for Taxa 2.
#' @param relativecovariance The relative covariance.
#' @return A numeric vector representing the gradient.
gradient_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
  taxa1scalecorrelation <- params[1]
  taxa2scalecorrelation <- params[2]
  scalestdev <- params[3]
  
  grad_rho1 <- taxa1relativesd * scalestdev
  grad_rho2 <- taxa2relativesd * scalestdev
  grad_x <- taxa1relativesd * taxa1scalecorrelation + taxa2relativesd * taxa2scalecorrelation + 2 * scalestdev
  
  return(c(grad_rho1, grad_rho2, grad_x))
}

#' Constraint function for Symmetric Positive Semi-Definite case
#'
#' This function checks whether the covariance matrix is Symmetric Positive 
#' Semi-Definite (SPSD) based on given parameters.
#'
#' @param params A numeric vector with correlation and standard deviation parameters.
#' @param taxa1relativesd Relative standard deviation for Taxa 1.
#' @param taxa2relativesd Relative standard deviation for Taxa 2.
#' @param relativecovariance The relative covariance.
#' @return A numeric value representing the constraint evaluation.
constraint_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
  
  taxa1scalecorrelation <- params[1]
  taxa2scalecorrelation <- params[2]
  scalestdev <- params[3]
  
  term1 <- taxa1relativesd * taxa1scalecorrelation + taxa2relativesd * taxa2scalecorrelation
  term2 <- sqrt(2) * sqrt((taxa1relativesd^2 * taxa1scalecorrelation^2) + (taxa2relativesd^2 * taxa2scalecorrelation^2))
  
  term3 <- (1 / (2 * scalestdev)) * ((taxa1relativesd^2 + taxa2relativesd^2) -
                                       sqrt((taxa1relativesd^2 - taxa2relativesd^2)^2 + 4 * relativecovariance^2) + 4 * scalestdev^2)
  
  g_1 <- term1 - term2 + term3
  
  return(g_1)
}

#' Vectorized constraint function for covariance optimization
#'
#' This function computes the constraint for the SPSD condition in a vectorized 
#' manner, enabling efficient computation over multiple parameter sets.
#'
#' @param rho1 A numeric vector for correlation values for Taxa 1.
#' @param rho2 A numeric vector for correlation values for Taxa 2.
#' @param scalestdevstep A numeric vector for standard deviation steps.
#' @param taxa1relativesd Relative standard deviation for Taxa 1.
#' @param taxa2relativesd Relative standard deviation for Taxa 2.
#' @param relativecovariance The relative covariance.
#' @return A numeric vector representing the constraint evaluation for each set of parameters.
vectorized_constraint_function <- function(rho1, rho2, scalestdevstep, taxa1relativesd, taxa2relativesd, relativecovariance) {
  n <- length(rho1)
  if (length(taxa1relativesd) == 1) taxa1relativesd <- rep(taxa1relativesd, n)
  if (length(taxa2relativesd) == 1) taxa2relativesd <- rep(taxa2relativesd, n)
  if (length(relativecovariance) == 1) relativecovariance <- rep(relativecovariance, n)
  
  term1 <- taxa1relativesd * rho1 + taxa2relativesd * rho2
  term2 <- sqrt(2) * sqrt(taxa1relativesd^2 * rho1^2 + taxa2relativesd^2 * rho2^2)
  term3 <- (1 / (2 * scalestdevstep)) * ((taxa1relativesd^2 + taxa2relativesd^2) -
                                           sqrt((taxa1relativesd^2 - taxa2relativesd^2)^2 + 4 * relativecovariance^2) + 4 * scalestdevstep^2)
  
  g_1 <- term1 - term2 + term3
  
  return(g_1)
}

#' Gradient of the constraint function for covariance optimization
#'
#' This function calculates the gradient of the constraint function for SPSD 
#' covariance matrices with respect to the correlation and standard deviation parameters.
#'
#' @param params A numeric vector with correlation and standard deviation parameters.
#' @param taxa1relativesd Relative standard deviation for Taxa 1.
#' @param taxa2relativesd Relative standard deviation for Taxa 2.
#' @param relativecovariance The relative covariance.
#' @return A numeric vector representing the gradient of the constraint function.
constraint_gradient_function <- function(params, taxa1relativesd, taxa2relativesd, relativecovariance) {
  taxa1scalecorrelation <- params[1]
  taxa2scalecorrelation <- params[2]
  scalestdev <- params[3]
  
  epsilon <- 1e-8
  denominator <- sqrt(taxa1relativesd^2 * taxa1scalecorrelation^2 + taxa2relativesd^2 * taxa2scalecorrelation^2 + epsilon)
  
  grad_taxa1scalecorrelation <- taxa1relativesd - sqrt(2) * (taxa1relativesd^2 * taxa1scalecorrelation) / denominator
  grad_taxa2scalecorrelation <- taxa2relativesd - sqrt(2) * (taxa2relativesd^2 * taxa2scalecorrelation) / denominator
  
  term_to_simplify <- (taxa1relativesd^2 + taxa2relativesd^2) - sqrt((taxa1relativesd^2 - taxa2relativesd^2)^2 + 4 * relativecovariance^2)
  grad_scalestdev <- -1 / (2 * (scalestdev + epsilon)^2) * term_to_simplify + 2 * scalestdev
  
  return(c(grad_taxa1scalecorrelation, grad_taxa2scalecorrelation, grad_scalestdev))
}

#' Custom Cumulative Variance Function
#'
#' This function calculates the cumulative variance for a vector of data.
#'
#' @param x A numeric vector for which to calculate cumulative variance.
#'
#' @return A numeric vector containing the cumulative variance of `x`.
#' @keywords internal
cumvar <- function(x) {
  n <- length(x)
  cumsum((x - cumsum(x) / seq_along(x))^2) / seq_along(x)
}

#' Calculate Effective Sample Size (ESS)
#'
#' This function computes the effective sample size for each bootstrap estimate.
#'
#' @param bootstrap_estimates A data frame of bootstrap estimates containing columns `minsigma_absolute_minimum_covariance` and `maxsigma_absolute_maximum_covariance`.
#'
#' @return A data frame with the effective sample size for both the minimum and maximum covariance estimates.
#' @keywords internal
calculate_effective_sample_size <- function(bootstrap_estimates) {
  ess <- function(x) {
    n <- length(x)
    acf_x <- acf(x, plot = FALSE)
    return(n / (1 + 2 * sum(acf_x$acf[-1])))
  }
  bootstrap_estimates %>%
    group_by(S, comparison) %>%
    summarise(
      ess_min = ess(minsigma_absolute_minimum_covariance),
      ess_max = ess(maxsigma_absolute_maximum_covariance)
    )
}

#' Calculate Rhat (Gelman-Rubin Diagnostic)
#'
#' This function calculates the Rhat diagnostic to assess convergence using the Gelman-Rubin diagnostic.
#'
#' @param bootstrap_estimates A data frame of bootstrap estimates containing columns `minsigma_absolute_minimum_covariance` and `maxsigma_absolute_maximum_covariance`.
#'
#' @return A list containing the Rhat values for both minimum and maximum covariance estimates.
#' @keywords internal
calculate_rhat <- function(bootstrap_estimates) {
  chains_min <- split(bootstrap_estimates$minsigma_absolute_minimum_covariance, bootstrap_estimates$S)
  chains_max <- split(bootstrap_estimates$maxsigma_absolute_maximum_covariance, bootstrap_estimates$S)
  mcmc_chains_min <- coda::mcmc.list(lapply(chains_min, mcmc))
  mcmc_chains_max <- coda::mcmc.list(lapply(chains_max, mcmc))
  list(
    rhat_min = coda::gelman.diag(mcmc_chains_min)$psrf,
    rhat_max = coda::gelman.diag(mcmc_chains_max)$psrf
  )
}

#' Jackknife-after-Bootstrap
#'
#' This function calculates the jackknife-after-bootstrap diagnostic to estimate bias and variability.
#'
#' @param bootstrap_estimates A data frame of bootstrap estimates containing columns `minsigma_absolute_minimum_covariance` and `maxsigma_absolute_maximum_covariance`.
#'
#' @return A list containing jackknife estimates for both minimum and maximum covariance estimates.
#' @keywords internal
jackknife_after_bootstrap <- function(bootstrap_estimates) {
  n <- length(bootstrap_estimates$minsigma_absolute_minimum_covariance)
  list(
    jackknife_min = sapply(1:n, function(i) {
      mean(bootstrap_estimates$minsigma_absolute_minimum_covariance[-i])
    }),
    jackknife_max = sapply(1:n, function(i) {
      mean(bootstrap_estimates$maxsigma_absolute_maximum_covariance[-i])
    })
  )
}

#' Bootstrap-after-Bootstrap
#'
#' This function performs bootstrap resampling on the bootstrap estimates to calculate variability.
#'
#' @param bootstrap_estimates A data frame of bootstrap estimates containing columns `minsigma_absolute_minimum_covariance` and `maxsigma_absolute_maximum_covariance`.
#' @param num_resamples The number of bootstrap resamples to perform (default is 1000).
#'
#' @return A list containing bootstrap resampling results for both minimum and maximum covariance estimates.
#' @keywords internal
bootstrap_after_bootstrap <- function(bootstrap_estimates, num_resamples = 1000) {
  list(
    bootstrap_min = replicate(num_resamples, {
      resample_indices <- sample(seq_along(bootstrap_estimates$minsigma_absolute_minimum_covariance), replace = TRUE)
      mean(bootstrap_estimates$minsigma_absolute_minimum_covariance[resample_indices])
    }),
    bootstrap_max = replicate(num_resamples, {
      resample_indices <- sample(seq_along(bootstrap_estimates$maxsigma_absolute_maximum_covariance), replace = TRUE)
      mean(bootstrap_estimates$maxsigma_absolute_maximum_covariance[resample_indices])
    })
  )
}

#' Calculate Monte Carlo Standard Error (MCSE)
#'
#' This function calculates the Monte Carlo standard error for each bootstrap estimate.
#'
#' @param bootstrap_estimates A data frame of bootstrap estimates containing columns `minsigma_absolute_minimum_covariance` and `maxsigma_absolute_maximum_covariance`.
#'
#' @return A data frame with the Monte Carlo standard error for both minimum and maximum covariance estimates.
#' @keywords internal
calculate_mcse <- function(bootstrap_estimates) {
  mcse <- function(x) {
    sd(x) / sqrt(length(x))
  }
  bootstrap_estimates %>%
    group_by(S, comparison) %>%
    summarise(
      mcse_min = mcse(minsigma_absolute_minimum_covariance),
      mcse_max = mcse(maxsigma_absolute_maximum_covariance)
    )
}

#' Simulate sparse correlated microbiome data with Poisson-distributed true abundances and flow cytometry data
#'
#' This function generates simulated microbiome count data with user-defined 
#' sparsity in positive and negative correlations between taxa. The true abundances 
#' are modeled using latent variables from a multivariate normal distribution, 
#' and the correlations between taxa can be controlled by the `sparsity`.
#'
#' The simulated true abundances can be resampled into sequencing counts using a multinomial distribution.
#' Additionally, flow cytometry measurements for total cell counts can be simulated with user-specified standard deviations.
#'
#' @param n_taxa Integer. The number of taxa (features) to simulate.
#' @param n_samples Integer. The number of samples per condition (Pre and Post) to simulate.
#' @param rare_pct Numeric. The proportion of rare taxa (between 0 and 1).
#' @param medium_pct Numeric. The proportion of medium-abundance taxa (between 0 and 1).
#' @param seq_depth Integer. The total sequencing depth for resampling the simulated data.
#' @param sparsity Numeric. The proportion of taxa pairs to introduce positive or negative correlations (from 0 to 100, where 0 means no correlations and 100 means fully correlated).
#' @param flow_sd Numeric. The standard deviation for flow cytometry measurements (default = 300).
#' @param replicates Integer. Number of replicates for each flow cytometry sample (default = 1).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{W}{Matrix of true abundances for both Pre and Post conditions.}
#'   \item{W.para}{Matrix of proportions of taxa in each sample.}
#'   \item{W.perp}{Vector of total abundances (per sample) for all samples.}
#'   \item{Y}{Matrix of resampled sequencing counts using a multinomial distribution.}
#'   \item{rdat_flow}{Flow cytometry measurements with replicates.}
#'   \item{rdat_flow_collapse}{Flow cytometry data collapsed to mean and standard deviation (only if replicates > 1).}
#'   \item{true_abundances}{Dataframe of true abundances with Pre/Post condition labels.}
#'   \item{corr_matrix}{The correlation matrix used for simulating latent variables.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' simulated_data <- prism.simulate_prepost(
#'   n_taxa = 20, n_samples = 50, rare_pct = 0.2, medium_pct = 0.3, seq_depth = 1000, sparsity = 50
#' )
#' }
#'
#' @export
# Load necessary libraries
prism.simulate_prepost <- function(
  n_taxa = 20,
  n_samples = 1000,
  rare_pct = 0.2,
  medium_pct = 0.4,
  seq_depth = 10000,
  sparsity = 20,  # Percentage of non-zero correlations
  flow_sd = 0.55,
  replicates = 1,
  post_scale_factor = 0.8,
  taxa_index = 8,
  total_abundance_scale = 1e7,
  controlledcounts = FALSE,
  iterations = 1000,
  seed = NULL
) {
  
  # Load necessary libraries
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    install.packages("Matrix")
  }
  
  library(Matrix)     # For sparse matrices
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## 1. Assign Taxa Categories and use total abundance scale specified
  if (controlledcounts) {
    # Define the base vector d
    d <- c(4000, 4000, 4000, 4000, 4000, 400, 400, 400, 400, 4000, 400, 500, 500, 500, 400, 400, 400, 400, 400, 400)
    taxa_means_pre <- log(approx(x = seq_along(d), y = d, n = n_taxa)$y)
    log_scale_sds <- rep(0.01, n_taxa)
    cat("Adjust sequencing depth accordingly between ranges [100:5000]\n")
  } else {
    n_rare <- round(n_taxa * rare_pct)
    n_medium <- round(n_taxa * medium_pct)
    n_frequent <- n_taxa - n_rare - n_medium
    
    rare_taxa_indices <- 1:n_rare
    medium_taxa_indices <- (n_rare + 1):(n_rare + n_medium)
    frequent_taxa_indices <- (n_rare + n_medium + 1):n_taxa
    
    ## 2. Assign Means Based on Categories
    means_rare <- runif(n_rare, .0001, .005)
    means_medium <- runif(n_medium, .006, .07)
    means_frequent <- runif(n_frequent, .5, .7)
    
    means_rare <- means_rare * total_abundance_scale
    means_medium <- means_medium * total_abundance_scale
    means_frequent <- means_frequent * total_abundance_scale
    
    taxa_means_pre <- c(means_rare, means_medium, means_frequent)
    taxa_means_pre <- log(taxa_means_pre)
    
    # Define the standard deviations for each taxon on the log scale
    log_scale_sds <- runif(n_taxa, 0.05, 0.3)  # Adjust range based on biological expectations
  }

  # Function to create a sparse, symmetric, PSD correlation matrix using Matrix package
  create_sparse_psd_corr_matrix <- function(n_taxa, sparsity, taxa_index, corr_min = -0.9, corr_max = 0.9, max_attempts = 10000) {
    n_pairs <- n_taxa * (n_taxa - 1) / 2
    n_correlations <- ifelse(sparsity >= 100, n_pairs, round((sparsity / 100) * n_pairs))
    n_correlations <- max(n_correlations, 1)  # Ensure at least one correlation
    
    # Create all possible unique pairs (excluding diagonal)
    all_pairs <- combn(1:n_taxa, 2, simplify = FALSE)
    
    # Ensure the first pair includes taxa_index
    remaining_taxa <- setdiff(1:n_taxa, taxa_index)
    first_pair <- c(taxa_index, sample(remaining_taxa, 1))
    all_pairs <- all_pairs[!sapply(all_pairs, function(x) identical(x, first_pair))]
    
    # Initialize attempt counter
    attempt <- 1
    
    while (attempt <= max_attempts) {
      # Select (n_correlations - 1) random pairs
      if (n_correlations > 1) {
        selected_pairs <- sample(all_pairs, n_correlations - 1, replace = FALSE)
        selected_pairs <- c(list(first_pair), selected_pairs)
      } else {
        selected_pairs <- list(first_pair)
      }
      
      # Assign random correlation values within [corr_min, corr_max]
      corr_values <- runif(length(selected_pairs), min = corr_min, max = corr_max)
      
      # Create a sparse matrix in triplet form
      triplet <- do.call(rbind, selected_pairs)
      sparse_triplet <- sparseMatrix(
        i = triplet[,1],
        j = triplet[,2],
        x = corr_values,
        dims = c(n_taxa, n_taxa),
        symmetric = FALSE  # Initially not symmetric
      )
      
      # Now make the matrix symmetric
      sparse_sym <- forceSymmetric(sparse_triplet, uplo = "U")  # Use only the upper triangle
      
      # Set diagonal to 1
      diag(sparse_sym) <- 1
      
      # Check for positive semi-definiteness
      eigenvalues <- eigen(as.matrix(sparse_sym), symmetric = TRUE, only.values = TRUE)$values
      if (all(eigenvalues >= 0)) {
        # Valid PSD matrix found
        actual_nonzero <- length(corr_values)
        actual_sparsity <- (actual_nonzero / n_pairs) * 100
        cat(sprintf("PSD matrix achieved on attempt %d with sparsity %.2f%%\n", attempt, actual_sparsity))
        return(sparse_sym)
      }
      
      # Increment attempt counter
      attempt <- attempt + 1
    }

    # If maximum attempts reached without finding a PSD matrix
    message(sprintf("Failed to generate a PSD correlation matrix after %d attempts. Using nearPD to approximate.", max_attempts))
    calculate_sparsity <- function(matrix) {
      total_elements <- length(matrix)
      zero_elements <- sum(matrix == 0)
      sparsity <- (zero_elements / total_elements) * 100
      return(sparsity)
    }

    nearPD_result <- nearPD(sparse_sym, corr = TRUE, keepDiag = TRUE)
    sparsity <- calculate_sparsity(as.matrix(nearPD_result$mat))
    cat(sprintf("Sparsity of the returned matrix is: %.2f%%\n", sparsity))
    return(as.matrix(nearPD_result$mat))
                     
  }
  
  # Generate the sparse, PSD correlation matrix
  corr_matrix_sparse <- create_sparse_psd_corr_matrix(n_taxa, sparsity, taxa_index)
  
  # Convert sparse matrix to dense for covariance computation
  corr_matrix_dense <- as.matrix(corr_matrix_sparse)
  
  # log-scale covariance matrix
  # Create a diagonal matrix with log-scale standard deviations
  D <- diag(log_scale_sds)
  
  # Multiply the correlation matrix by the standard deviation matrix to get the covariance matrix
  log_cov_matrix <- D %*% corr_matrix_dense %*% D
  
  ## 4. Simulate Latent Variables
  generate_latent_variables <- function(n_samples, cov_matrix, taxa_means_pre) {
    latent_vars = MASS::mvrnorm(n_samples, mu = taxa_means_pre, Sigma = cov_matrix)

    # Step 2: Exponentiate to get Poisson rate parameters (lambda)
    lambda <- exp(latent_vars)

    # Step 3: Sample from Poisson distribution
    W_sampled_poisson <- matrix(rpois(n = n_samples * ncol(lambda), lambda = lambda), nrow = n_samples, ncol = ncol(lambda))

    return(W_sampled_poisson)
  }
  
  W_sampled <- generate_latent_variables(n_samples, log_cov_matrix, taxa_means_pre)
  
  ## 5. Exponentiate Latent Variables to Obtain Positive Values
  half_samples <- floor(n_samples / 2)
  W_pre <- W_sampled[1:half_samples, ]
  W_post <- W_sampled[(half_samples + 1):n_samples, ]
  
  ## 6. No Scaling of Pre-Treatment Data to represent the baseline abundances
  
  ## 7. Introduce Narrow Spectrum Antibiotic Effect on Post-Treatment Data
  # The specified taxa's mean is reduced by post_scale_factor
  # Other taxa are adjusted based on the correlation matrix
  
  # Apply antibiotic effect only to post-treatment 
  W_post[, taxa_index] <- W_post[, taxa_index] * post_scale_factor
  
  # Adjust other taxa based on correlation
  correlations <- corr_matrix_dense[taxa_index, ]
  correlations[taxa_index] <- 0  # Exclude self-correlation
  
  adjustment_proportion <- 1 - (1 - post_scale_factor) * correlations
  adjustment_proportion[adjustment_proportion < 0] <- 0
  
  W_post <- sweep(W_post, 2, adjustment_proportion, FUN = "*")
  
  # Record Post-treatment Means
  taxa_means_post <- colMeans(W_post)
  
  ## 8. Combine Pre and Post Data
  W <- rbind(W_pre, W_post)
  
  # Compute correlation matrix of W
  W_corr_matrix <- cor(W)
  
  Condition <- factor(rep(c("Pre", "Post"), each = half_samples), levels = c("Pre", "Post"))
  
  ## 9. Normalize to Get Relative Abundances
  W_para <- sweep(W, 1, rowSums(W), "/")
  
  ## 10. Simulate Sequencing Counts Using Multinomial Distribution
  resample_data <- function(W_para, seq_depth) {
    t(apply(W_para, 1, function(p) rmultinom(1, size = seq_depth, prob = p)))
  }
  
  Y <- resample_data(W_para, seq_depth)
  
  ## 11. Simulate Flow Cytometry Data
  flow_cytometry <- function(totals, replicates, flow_sd) {
    flow_vals <- sapply(totals, function(total) {
      rlnorm(replicates, meanlog = log(total), sdlog = flow_sd)
    })
    flow_data <- data.frame(
      sample = rep(1:length(totals), each = replicates),
      flow = as.vector(flow_vals)
    )
    return(flow_data)
  }
  
  W_perp <- rowSums(W)
  flow_data <- flow_cytometry(W_perp, replicates, flow_sd)
  
  ## 12. Collapse Flow Data if Replicates > 1
  if (replicates > 1) {
    flow_data_collapse <- flow_data %>%
      group_by(sample) %>%
      summarise(mean_flow = mean(flow), stdev_flow = sd(flow)) %>%
      ungroup()
  }
  
  # Format dataframes
  colnames(W) <- paste0("Taxa", 1:ncol(W)) 
  dummy <- as.data.frame(W)
  colnames(dummy) <- paste0("Taxa", 1:ncol(W)) 
  dummy$Condition <- Condition
  colnames(Y) <- paste0("Taxa", 1:ncol(Y)) 
  
  names(taxa_means_pre) <- paste0("Taxa", 1:n_taxa)
  names(taxa_means_post) <- paste0("Taxa", 1:n_taxa)
            
  ## 13. Compile Results
  results <- list(
    W.pre = W_pre,
    W.post = W_post,
    W = W,
    W.para = W_para,
    W.perp = W_perp,
    Y = Y,
    W.condition = dummy,
    taxa_means_pre = taxa_means_pre,
    taxa_means_post = taxa_means_post,
    flow = flow_data,
    corr_matrix = corr_matrix_sparse, 
    cov_matrix = log_cov_matrix,
    W_corr_matrix = W_corr_matrix
  )
  
  if (replicates > 1) {
    results$flow_collapse <- flow_data_collapse
  }
  
  return(results)
}


#' Optimizes sigma based on different algorithms
#' 
#' This internal function performs optimization using a variety of algorithms to minimize and maximize the sigma parameter. 
#' It supports several algorithms, including COBYLA, MMA, SLSQP, AUGLAG with SLSQP or LBFGS as the local optimizer, and 
#' a grid search approach. The function is flexible in handling multiple optimization methods and evaluates both 
#' the objective and constraint functions during optimization.
#' 
#' @param algorithm The optimization algorithm to use. Supported values are "COBYLA", "MMA", "SLSQP", "AUGLAG_SLSQP", 
#' "AUGLAG_LBFGS", and "GRID_SEARCH".
#' @param initialparameters The initial parameters for optimization.
#' @param taxa1relativesd Standard deviation of taxa 1.
#' @param taxa2relativesd Standard deviation of taxa 2.
#' @param relativecovariance The relative covariance between taxa 1 and 2.
#' @param rho_lower_bound_1 Lower bound for the first rho parameter.
#' @param rho_upper_bound_1 Upper bound for the first rho parameter.
#' @param rho_lower_bound_2 Lower bound for the second rho parameter.
#' @param rho_upper_bound_2 Upper bound for the second rho parameter.
#' @param lowerscalestdev Lower bound for the scale standard deviation.
#' @param upperscalestdev Upper bound for the scale standard deviation.
#' @param outputdirectory Optional. The directory to output grid search results.
#' @param d1 Optional. Identifier for the first dimension of the optimization process.
#' @param d2 Optional. Identifier for the second dimension of the optimization process.
#' @param s Optional. The current iteration of the optimization process.
#' 
#' @return A list containing the results of the optimization, including minimum and maximum values for sigma and corresponding parameter values.
#' @keywords internal
#' @examples
#' # Example usage (for internal purposes):
#' result <- optimize_sigma("COBYLA", initialparameters, taxa1relativesd, taxa2relativesd, relativecovariance, 
#'                          rho_lower_bound_1, rho_upper_bound_1, rho_lower_bound_2, rho_upper_bound_2, 
#'                          lowerscalestdev, upperscalestdev)
optimize_sigma <- function(algorithm, initialparameters, taxa1relativesd, taxa2relativesd, relativecovariance, 
                           rho_lower_bound_1, rho_upper_bound_1, rho_lower_bound_2, rho_upper_bound_2, 
                           lowerscalestdev, upperscalestdev, outputdirectory = NULL, d1 = NULL, d2 = NULL, s = NULL) {
  
  # Helper function to perform optimization using nloptr
  run_optimization <- function(initialparams, eval_f, eval_g_ineq, lb, ub, opts, eval_grad_f = NULL, eval_jac_g_ineq = NULL) {
    return(nloptr(
      x0 = initialparams,
      eval_f = eval_f,
      eval_grad_f = eval_grad_f,
      eval_g_ineq = eval_g_ineq,
      eval_jac_g_ineq = eval_jac_g_ineq,
      lb = lb,
      ub = ub,
      opts = opts
    ))
  }
  
  # Switch block for selecting algorithm
  result <- switch(algorithm,
                   
                   "COBYLA" = {
                     # Define objective and constraint functions
                     objective_f <- function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     constraint_f <- function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     
                     # Optimization for minimum
                     res_min <- run_optimization(
                       initialparameters,
                       eval_f = objective_f,
                       eval_g_ineq = constraint_f,
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
                     )
                     
                     # Optimization for maximum (negating objective)
                     res_max <- run_optimization(
                       initialparameters,
                       eval_f = function(params) -objective_f(params),
                       eval_g_ineq = constraint_f,
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1000000, "xtol_rel" = 1e-5)
                     )
                     
                     list(res_min = res_min, res_max = res_max)
                   },
                   
                   "MMA" = {
                     # Define objective, gradient, and constraint functions
                     objective_f <- function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     gradient_f <- function(params) gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     constraint_f <- function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     constraint_grad_f <- function(params) constraint_gradient_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance)
                     
                     # Optimization for minimum
                     res_min <- run_optimization(
                       initialparameters,
                       eval_f = objective_f,
                       eval_g_ineq = constraint_f,
                       eval_grad_f = gradient_f,
                       eval_jac_g_ineq = constraint_grad_f,
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = list("algorithm" = "NLOPT_LD_MMA", "maxeval" = 10000, "ftol_rel" = 1e-4)
                     )
                     
                     # Optimization for maximum (negating the objective and gradient)
                     res_max <- run_optimization(
                       initialparameters,
                       eval_f = function(params) -objective_f(params),
                       eval_g_ineq = constraint_f,
                       eval_grad_f = function(params) -gradient_f(params),
                       eval_jac_g_ineq = constraint_grad_f,
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = list("algorithm" = "NLOPT_LD_MMA", "maxeval" = 10000, "ftol_rel" = 1e-4)
                     )
                     
                     list(res_min = res_min, res_max = res_max)
                   },
                   
                   "SLSQP" = {
                     # Define optimization options
                     opts <- list(
                       "algorithm" = "NLOPT_LD_SLSQP",
                       "xtol_rel" = 1e-4,
                       "ftol_rel" = 1e-4,
                       "maxeval" = 10000
                     )
                     
                     # Perform the optimization to find the minimum sigma
                     res_min <- run_optimization(
                       initialparameters,
                       eval_f = function(params) objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                       eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = opts
                     )
                     
                     # Optimization for maximum (negating the objective)
                     res_max <- run_optimization(
                       initialparameters,
                       eval_f = function(params) -objective_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                       eval_g_ineq = function(params) constraint_function_wrapper(params, taxa1relativesd, taxa2relativesd, relativecovariance),
                       lb = c(rho_lower_bound_1, rho_lower_bound_2, lowerscalestdev),
                       ub = c(rho_upper_bound_1, rho_upper_bound_2, upperscalestdev),
                       opts = opts
                     )
                     
                     list(res_min = res_min, res_max = res_max)
                   },
                   
                   "AUGLAG_SLSQP" = {
                     # Find the minimum sigma using AUGLAG with SLSQP as the inner algorithm
                     res_min <- run_optimization(
                       initialparameters,
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
                     res_max <- run_optimization(
                       initialparameters,
                       eval_f = function(params) -objective_function_wrapper(params),
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
                     res_min <- run_optimization(
                       initialparameters,
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
                     res_max <- run_optimization(
                       initialparameters,
                       eval_f = function(params) -objective_function_wrapper(params),
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
                    # Define parameter steps (unchanged for now)
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
                      # Return a data frame with NA values
                      res_min <- list(
                        objective = NA,
                        solution = c(NA, NA, NA),
                        message = "GRIDSEARCH_NO_VALID_ROWS",
                        status = "GRIDSEARCH_NO_VALID_ROWS",
                        iterations = NA
                      )
                      
                      res_max <- list(
                        objective = NA,
                        solution = c(NA, NA, NA),
                        message = "GRIDSEARCH_NO_VALID_ROWS",
                        status = "GRIDSEARCH_NO_VALID_ROWS",
                        iterations = NA
                      )
                    
                      message(sprintf("No valid rows found for d1 %s : d2 %s for bootstrap %s | returning NA values", d1, d2, s))
                    
                    } else {
                      # Find min and max
                      min_sigma_row <- rpars[which.min(rpars$sigma), , drop = FALSE]
                      max_sigma_row <- rpars[which.max(rpars$sigma), , drop = FALSE]
                      
                      # Extract values from min_sigma_row
                      res_min <- list(
                        objective = ifelse(is.null(min_sigma_row$sigma), {message(sprintf("min_sigma_row$sigma is NULL for d1: %s d2: %s s: %s", d1, d2, s)); NA}, min_sigma_row$sigma),
                        solution = c(ifelse(is.null(min_sigma_row$rho1), {message(sprintf("min_sigma_row$rho1 is NULL for d1: %s d2: %s s: %s", d1, d2, s)); NA}, min_sigma_row$rho1),
                                     ifelse(is.null(min_sigma_row$rho2), {message(sprintf("min_sigma_row$rho2 is NULL for d1: %s d2: %s s: %s", d1, d2, s)); NA}, min_sigma_row$rho2),
                                     ifelse(is.null(min_sigma_row$scalestdevstep), {message(sprintf("min_sigma_row$scalestdevstep is NULL for d1: %s d2: %s s: %s", d1, d2, s)); NA}, min_sigma_row$scalestdevstep)),
                        message = "GRIDSEARCH_SUCCESS",
                        status = "GRIDSEARCH_SUCCESS",
                        iterations = nrow(pars),
                        spsd = nrow(rpars)
                      )
                      
                      # Extract values from max_sigma_row
                      res_max <- list(
                        objective = ifelse(is.null(max_sigma_row$sigma), {message(sprintf("max_sigma_row$sigma is NULL for d1 %s : d2 %s for bootstrap %s", d1, d2, s)); NA}, -max_sigma_row$sigma),  # Negative for maximization
                        solution = c(ifelse(is.null(max_sigma_row$rho1), {message(sprintf("max_sigma_row$rho1 is NULL for d1 %s : d2 %s for bootstrap %s", d1, d2, s)); NA}, max_sigma_row$rho1),
                                     ifelse(is.null(max_sigma_row$rho2), {message(sprintf("max_sigma_row$rho2 is NULL for d1 %s : d2 %s for bootstrap %s", d1, d2, s)); NA}, max_sigma_row$rho2),
                                     ifelse(is.null(max_sigma_row$scalestdevstep), {message(sprintf("max_sigma_row$scalestdevstep is NULL for d1 %s : d2 %s for bootstrap %s", d1, d2, s)); NA}, max_sigma_row$scalestdevstep)),
                        message = "GRIDSEARCH_SUCCESS",
                        status = "GRIDSEARCH_SUCCESS",
                        iterations = nrow(pars),
                        spsd = nrow(rpars)
                      )
                      
                      message(sprintf("Valid rows found, returning min and max results for d1 %s : d2 %s for bootstrap %s", d1, d2, s))
                    }
                    list(res_min = res_min, res_max = res_max)
                  },
                  stop("Invalid algorithm selected: ", algorithm)  # Default case if no match is found
  )
  
  return(result)
}


#' Estimate Rho and Scale Standard Deviation (SD) for Taxa
#' 
#' This internal function estimates the correlation (rho) and scale standard deviation (SD) of external scale measurements 
#' with taxa using a bootstrap-based approach. It performs this estimation in parallel and computes confidence intervals 
#' for both rho and SD based on Fisher Z-transformation for rho and chi-squared distribution for variance.
#' 
#' @param externalscalemeasurements A matrix of external scale measurements. The number of columns should match the number of rows in `Y`.
#' @param Y A matrix representing the data for taxa.
#' @param S The number of bootstrap samples.
#' @param D The number of taxa.
#' @param rWparaoriginal A 3D array representing taxa-specific values across bootstrap samples.
#' @param bootstrap_samples A matrix where each column corresponds to the indices of bootstrap samples.
#' @param alpha The significance level for confidence intervals (e.g., 0.05 for 95% confidence intervals).
#' @param prefix A string used as the prefix for filenames when saving the plot results.
#' 
#' @return A list with two components:
#' \describe{
#'   \item{scalestdev}{A matrix of scale standard deviations [S x 2], where each row contains the lower and upper bounds.}
#'   \item{rhobounds}{A 3D array of rho confidence bounds [D x 2 x S], representing the lower and upper bounds for each taxa across bootstrap samples.}
#' }
#' @keywords internal
#' @examples
#' # Example usage (for internal purposes):
#' result <- estimate_rho_and_sd(externalscalemeasurements, Y, S, D, rWparaoriginal, bootstrap_samples, alpha = 0.05, prefix = "output")
estimate_rho_and_sd <- function(externalscalemeasurements, Y, S, D, rWparaoriginal, bootstrap_samples, alpha, prefix) {
  if (!is.null(externalscalemeasurements) && is.matrix(externalscalemeasurements) && ncol(Y) == nrow(externalscalemeasurements)) {
    
    # Perform rho and SD estimation using parallel processing
    rhoandsd_list <- foreach(s = 1:S, .packages = c('stats')) %dopar% {
      n <- length(externalscalemeasurements)
      sample_indices <- bootstrap_samples[, s]
      
      # Variance of the external scale measurements
      S2 <- var(externalscalemeasurements[sample_indices])
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
        r <- cor(rWparaoriginal[taxa, sample_indices, s], externalscalemeasurements[sample_indices])
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
    
    # Clean up
    rm(rhoandsd_list)
    
    # Plot the scale SD
    prism.scalesdhistogram(scalestdev, S, filename = prefix)
    
    # Plot the rho
    prism.rhoridges(rhobounds, D, S, filename = prefix)
    
    return(list(scalestdev = scalestdev, rhobounds = rhobounds))
    
  } else {
    return(NULL)
  }
}



