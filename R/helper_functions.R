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


#' Simulate Pre and Post Abundance Data with Flow Cytometry Measurements
#'
#' This function allows users to simulate three scenarios:
#' 1. Apply mean shifts between pre and post conditions with no covariance structure.
#' 2. Apply a covariance structure with no change in means between pre and post conditions.
#' 3. Apply both mean shifts and a covariance structure together.
#'
#' @param n Integer. Number of samples to simulate. Default is `50`.
#' @param d Integer. Number of taxa to simulate. Default is `20`.
#' @param seq_depth Integer. Sequencing depth for resampling. Default is `5000`.
#' @param replicate Integer. Number of flow cytometry replicates per sample. Default is `1`.
#' @param corr_strengths Matrix. Covariance matrix for simulating correlated taxa. Default is `NULL`.
#' @param scenario Character. The scenario to simulate. Options are:
#'   \describe{
#'     \item{"means_only"}{Apply different means between pre and post conditions without a covariance structure.}
#'     \item{"covariance_only"}{Apply a covariance structure without changing the means.}
#'     \item{"both"}{Apply both different means and a covariance structure.}
#'   } Default is `"both"`.
#' @param coverage Numeric. Percentage of sequencing depth used to resample the true values. Must be between 0 and 1. Default is `1` (i.e., 100% of sequencing depth).
#' @param perfect_resolution Logical. If `TRUE`, uses near-perfect resolution, allowing multinomial resampling but making `seq_depth` very large. Default is `FALSE`.
#' @param flow_sd Numeric. Standard deviation for flow cytometry measurement error. Default is `300`.
#'
#' @return A list containing:
#' \describe{
#'   \item{dat}{Data frame of true abundances with pre and post conditions. Columns include `Condition` (Pre/Post) and taxa counts.}
#'   \item{rdat}{Data frame of resampled abundances with specified sequencing depth. Columns mirror those in `dat`.}
#'   \item{flow_data}{Data frame of flow cytometry measurements. Includes `sample` and `flow` columns.}
#'   \item{flow_data_collapse}{(Optional) Data frame of collapsed flow cytometry measurements with mean and standard deviation per sample. Included only if `replicate > 1`.}
#' }
#'
#' @details
#' - **Means Only** (`means_only`): Simulates abundance changes without covariances.
#' - **Covariance Only** (`covariance_only`): Simulates covariances with the same means in both conditions.
#' - **Both** (`both`): Simulates both changes in mean abundance and a covariance structure.
#'
#' @importFrom dplyr mutate select filter group_by ungroup
#' @importFrom purrr map
#' @importFrom stats rpois rmultinom rnorm cov
#' @export
prism.simulate_prepost <- function(n = 50, d = 20, seq_depth = 5000, replicate = 1, corr_strengths = NULL, scenario = "both", coverage = 1, perfect_resolution = FALSE, flow_sd = 300) {
  
  # Check if corr_strengths is provided and if its dimensions match the number of taxa (d)
  if (!is.null(corr_strengths)) {
    if (nrow(corr_strengths) != d || ncol(corr_strengths) != d) {
      stop("The dimensions of the covariance matrix `corr_strengths` must match the number of taxa `d`.")
    }
  }

  # Check if coverage is within the valid range
  if (coverage < 0 || coverage > 1) {
    stop("Coverage must be a numeric value between 0 and 1.")
  }

  ## Helper Function to Create Abundances
  create_abundances <- function(d_pre, d_post, n, corr_strengths = NULL, apply_covariance = FALSE) {
    if (apply_covariance && !is.null(corr_strengths)) {
      # Generate correlated log-abundances using a multivariate normal distribution
      pre_log_abundances <- rmvnorm(n, mu = log(d_pre), Sigma = corr_strengths)
      post_log_abundances <- rmvnorm(n, mu = log(d_post), Sigma = corr_strengths)

      W_pre <- exp(pre_log_abundances)
      W_post <- exp(post_log_abundances)
      
      # Convert log-abundances back to abundance counts using exp() and Poisson resampling
      pre_abundances <- matrix(rpois(n * d, W_pre), nrow = n)
      post_abundances <- matrix(rpois(n * d, W_post), nrow = n)
    } else {
      # If no covariance is applied, use independent Poisson resampling
      pre_abundances <- matrix(rpois(n * d, d_pre), nrow = n)
      post_abundances <- matrix(rpois(n * d, d_post), nrow = n)
    }
    
    # Combine into a data frame with conditions
    dat <- rbind(data.frame(Condition = "Pre", pre_abundances),
                 data.frame(Condition = "Post", post_abundances))
    
    colnames(dat)[-1] <- paste0("Taxa", 1:d)
    return(dat)
  }
  
  ## Resample Data with Coverage or Near-Perfect Resolution
  resample_data <- function(dat, seq_depth, coverage, perfect_resolution) {
    if (perfect_resolution) {
      # Set a very large seq_depth for near-perfect multinomial sampling
      seq_depth <- 1e6
    }
    
    # Apply sequencing depth and coverage
    ddat <- as.matrix(dat[,-1]) / rowSums(as.matrix(dat[,-1]))
    for (i in 1:nrow(dat)) {
      adjusted_depth <- seq_depth * coverage
      dat[i, -1] <- stats::rmultinom(1, size = adjusted_depth, prob = ddat[i, ])
    }
    return(dat)
  }
  
  ## Simulate Flow Cytometry Measurements with Specified SD
  flow_cytometry <- function(totals, samp_names, replicates, flow_sd) {
    samp_names <- rep(samp_names, each = replicates)
    flow_vals <- sapply(totals, FUN = function(total, replicates) {
      stats::rnorm(replicates, mean = total, sd = flow_sd)
    }, replicates = replicates, simplify = TRUE)
    flow_data <- data.frame("sample" = samp_names, "flow" = c(flow_vals))
    return(flow_data)
  }
  
  ### Set the Means for Pre and Post Conditions Based on Scenario
  if (scenario == "means_only") {
    d_pre <- runif(d, 400, 5000)
    d_post <- d_pre
    num_changed_taxa <- round(0.4 * d)  # Change ~40% of taxa post-condition
    d_post[1:num_changed_taxa] <- runif(num_changed_taxa, 50, 3000)
    
    # Generate true abundances (without covariance structure)
    dat <- create_abundances(d_pre, d_post, n, apply_covariance = FALSE)
    
  } else if (scenario == "covariance_only") {
    d_pre <- runif(d, 400, 5000)
    d_post <- d_pre
    
    # Generate true abundances with the covariance structure
    dat <- create_abundances(d_pre, d_post, n, corr_strengths = corr_strengths, apply_covariance = TRUE)
    
  } else if (scenario == "both") {
    d_pre <- runif(d, 400, 5000)
    d_post <- d_pre
    num_changed_taxa <- round(0.4 * d)
    d_post[1:num_changed_taxa] <- runif(num_changed_taxa, 50, 3000)

    # Generate true abundances with the covariance structure
    dat <- create_abundances(d_pre, d_post, n, corr_strengths = corr_strengths, apply_covariance = TRUE)
  } else {
    stop("Invalid scenario. Choose from 'means_only', 'covariance_only', or 'both'.")
  }
  
  ## Resample Data with Coverage or Near-Perfect Resolution
  rdat <- resample_data(dat, seq_depth = seq_depth, coverage = coverage, perfect_resolution = perfect_resolution)
  
  ## Finding Sample Totals
  totals <- rowSums(dat[,-1])
  
  ## Simulate Flow Cytometry Measurements with User-Specified SD
  flow_data <- flow_cytometry(totals, samp_names = rownames(dat), replicates = replicate, flow_sd = flow_sd)
  
  ## If Replicates > 1, Collapse Flow Data to Compute Mean and SD
  if (replicate > 1) {
    flow_data_collapse <- flow_data %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(mean = mean(flow)) %>%
      dplyr::mutate(stdev = stats::sd(flow)) %>%
      dplyr::select(-flow) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()
  }
  
  ## Compile Results
  if (replicate > 1) {
    return(list(
      dat = dat,
      rdat = rdat,
      rdat_flow = flow_data,
      rdat_flow_collapse = flow_data_collapse, 
      dat_scale = totals
    ))
  } else {
    return(list(
      dat = dat,
      rdat = rdat,
      rdat_flow = flow_data,
      dat_scale = totals
    ))
  }
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
