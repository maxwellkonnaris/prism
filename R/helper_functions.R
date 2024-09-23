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
  
  # Calculate p-values for each row in the results dataframe
  results$p_value <- apply(results, 1, function(row) {
    lower <- as.numeric(row["ninetyfive_ci_lower"])
    upper <- as.numeric(row["ninetyfive_ci_upper"])
    get_p_value(lower, upper)
  })
  
  # Number of hypotheses/tests
  m <- nrow(results)
  
  # Bonferroni correction
  results$bonferroni_p_value <- pmin(results$p_value * m, 1)
  
  # Benjamini-Hochberg FDR method
  results <- results[order(results$p_value), ]
  results$bh_p_value <- p.adjust(results$p_value, method = "BH")
  results <- results[order(as.numeric(rownames(results))), ]
  
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
#' @export
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
#' @export
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
#' summary_stats <- calculate_summary_stats(results$all_inner_results, group_col = "comparison", exclude_cols = c("d1", "d2"), save_as_csv = TRUE, csv_path = "bootstrap_summary_stats.csv")
#' print(summary_stats)
calculate_bootstrap_summary <- function(data, group_col = "comparison", exclude_cols = c("d1", "d2", "s"), save_as_csv = FALSE, csv_path = "bootstrap_summary.csv") {
  
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


#' Simulate Data
#'
#' This function simulates data for given dimensions and sample size, with specified sequencing depth.
#'
#' @param D An integer specifying the number of dimensions.
#' @param N An integer specifying the number of samples.
#' @param seq.depth A numeric value specifying the sequencing depth.
#' @return A list containing the simulated data matrix \code{Y}, the log-transformed matrix \code{logW}, and the covariance matrix \code{Sigma}.
#' @examples
#' # Example usage:
#' simulated_data <- simulate_data(3, 100, 1000)
#' @export
simulate_data <- function(D, N, seq.depth) {
  Sigma <- diag(D)
  Sigma[2, 1] <- Sigma[1, 2] <- -0.5
  Sigma[3, 1] <- Sigma[1, 3] <- +0.5

  logW <- rmvnorm(N, rep(0, D), Sigma)
  W <- exp(logW)
  Wpara <- t(miniclo(t(W)))
  Y <- Wpara * seq.depth

  return(list(Y = Y, logW = logW, Sigma = Sigma))
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
#' rho_values <- calculate_proportionality_metrics(relative_abundances)
#' print(rho_values)
#' @export
calculate_proportionality <- function(data) {
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
    lock <- filelock::lock(lock_file)
    
    write.table(results, file = pair_file_name, append = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }, error = function(e) {
    message("Error while writing to file: ", pair_file_name, "\n", e)
  }, finally = {
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
