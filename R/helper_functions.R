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
#' @export
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
#' @export
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
#' @export
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

#' Calculate Comprehensive Summary Statistics
#'
#' This function calculates a comprehensive set of summary statistics (mean, standard deviation, variance, median, minimum, 1st quartile, 3rd quartile, and maximum) for each numeric column in a dataframe, grouped by a specified column. Optionally, the summary statistics can be saved as a CSV file.
#' This is meant to provide summary statistics for the bootstrap results of the estimate_covariance() function in PRISM.
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
#'
calculate_bootstrap_summary <- function(data, group_col = "comparison", exclude_cols = c("d1", "d2", "s"), save_as_csv = FALSE, csv_path = "bootstrap_summary.csv") {
  
  # Ensure the group column and exclude columns are character vectors
  group_col <- as.character(group_col)
  exclude_cols <- as.character(exclude_cols)
  
  # Select columns to include in the summary statistics calculation
  data <- data %>%
    select(-all_of(exclude_cols))
  
  # Define a custom summary function to include a comprehensive set of summary statistics
  summary_stats_fn <- function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      variance = var(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      q1 = quantile(x, probs = 0.25, na.rm = TRUE),
      q3 = quantile(x, probs = 0.75, na.rm = TRUE),
      max = max(x, na.rm = TRUE))
  }
  
  # Create summary statistics grouped by the specified column
  summary_stats <- data %>%
    dplyr::group_by(across(all_of(group_col))) %>%
    dplyr::summarise(across(dplyr::where(is.numeric), summary_stats_fn, .names = "{col}_{fn}")) %>%
    tidyr::pivot_longer(cols = -all_of(group_col), 
                        names_to = c(".value", "stat"), 
                        names_sep = "_") %>%
    dplyr::arrange(across(all_of(group_col)))
  
  # Save the summary statistics as a CSV file if required
  if (save_as_csv) {
    write.csv(summary_stats, csv_path, row.names = FALSE)
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
