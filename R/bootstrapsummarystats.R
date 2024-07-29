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
    group_by(across(all_of(group_col))) %>%
    summarise(across(where(is.numeric), summary_stats_fn, .names = "{col}_{fn}"))
  
  # Save the summary statistics as a CSV file if required
  if (save_as_csv) {
    write.csv(summary_stats, csv_path, row.names = FALSE)
  }
  
  return(summary_stats)
}
