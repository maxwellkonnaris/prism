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
  required_columns <- c("cilower", "ciupper", "comparison")
  if (!all(required_columns %in% colnames(results))) {
    stop("The input results dataframe must contain the following columns: cilower, ciupper, comparison")
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
    lower <- as.numeric(row["cilower"])
    upper <- as.numeric(row["ciupper"])
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

