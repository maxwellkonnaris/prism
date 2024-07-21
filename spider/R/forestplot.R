#' Forest Plot of Confidence Intervals
#'
#' This function creates a forest plot of confidence intervals ordered by the largest range.
#'
#' @param data A data frame containing the comparison names, lower confidence intervals, and upper confidence intervals.
#' @return A ggplot object representing the forest plot.
#' @examples
#' \dontrun{
#'   results <- data.frame(
#'     comparison = c("A:B", "A:C", "B:C"),
#'     cilower = c(0.1, 0.2, 0.3),
#'     ciupper = c(0.4, 0.5, 0.6)
#'   )
#'   forest_plot(results)
#' }
#' @import ggplot2
#' @import dplyr
#' @export
forest_plot <- function(data) {
  # Ensure the data has the necessary columns
  if (!all(c("comparison", "cilower", "ciupper") %in% colnames(data))) {
    stop("Data must contain 'comparison', 'cilower', and 'ciupper' columns")
  }
  
  # Calculate the range of the confidence intervals
  data <- data %>%
    mutate(range = ciupper - cilower)
  
  # Reorder the comparison names by the range
  data$comparison <- factor(data$comparison, levels = data$comparison[order(data$range, decreasing = TRUE)])
  
  # Create the forest plot with only confidence intervals
  plot <- ggplot(data, aes(x = comparison, ymin = cilower, ymax = ciupper)) +
    geom_errorbar(width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Forest Plot of Confidence Intervals",
      x = "Taxa Comparison",
      y = "Confidence Interval of True Covariance/Variance (Log Scale)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(plot)
}
