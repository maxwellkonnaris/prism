#' Forest Plot of Confidence Intervals
#'
#' This function creates a forest plot of confidence intervals ordered by the largest range.
#'
#' @param data A data frame containing the comparison names, lower confidence intervals, and upper confidence intervals.
#' @return A ggplot object representing the forest plot.
#' @import ggplot2
#' @import dplyr
#' @export
forest_plot <- function(data, bg=transparent) {
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
    geom_errorbar(width = 0.3, size = 1) +
    coord_flip() +
    theme_minimal(base_size = 15) +
    labs(
      title = "Forest Plot of Confidence Intervals",
      x = "Taxa Comparison",
      y = "Confidence Interval of Estimated Covariance/Variance (log scale)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15)
    )

  # Customize the background based on the bg parameter
  if (bg == "transparent") {
    # Do nothing as the default is already minimal with transparent background
  } else if (bg == "white") {
    plot <- plot + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
  } else {
    stop("bg parameter must be 'transparent' or 'white'")
  }
  
  return(plot)
}

# Example usage (uncomment and run when using actual data)
# results <- data.frame(
#   comparison = c("A:B", "A:C", "B:C"),
#   cilower = c(0.1, 0.2, 0.3),
#   ciupper = c(0.4, 0.5, 0.6)
# )
# plot <- forest_plot(results)
# ggsave("forest_plot.png", plot, width = 12, height = 8, dpi = 300)
