#' Forest Plot of Confidence Intervals
#'
#' This function creates a forest plot of confidence intervals ordered by the largest range.
#'
#' @param data A data frame containing the comparison names, lower confidence intervals, upper confidence intervals, and optionally minsigma and maxsigma columns.
#' @param bg A character string indicating the background color of the plot. Options are "transparent" (default) or "white".
#' @param save A character string indicating the file format to save the plot. Options are "png", "jpg", "svg", "pdf". Default is NULL, which means the plot is not saved.
#' @return A ggplot object representing the forest plot.
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples
#' # Example usage:
#' results <- data.frame(
#'   comparison = c("A:B", "A:C", "B:C"),
#'   cilower = c(0.1, 0.2, 0.3),
#'   ciupper = c(0.4, 0.5, 0.6),
#'   minsigma = c(0.05, 0.15, 0.25),
#'   maxsigma = c(0.45, 0.55, 0.65)
#' )
#' plot <- forest_plot(results, save = "png")
#' plot <- forest_plot(results, save = "jpg")
#' plot <- forest_plot(results, save = "svg")
#' plot <- forest_plot(results, save = "pdf")
forest_plot <- function(data, bg = "white", save = NULL, filename = NULL) {
  # Ensure the data has the necessary columns
  if (!all(c("comparison", "cilower", "ciupper") %in% colnames(data))) {
    stop("Data must contain 'comparison', 'cilower', and 'ciupper' columns")
  }
  
  # Check if the data has minsigma and maxsigma columns
  has_sigma <- all(c("minsigma", "maxsigma") %in% colnames(data))
  
  # Calculate the range of the confidence intervals
  data <- data %>%
    mutate(range = ciupper - cilower)
  
  # Reorder the comparison names by the range
  data$comparison <- factor(data$comparison, levels = data$comparison[order(data$range, decreasing = TRUE)])
  
  # Highlight intervals that do not cover 0
  data <- data %>%
    mutate(highlight = ifelse((cilower > 0 & ciupper > 0) | (cilower < 0 & ciupper < 0), "95% CI Doesnt Cover Zero", "95% CI Covers Zero"))
  
  # Create the forest plot with confidence intervals and sigma ranges
  plot <- ggplot(data, aes(x = comparison)) +
    coord_flip() +
    theme_minimal(base_size = 15) +
    labs(
      title = "Forest Plot of Confidence Intervals",
      x = "Taxa Comparison",
      y = "Confidence Interval of Estimated Covariance/Variance (log scale)",
      color = "Legend"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15)
    ) +
    scale_color_manual(values = c("95% CI Covers Zero" = "#FF00FF", "95% CI Doesnt Cover Zero" = "green", "Sigma Range" = "#000000"))
  
  # Add sigma ranges if available
  if (has_sigma) {
    plot <- plot + geom_errorbar(aes(ymin = minsigma, ymax = maxsigma, color = "Sigma Range"), width = 0.3, size = 1)
  }
  
  # Add the 95% confidence intervals and highlight those not covering 0
  plot <- plot + geom_errorbar(aes(ymin = cilower, ymax = ciupper, color = highlight), width = 0.3, size = 1)
  
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
  
  # Ensure there is a minimum tick mark
  plot <- plot + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  
  # Save the plot if save is not NULL
  if (!is.null(save)) {
    file_name <- paste0("forest_plot.", save)
    if (!is.null(filename)) {
      file_name <- paste0(filename, ".", save  
    }
    ggsave(file_name, plot, width = 12, height = 15, dpi = 300, device = save)
  }
  
  return(plot)
}
