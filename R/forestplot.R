#' Forest Plot of Confidence Intervals
#'
#' This function creates a forest plot of confidence intervals ordered by the largest range.
#'
#' @param data A data frame containing the covariance comparison names, lower bound for the range of sigma, upper bound of the range of sigma, 2.5 percent quartile of the sorted minimized sigmas, and the 97.5 percent quartile of the sorted maximized sigmas.
#' @param bg A character string indicating the background color of the plot. Options are "white" (default) or "transparent".
#' @param save A character string indicating the file format to save the plot. Options are "png", "jpg", "svg", "pdf". Default is NULL, which means the plot is not saved.
#' @param filename A character string indicating the file name when saving the plot. Default is NULL, which means the plot is saved as forest_plot if save format is indicated.
#' @param dir_path A character string indicating the directory to store the plot. Default is \code{"./plots/"} which creates the plots directory in the current directory.
#' @return A ggplot object representing the forest plot.
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples
#' # Example usage:
#' results <- data.frame(
#'   comparison = c("A:B", "A:C", "B:C"),
#'   ninetyfive_ci_lower = c(0.1, 0.2, 0.3),
#'   ninetyfive_ci_upper = c(0.4, 0.5, 0.6),
#'   minsigma_absolute_minimum_covariance = c(0.05, 0.15, 0.25),
#'   maxsigma_absolute_maximum_covariance = c(0.45, 0.55, 0.65),
#'   p_value = c(0.01, 0.05, 0.10), # Example p-values
#'   bonferroni_p_value = c(0.03, 0.15, 0.30), # Example adjusted p-values
#'   bh_p_value = c(0.02, 0.10, 0.25) # Example adjusted p-values
#' )
#' plot <- forest_plot(results, save = "png", filename = "sampledataset")
#' plot <- forest_plot(results, save = "jpg", filename = "sampledataset")
#' plot <- forest_plot(results, save = "svg", filename = "sampledataset")
#' plot <- forest_plot(results, save = "pdf", filename = "sampledataset")
forest_plot <- function(data, bg = "white", save = NULL, filename = NULL, dir_path = "./plots/") {
  # Ensure the data has the necessary columns
  if (!all(c("comparison", "ninetyfive_ci_lower", "ninetyfive_ci_upper", "minsigma_absolute_minimum_covariance", "maxsigma_absolute_maximum_covariance", "p_value", "bonferroni_p_value", "bh_p_value") %in% colnames(data))) {
    stop("Data must contain 'comparison', 'ninetyfive_ci_lower', 'ninetyfive_ci_upper', 'minsigma_absolute_minimum_covariance', 'maxsigma_absolute_maximum_covariance', 'p_value', 'bonferroni_p_value', and 'bh_p_value' columns")
  }

  # Check if the directory exists
  if (!dir.exists(dir_path)) {
    # Create the directory
    dir.create(dir_path)
    if (dir.exists(dir_path)) {
      cat("Directory created successfully!")
    } else {
      cat("Failed to create directory.")
    }
  } else {
    cat("Directory already present")
  }
  
  # Calculate the range of the confidence intervals
  data <- data %>%
    mutate(range = ninetyfive_ci_upper - ninetyfive_ci_lower)
  
  # Reorder the comparison names by the range
  data$comparison <- factor(data$comparison, levels = data$comparison[order(data$range, decreasing = TRUE)])
  
  # Highlight intervals that do not cover 0
  data <- data %>%
    mutate(highlight = ifelse((ninetyfive_ci_lower > 0 & ninetyfive_ci_upper > 0) | (ninetyfive_ci_lower < 0 & ninetyfive_ci_upper < 0), "95% CI Doesnt Cover Zero", "95% CI Covers Zero"))
  
  # Create the forest plot with confidence intervals and sigma ranges
  plot <- ggplot(data, aes(x = comparison)) +
    coord_flip() +
    theme_minimal(base_size = 15) +
    labs(
      title = "Covariance Intervals",
      x = "Taxa Comparison",
      y = "Estimated Covariance/Variance Range (log scale)",
      color = "Legend",
      size = "p-value"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15)
    ) +
    scale_color_manual(values = c("95% CI Covers Zero" = "#FF00FF", "95% CI Doesnt Cover Zero" = "green", "Covariance Range" = "#000000")) + 
    scale_size_continuous(range = c(1, 10), breaks = c(1, 2, 3), labels = c("0.1", "0.01", "0.001"))
  
  # Add the range from minimum to maximum sigma value
  plot <- plot + geom_errorbar(aes(ymin = minsigma_absolute_minimum_covariance, ymax = maxsigma_absolute_maximum_covariance, color = "Covariance Range"), width = 0.5, size = 1)

  # Add points for the p-values
  plot <- plot + geom_point(aes(y = (ninetyfive_ci_lower + ninetyfive_ci_upper) / 2, size = -log10(p_value)), color = "black")

  # Add the 95% confidence intervals and highlight those not covering 0
  plot <- plot + geom_errorbar(aes(ymin = ninetyfive_ci_lower, ymax = ninetyfive_ci_upper, color = highlight), width = 0.5, size = 1)
  
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
    file_name <- paste0(dir_path,"forest_plot.", save)
    if (!is.null(filename)) {
      file_name <- paste0(dir_path, filename, ".", save)
    }
    ggsave(file_name, plot, width = 12, height = 15, dpi = 300, device = save)
  }
  
  return(plot)
}
