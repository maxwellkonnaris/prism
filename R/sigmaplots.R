#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param all_inner_results A data frame containing sigma values and corresponding parameters.
#' @param save_format A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @param plot_type A character string specifying the type of plot ("point" or "average"). Defaults to "point".
#' @return Plots visualizing the impact of the parameters on sigma.
#' @import ggplot2
#' @import gridExtra
#' @examples
#' # Example usage:
#' # Assuming 'results' is the output from run_analysis function
#' # results <- run_analysis(Y)
#' # all_inner_results <- results$all_inner_results
#' # sigmaplot(all_inner_results)  # Without saving
#' # sigmaplot(all_inner_results, save_format = "png", plot_type = "average")  # Save as PNG with average lines
#' @export
sigmaplot <- function(all_inner_results, save_format = NULL, plot_type = "point") {
  
  # Convert character columns to numeric
  all_inner_results$s <- as.numeric(all_inner_results$s)
  all_inner_results$minsigma <- as.numeric(all_inner_results$minsigma)
  all_inner_results$maxsigma <- as.numeric(all_inner_results$maxsigma)
  all_inner_results$min_rho1 <- as.numeric(all_inner_results$min_rho1)
  all_inner_results$min_rho2 <- as.numeric(all_inner_results$min_rho2)
  all_inner_results$min_x <- as.numeric(all_inner_results$min_x)
  all_inner_results$max_rho1 <- as.numeric(all_inner_results$max_rho1)
  all_inner_results$max_rho2 <- as.numeric(all_inner_results$max_rho2)
  all_inner_results$max_x <- as.numeric(all_inner_results$max_x)
  
  # Custom theme to adjust legend position and size
  custom_theme <- theme_minimal() + 
    theme(legend.position = "right", 
          legend.key.size = unit(0.5, "lines"), 
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  # Function to create plots
  create_plot <- function(xvar, yvar, title) {
    if (plot_type == "average") {
      p <- ggplot(all_inner_results, aes_string(x = xvar, y = yvar, color = "comparison", fill = "comparison")) +
        stat_summary(aes_string(group = "comparison"), 
                     fun.data = mean_se, 
                     geom = "ribbon", 
                     alpha = 0.2) +  # Error zone
        stat_summary(aes_string(group = "comparison"), 
                     fun = mean, 
                     geom = "line", 
                     size = 1.5) +  # Thicker line
        labs(title = title, x = xvar, y = yvar) +
        custom_theme
    } else {
      p <- ggplot(all_inner_results, aes_string(x = xvar, y = yvar, color = "comparison")) +
        geom_point(shape = 16) +
        labs(title = title, x = xvar, y = yvar) +
        custom_theme
    }
    return(p)
  }
  
  # Min and Max Sigma vs rho1
  p1_min <- create_plot("min_rho1", "minsigma", "Min Sigma vs rho1")
  p1_max <- create_plot("max_rho1", "maxsigma", "Max Sigma vs rho1")
  
  # Min and Max Sigma vs rho2
  p2_min <- create_plot("min_rho2", "minsigma", "Min Sigma vs rho2")
  p2_max <- create_plot("max_rho2", "maxsigma", "Max Sigma vs rho2")
  
  # Min and Max Sigma vs x
  p3_min <- create_plot("min_x", "minsigma", "Min Sigma vs x")
  p3_max <- create_plot("max_x", "maxsigma", "Max Sigma vs x")
  
  # Combined
  combined_plot <- gridExtra::grid.arrange(p1_min, p1_max, p2_min, p2_max, p3_min, p3_max, ncol = 2)
  
  # Display the plots
  gridExtra::grid.arrange(p1_min, p1_max, p2_min, p2_max, p3_min, p3_max, ncol = 2)
  
  # Save plots to files if save_format is specified
  if (!is.null(save_format)) {
    ggsave(filename = paste0("Min_Sigma_vs_rho1.", save_format), plot = p1_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Max_Sigma_vs_rho1.", save_format), plot = p1_max, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Min_Sigma_vs_rho2.", save_format), plot = p2_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Max_Sigma_vs_rho2.", save_format), plot = p2_max, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Min_Sigma_vs_x.", save_format), plot = p3_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Max_Sigma_vs_x.", save_format), plot = p3_max, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0("Combined.", save_format), plot = combined_plot, width = 15, height = 25, dpi = 300)
  }
}
