#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param all_inner_results A data frame containing sigma values and corresponding parameters.
#' @param save_format A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @param plot_type A character string specifying the type of plot ("point" or "line"). Defaults to "point".
#' @return Plots visualizing the impact of the parameters on sigma.
#' @import ggplot2
#' @import gridExtra
#' @examples
#' # Example usage:
#' # Assuming 'results' is the output from run_analysis function
#' # results <- run_analysis(Y)
#' # all_inner_results <- results$all_inner_results
#' # sigmaplot(all_inner_results)  # Without saving
#' # sigmaplot(all_inner_results, save = "png", plot_type = "line")  # Save as PNG with best fit lines
#' @export
sigmaplot <- function(all_inner_results, filename = NULL, save = NULL, plot_type = "line") {
  
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

    # Ensure 'comparison' column exists
  if (!"comparison" %in% colnames(all_inner_results)) {
    all_inner_results$comparison <- "default"  # Add a default value if the column doesn't exist
  }
  
  # Custom theme to adjust legend position and size
  custom_theme <- theme_minimal() + 
    theme(legend.position = "bottom", 
          legend.key.size = unit(0.5, "lines"), 
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  # Function to create plots
  create_plot <- function(xvar, yvar, comparison, title) {
    if (plot_type == "line") {
      p <- ggplot(all_inner_results, aes_string(x = xvar, y = yvar, color = comparison)) +
        geom_smooth(method = "lm", se = TRUE) +  
        labs(title = title, x = xvar, y = yvar) +
        custom_theme
    } else {
      p <- ggplot(all_inner_results, aes_string(x = xvar, y = yvar, color = comparison)) +
        geom_point(shape = 16) +  # Use circle shape for points
        labs(title = title, x = xvar, y = yvar) +
        custom_theme
    }
    return(p)
  }
  
  # Min and Max Sigma vs rho1
  p1_min <- create_plot("min_rho1", "minsigma", "comparison", "Min Sigma vs rho1")
  p1_max <- create_plot("max_rho1", "maxsigma", "comparison", "Max Sigma vs rho1")
  
  # Min and Max Sigma vs rho2
  p2_min <- create_plot("min_rho2", "minsigma", "comparison", "Min Sigma vs rho2")
  p2_max <- create_plot("max_rho2", "maxsigma", "comparison", "Max Sigma vs rho2")
  
  # Min and Max Sigma vs x
  p3_min <- create_plot("min_x", "minsigma", "comparison", "Min Sigma vs x")
  p3_max <- create_plot("max_x", "maxsigma", "comparison", "Max Sigma vs x")
  
  # # Combined
  # combined_plot <- gridExtra::grid.arrange(p1_min, p1_max, p2_min, p2_max, p3_min, p3_max, ncol = 2)
  
  # # Display the plots with combined legend at the bottom
  # gridExtra::grid.arrange(p1_min + theme(legend.position="none"), 
  #                         p1_max + theme(legend.position="none"), 
  #                         p2_min + theme(legend.position="none"), 
  #                         p2_max + theme(legend.position="none"), 
  #                         p3_min + theme(legend.position="none"), 
  #                         p3_max + theme(legend.position="none"), 
  #                         ncol = 2,
  #                         top = ggplot2::textGrob("Sigma Plots", gp = grid::gpar(fontsize = 16, fontface = "bold")),
  #                         bottom = ggplot2::textGrob("Combined Legend",
  #                                                    gp = grid::gpar(fontsize = 10, fontface = "bold"),
  #                                                    vp = grid::viewport(y = unit(1, "npc"))))
  
  # Save plots to files if save_format is specified
  if (!is.null(save)) {
    ggsave(filename = paste0(filename,"Min_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0(filename,"Max_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_max, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0(filename,"Min_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0(filename,"Max_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_max, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0(filename,"Min_Sigma_vs_x_", plot_type, ".", save), plot = p3_min, width = 10, height = 10, dpi = 300)
    ggsave(filename = paste0(filename,"Max_Sigma_vs_x_", plot_type, ".", save), plot = p3_max, width = 10, height = 10, dpi = 300)
    #ggsave(filename = paste0(filename, "Combined_", plot_type, ".", save), plot = combined_plot, width = 15, height = 30, dpi = 300)
  }
}
