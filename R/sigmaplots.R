#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param all_inner_results A data frame containing sigma values and corresponding parameters.
#' @param save_format A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @return Plots visualizing the impact of the parameters on sigma.
#' @import ggplot2
#' @import gridExtra
#' @examples
#' # Example usage:
#' # Assuming 'results' is the output from run_analysis function
#' # results <- run_analysis(Y)
#' # all_inner_results <- results$all_inner_results
#' # sigmaplot(all_inner_results)  # Without saving
#' # sigmaplot(all_inner_results, save_format = "png")  # Save as PNG
#' @export
sigmaplot <- function(all_inner_results, save_format = NULL) {
  
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
  
  # Min and Max Sigma vs rho1
  p1_min <- ggplot(all_inner_results, aes(x = min_rho1, y = minsigma, color = comparison)) +
    geom_point() +
    labs(title = "Min Sigma vs rho1", x = "rho1", y = "Min Sigma") +
    theme_minimal()
  
  p1_max <- ggplot(all_inner_results, aes(x = max_rho1, y = maxsigma, color = comparison)) +
    geom_point() +
    labs(title = "Max Sigma vs rho1", x = "rho1", y = "Max Sigma") +
    theme_minimal()
  
  # Min and Max Sigma vs rho2
  p2_min <- ggplot(all_inner_results, aes(x = min_rho2, y = minsigma, color = comparison)) +
    geom_point() +
    labs(title = "Min Sigma vs rho2", x = "rho2", y = "Min Sigma") +
    theme_minimal()
  
  p2_max <- ggplot(all_inner_results, aes(x = max_rho2, y = maxsigma, color = comparison)) +
    geom_point() +
    labs(title = "Max Sigma vs rho2", x = "rho2", y = "Max Sigma") +
    theme_minimal()
  
  # Min and Max Sigma vs x
  p3_min <- ggplot(all_inner_results, aes(x = min_x, y = minsigma, color = comparison)) +
    geom_point() +
    labs(title = "Min Sigma vs x", x = "x", y = "Min Sigma") +
    theme_minimal()
  
  p3_max <- ggplot(all_inner_results, aes(x = max_x, y = maxsigma, color = comparison)) +
    geom_point() +
    labs(title = "Max Sigma vs x", x = "x", y = "Max Sigma") +
    theme_minimal()

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
    ggsave(filename = paste0("Combined", save_format), plot = combined_plot, width = 15, heigt = 25, dpi = 300)
  }
}
