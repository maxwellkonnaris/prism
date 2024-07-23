#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param sigma_values A data frame containing sigma values and corresponding parameters.
#' @param save_format A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @return Plots visualizing the impact of the parameters on sigma.
#' @import ggplot2
#' @examples
#' # Example usage:
#' # Assuming 'results' is the output from run_analysis function
#' # results <- run_analysis(Y)
#' # sigma_values <- results$sigma_values
#' # sigmaplot(sigma_values)  # Without saving
#' # sigmaplot(sigma_values, save_format = "png")  # Save as PNG
#' @export
sigmaplot <- function(sigma_values, save_format = NULL) {
  comparisons <- unique(sigma_values$comparison)
  
  for (comparison in comparisons) {
    subset_data <- sigma_values[sigma_values$comparison == comparison, ]
    
    # Find the minimum and maximum sigma points
    min_sigma_point <- subset_data[which.min(subset_data$sigma), ]
    max_sigma_point <- subset_data[which.max(subset_data$sigma), ]
    
    # Plot Sigma vs rho1
    p1 <- ggplot(subset_data, aes(x = rho1, y = sigma)) +
      geom_point(alpha = 0.5) +
      geom_point(data = min_sigma_point, aes(x = rho1, y = sigma), color = "blue", size = 3, shape = 17) +
      geom_point(data = max_sigma_point, aes(x = rho1, y = sigma), color = "red", size = 3, shape = 17) +
      labs(title = paste("Sigma vs rho1 for Comparison", comparison), x = "rho1", y = "Sigma") +
      theme_minimal()
    
    # Plot Sigma vs rho2
    p2 <- ggplot(subset_data, aes(x = rho2, y = sigma)) +
      geom_point(alpha = 0.5) +
      geom_point(data = min_sigma_point, aes(x = rho2, y = sigma), color = "blue", size = 3, shape = 17) +
      geom_point(data = max_sigma_point, aes(x = rho2, y = sigma), color = "red", size = 3, shape = 17) +
      labs(title = paste("Sigma vs rho2 for Comparison", comparison), x = "rho2", y = "Sigma") +
      theme_minimal()
    
    # Plot Sigma vs x
    p3 <- ggplot(subset_data, aes(x = x, y = sigma)) +
      geom_point(alpha = 0.5) +
      geom_point(data = min_sigma_point, aes(x = x, y = sigma), color = "blue", size = 3, shape = 17) +
      geom_point(data = max_sigma_point, aes(x = x, y = sigma), color = "red", size = 3, shape = 17) +
      labs(title = paste("Sigma vs x for Comparison", comparison), x = "x", y = "Sigma") +
      theme_minimal()
    
    # Save plots to files if save_format is specified
    if (!is.null(save_format)) {
      ggsave(filename = paste0("Sigma_vs_rho1_", comparison, ".", save_format), plot = p1, width = 10, height = 6, dpi = 300)
      ggsave(filename = paste0("Sigma_vs_rho2_", comparison, ".", save_format), plot = p2, width = 10, height = 6, dpi = 300)
      ggsave(filename = paste0("Sigma_vs_x_", comparison, ".", save_format), plot = p3, width = 10, height = 6, dpi = 300)
    }
  }
}
