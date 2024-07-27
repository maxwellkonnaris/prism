#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param all_inner_results A data frame containing sigma values and corresponding parameters.
#' @param save A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @param plot_type A character string specifying the type of plot ("point" or "line"). Defaults to "line".
#' @param filename A character string specifying an alternate name of the file. Defaults to \code{NULL} which saves as parameter combination names.
#' @param individual A boolean type specifying whether you would like to save individual plots. Defaults to \code{FALSE}.
#' @param bg A character string indicating the background color of the plot. Options are "white" (default) or "transparent".
#' @return Plots visualizing the impact of the parameters on sigma.
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @examples
#' # Example usage:
#' # Assuming 'results' is the output from estimate_covariance function
#' # results <- estimate_covariance(Y)
#' # all_inner_results <- results$all_inner_results # OR 
#' # all_inner_results <- results$all_inner_results_moment
#' # sigmaplot(all_inner_results)  # Without saving
#' # sigmaplot(all_inner_results, save = "png", filename = "sampledataset", individual = FALSE, plot_type = "line")  # Save as PNG with best fit lines
#' @export
sigmaplot <- function(all_inner_results, bg="white", filename = NULL, save = NULL, individual = FALSE, plot_type = "line") {
  
  # Convert relevant columns to numeric
  all_inner_results <- all_inner_results %>%
    mutate(across(c(s, minsigma, maxsigma, min_rho1, min_rho2, min_x, max_rho1, max_rho2, max_x), as.numeric))

    # Ensure 'comparison' column exists
  if (!"comparison" %in% colnames(all_inner_results)) {
    all_inner_results$comparison <- "default"  # Add a default value if the column doesn't exist
  }
  all_inner_results$comparison <- as.factor(all_inner_results$comparison)
  
  
    # Customize the background based on the bg parameter
  if (bg == "transparent") {
    	# Custom theme to adjust legend position and size
	  custom_theme <- theme_minimal() + 
	    theme(legend.position = "bottom", 
		  legend.key.size = unit(0.5, "lines"), 
		  legend.text = element_text(size = 8),
		  plot.title = element_text(size = 14, face = "bold"),
		  axis.title = element_text(size = 12),
		  axis.text = element_text(size = 10))
  } else if (bg == "white") {
    	# Custom theme to adjust legend position and size
	  custom_theme <- theme_minimal() + 
	    theme(legend.position = "bottom", 
		  legend.key.size = unit(0.5, "lines"), 
		  legend.text = element_text(size = 8),
		  plot.background = element_rect(fill = "white", color = NA),
		  panel.background = element_rect(fill = "white", color = NA),
		  plot.title = element_text(size = 14, face = "bold"),
		  axis.title = element_text(size = 12),
		  axis.text = element_text(size = 10))
  } else {
    stop("bg parameter must be 'transparent' or 'white'")
  }
  
  
  # Function to create plots
  create_plot <- function(xvar, yvar, comparison, title) {
    if (plot_type == "line") {
      p <- ggplot(all_inner_results, aes_string(x = xvar, y = yvar, color = comparison)) +
      	geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +  
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
  
  # Combined
  combined_plot <- gridExtra::grid.arrange(p1_min, p1_max, p2_min, p2_max, p3_min, p3_max, ncol = 2)
  
  # Display the plots with combined legend at the bottom
  gridExtra::grid.arrange(p1_min + theme(legend.position="none"), 
                          p1_max + theme(legend.position="none"), 
                          p2_min + theme(legend.position="none"), 
                          p2_max + theme(legend.position="none"), 
                          p3_min + theme(legend.position="none"), 
                          p3_max + theme(legend.position="none"), 
                          ncol = 2,
                          top = grid::textGrob("Sigma Plots", gp = grid::gpar(fontsize = 16, fontface = "bold")),
                          bottom = grid::textGrob("Combined Legend",
                                                     gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                                     vp = grid::viewport(y = unit(1, "npc"))))
                                                   
  # Function to create combined histogram plot with facets for min and max parameters
  create_facet_histogram <- function(data, min_var, max_var, parameter) {
  	# Calculate quantiles for each comparison
  	quantiles <- data %>%
    		group_by(comparison) %>%
    		summarise(min_ci = quantile(!!sym(min_var), probs = 0.025, na.rm = TRUE),
	      	max_ci = quantile(!!sym(max_var), probs = 0.975, na.rm = TRUE))
  
	p <- ggplot(data) +
	    	geom_histogram(aes_string(x = min_var, fill = "'Min'"), alpha = 0.5, binwidth = 0.05, position = "identity") +
	    	geom_histogram(aes_string(x = max_var, fill = "'Max'"), alpha = 0.5, binwidth = 0.05, position = "identity") +
	    	geom_vline(data = quantiles, aes(xintercept = min_ci), color = "purple", linetype = "dashed", size = 1) +
	    	geom_vline(data = quantiles, aes(xintercept = max_ci), color = "purple", linetype = "dashed", size = 1) +
	    	labs(title = paste("Histogram of Min and Max", parameter, "by Comparison"),
			 x = paste(parameter, "Value"), y = "Frequency") +
	    	scale_fill_manual(name = "Type", values = c("Min" = "blue", "Max" = "red"), labels = c("Min", "Max")) +
	    	theme_minimal() +
	    	facet_wrap(~comparison, scales = "free")
	  
	return(p)
  }

  # Example usage: Create the facet histogram plots for the parameters sigma, rho1, rho2, and x
  facet_histogram_sigma <- create_facet_histogram(all_inner_results, "minsigma", "maxsigma", "Sigma")
  facet_histogram_rho1 <- create_facet_histogram(all_inner_results, "min_rho1", "max_rho1", "rho1")
  facet_histogram_rho2 <- create_facet_histogram(all_inner_results, "min_rho2", "max_rho2", "rho2")
  facet_histogram_x <- create_facet_histogram(all_inner_results, "min_x", "max_x", "x")

  # Plot the facet histograms
  print(facet_histogram_sigma)
  print(facet_histogram_rho1)
  print(facet_histogram_rho2)
  print(facet_histogram_x)
  
  # Save plots to files if save_format is specified
  if (!is.null(save)) {
  	if (!is.null(filename)) {
    		ggsave(filename = paste0(filename, "_Combined_", plot_type, ".", save), plot = combined_plot, width = 15, height = 30, dpi = 300)
    		ggsave(filename = paste0(filename, "_Histograms_sigma.", save), plot = facet_histogram_sigma, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(filename, "_Histograms_rho1.", save), plot = facet_histogram_rho1, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(filename, "_Histograms_rho2.", save), plot = facet_histogram_rho2, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(filename, "_Histograms_x.", save), plot = facet_histogram_x, width = 15, height = 30, dpi = 300) 
    	} else {
    		ggsave(filename = paste0("Combined_sigmaparameters_", plot_type, ".", save), plot = combined_plot, width = 15, height = 30, dpi = 300)
    		ggsave(filename = paste0("Histograms_sigma.", save), plot = facet_histogram_sigma, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0("Histograms_rho1.", save), plot = facet_histogram_rho1, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0("Histograms_rho2.", save), plot = facet_histogram_rho2, width = 15, height = 30, dpi = 300)  
    		ggsave(filename = paste0("Histograms_x.", save), plot = facet_histogram_x, width = 15, height = 30, dpi = 300) 
    	}
  }
  
  if (!is.null(save) && individual == TRUE) {
  	if (!is.null(filename)) {
	    ggsave(filename = paste0(filename,"Min_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(filename,"Max_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(filename,"Min_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(filename,"Max_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(filename,"Min_Sigma_vs_x_", plot_type, ".", save), plot = p3_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(filename,"Max_Sigma_vs_x_", plot_type, ".", save), plot = p3_max, width = 10, height = 10, dpi = 300)
	} else {
	    ggsave(filename = paste0("Min_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0("Max_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0("Min_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0("Max_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0("Min_Sigma_vs_x_", plot_type, ".", save), plot = p3_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0("Max_Sigma_vs_x_", plot_type, ".", save), plot = p3_max, width = 10, height = 10, dpi = 300)
	}
  }
}
