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
      size = "-log10(p-value)"
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

#' Plot Sigma Values Against Parameters
#'
#' This function creates plots to visualize the impact of parameters \code{rho1}, \code{rho2}, and \code{x} on sigma values.
#' It identifies and highlights the parameter points that minimize and maximize sigma for each comparison.
#' Optionally, it can save the plots in \code{jpg}, \code{png}, or \code{svg} format with 300 dpi resolution.
#'
#' @param all_inner_results A data frame containing sigma values and corresponding parameters.
#' @param save A character string specifying the format to save the plots (\code{NULL}, \code{"jpg"}, \code{"png"}, \code{"svg"}). Defaults to \code{NULL}.
#' @param plot_type A character string specifying the type of plot ("point" or "line"). Defaults to "line". "line" fits a quadratic line of best fit to the data points. "point" plots the individual points. 
#' @param filename A character string specifying an alternate name of the file. Defaults to \code{NULL} which saves as parameter combination names.
#' @param individual A boolean type specifying whether you would like to save individual plots. Defaults to \code{FALSE}.
#' @param bg A character string indicating the background color of the plot. Options are "white" (default) or "transparent".
#' @param dir_path A character string indicating the directory to store the plot. Default is \code{"./plots/"} which creates the plots directory in the current directory.
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
sigmaplot <- function(all_inner_results, bg="white", filename = NULL, save = NULL, individual = FALSE, dir_path="./plots/", plot_type = "line") {
  
  # Convert relevant columns to numeric
  all_inner_results <- all_inner_results %>%
    mutate(across(c(s, minsigma_absolute_minimum_covariance, maxsigma_absolute_maximum_covariance, minsigma_correlation_relativetaxa1_scale, minsigma_correlation_relativetaxa2_scale, minsigma_scale_variance, maxsigma_correlation_relativetaxa1_scale, maxsigma_correlation_relativetaxa2_scale, maxsigma_scale_variance), as.numeric))
          
    # Ensure 'comparison' column exists
  if (!"comparison" %in% colnames(all_inner_results)) {
    all_inner_results$comparison <- "default"  # Add a default value if the column doesn't exist
  }
  all_inner_results$comparison <- as.factor(all_inner_results$comparison)
  
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
  p1_min <- create_plot("minsigma_correlation_relativetaxa1_scale", "minsigma_absolute_minimum_covariance", "comparison", "Min Absolute Sigma vs Relativetaxa1_scale_correlation")
  p1_max <- create_plot("maxsigma_correlation_relativetaxa1_scale", "maxsigma_absolute_maximum_covariance", "comparison", "Max Absolute Sigma vs Relativetaxa1_scale_correlation")
  
  # Min and Max Sigma vs rho2
  p2_min <- create_plot("minsigma_correlation_relativetaxa2_scale", "minsigma_absolute_minimum_covariance", "comparison", "Min Absolute Sigma vs Relativetaxa2_scale_correlation")
  p2_max <- create_plot("maxsigma_correlation_relativetaxa2_scale", "maxsigma_absolute_maximum_covariance", "comparison", "Max Absolute Sigma vs Relativetaxa2_scale_correlation")
  
  # Min and Max Sigma vs x
  p3_min <- create_plot("minsigma_scale_variance", "minsigma_absolute_minimum_covariance", "comparison", "Min Absolute Sigma vs Scale Variance")
  p3_max <- create_plot("maxsigma_scale_variance", "maxsigma_absolute_maximum_covariance", "comparison", "Max Absolute Sigma vs Scale Variance")
  
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
                          top = grid::textGrob("Parameter Plots", gp = grid::gpar(fontsize = 16, fontface = "bold")),
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
	    	custom_theme +
	    	facet_wrap(~comparison, scales = "free")
	  
	return(p)
  }

  # Example usage: Create the facet histogram plots for the parameters sigma, rho1, rho2, and x
  facet_histogram_sigma <- create_facet_histogram(all_inner_results, "minsigma_absolute_minimum_covariance", "maxsigma_absolute_maximum_covariance", "Absolute Covariance")
  facet_histogram_rho1 <- create_facet_histogram(all_inner_results, "minsigma_correlation_relativetaxa1_scale", "maxsigma_correlation_relativetaxa1_scale", "Relativetaxa1_scale_correlation")
  facet_histogram_rho2 <- create_facet_histogram(all_inner_results, "minsigma_correlation_relativetaxa2_scale", "maxsigma_correlation_relativetaxa2_scale", "Relativetaxa1_scale_correlation")
  facet_histogram_x <- create_facet_histogram(all_inner_results, "minsigma_scale_variance", "maxsigma_scale_variance", "Scale Variance")

  # Plot the facet histograms
  print(facet_histogram_sigma)
  print(facet_histogram_rho1)
  print(facet_histogram_rho2)
  print(facet_histogram_x)
  
  # Save plots to files if save_format is specified
  if (!is.null(save)) {
  	if (!is.null(filename)) {
    		ggsave(filename = paste0(dir_path, filename, "_Combined_", plot_type, ".", save), plot = combined_plot, width = 15, height = 30, dpi = 300)
    		ggsave(filename = paste0(dir_path, filename, "_Histograms_abolutesigma.", save), plot = facet_histogram_sigma, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(dir_path, filename, "_Histograms_rho1.", save), plot = facet_histogram_rho1, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(dir_path, filename, "_Histograms_rho2.", save), plot = facet_histogram_rho2, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(dir_path, filename, "_Histograms_scalevariance.", save), plot = facet_histogram_x, width = 15, height = 30, dpi = 300) 
    	} else {
    		ggsave(filename = paste0(dir_path, "Combined_sigmaparameters_", plot_type, ".", save), plot = combined_plot, width = 15, height = 30, dpi = 300)
    		ggsave(filename = paste0(dir_path, "Histograms_absolutesigma.", save), plot = facet_histogram_sigma, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(dir_path, "Histograms_rho1.", save), plot = facet_histogram_rho1, width = 15, height = 30, dpi = 300) 
    		ggsave(filename = paste0(dir_path, "Histograms_rho2.", save), plot = facet_histogram_rho2, width = 15, height = 30, dpi = 300)  
    		ggsave(filename = paste0(dir_path, "Histograms_scalevariance.", save), plot = facet_histogram_x, width = 15, height = 30, dpi = 300) 
    	}
  }
  
  if (!is.null(save) && individual == TRUE) {
  	if (!is.null(filename)) {
	    ggsave(filename = paste0(dir_path, filename,"Min_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, filename,"Max_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, filename,"Min_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, filename,"Max_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, filename,"Min_Sigma_vs_scalevariance_", plot_type, ".", save), plot = p3_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, filename,"Max_Sigma_vs_scalevariance_", plot_type, ".", save), plot = p3_max, width = 10, height = 10, dpi = 300)
	} else {
	    ggsave(filename = paste0(dir_path, "Min_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, "Max_Sigma_vs_rho1_", plot_type, ".", save), plot = p1_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, "Min_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, "Max_Sigma_vs_rho2_", plot_type, ".", save), plot = p2_max, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, "Min_Sigma_vs_scalevariance_", plot_type, ".", save), plot = p3_min, width = 10, height = 10, dpi = 300)
	    ggsave(filename = paste0(dir_path, "Max_Sigma_vs_scalevariance_", plot_type, ".", save), plot = p3_max, width = 10, height = 10, dpi = 300)
	}
  }
}
