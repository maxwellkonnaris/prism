#' Forest Plot of Confidence Intervals for Multiple Data Frames
#'
#' This function creates a professional-looking forest plot of confidence intervals for multiple data frames. Each data frame is plotted side by side with shared y-axis labels, and each plot is labeled with the name of the data frame.
#'
#' @param data_list A named list of data frames. Each data frame should contain the necessary columns as specified below.
#' @param bg A character string indicating the background color of the plot. Options are "white" (default) or "transparent".
#' @param save A character string indicating the file format to save the plot. Options are "png", "jpg", "svg", "pdf". Default is NULL, which means the plot is not saved.
#' @param filename A character string indicating the file name when saving the plot. Default is NULL, which means the plot is saved as forest_plot if save format is indicated.
#' @param dir_path A character string indicating the directory to store the plot. Default is \code{"./plots/"} which creates the plots directory in the current directory.
#' @return A ggplot object representing the combined forest plot.
#' @import ggplot2
#' @import dplyr
#' @export
forest_plot <- function(data_list, bg = "white", save = NULL, filename = NULL, dir_path = "./plots/") {

  # Check if input is a data frame, convert to a named list if true
  if (is.data.frame(data_list)) {
    data_list <- list(PRISM = data_list)
  }
  
  # Check that data_list is a named list of data frames
  if (!is.list(data_list) || is.null(names(data_list))) {
    stop("data_list must be a named list of data frames.")
  }
  
  # Required columns
  required_columns <- c("comparison", "ninetyfive_ci_lower", "ninetyfive_ci_upper",
                        "minsigma_absolute_minimum_covariance", "maxsigma_absolute_maximum_covariance")
  
  # Initialize an empty list to store data frames with an added 'Dataset' column
  data_frames <- list()
  
  # First dataset (for reordering comparison levels)
  first_dataset_name <- names(data_list)[1]
  first_data <- data_list[[first_dataset_name]]
  
  # Ensure the first data has the necessary columns
  if (!all(required_columns %in% colnames(first_data))) {
    stop(paste("Data frame '", first_dataset_name, "' must contain columns:", paste(required_columns, collapse = ", ")))
  }
  
  # Reorder the comparison names in the first dataset by the lower bound of the 95% confidence interval
  comparison_order <- first_data %>%
    arrange(ninetyfive_ci_lower) %>%
    pull(comparison)
  
  # Loop over each data frame in data_list
  for (dataset_name in names(data_list)) {
    data <- data_list[[dataset_name]]
    
    # Ensure the data has the necessary columns
    if (!all(required_columns %in% colnames(data))) {
      stop(paste("Data frame '", dataset_name, "' must contain columns:", paste(required_columns, collapse = ", ")))
    }
    
    # Check if p_value is provided and add a flag column
    p_value_provided <- "p_value" %in% colnames(data)
    data$p_value_provided <- p_value_provided  # Add a flag column
    
    # Add a column for the dataset name
    data$Dataset <- dataset_name
    
    # Reorder the comparison factor based on the first dataset
    data$comparison <- factor(data$comparison, levels = comparison_order)
    
    # Append to the list
    data_frames[[dataset_name]] <- data
  }
  
  # Combine all data frames into one
  combined_data <- bind_rows(data_frames)
  
  # Highlight intervals that do not cover 0
  combined_data <- combined_data %>%
    mutate(highlight = ifelse(
      (ninetyfive_ci_lower > 0 & ninetyfive_ci_upper > 0) | (ninetyfive_ci_lower < 0 & ninetyfive_ci_upper < 0),
      "95% CI Doesn't Cover Zero", "95% CI Covers Zero"
    ))
  
  # Create the forest plot
  plot <- ggplot(combined_data, aes(x = comparison)) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Taxa Comparison",
      y = "Estimated Covariance",
      color = NULL
    ) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    ) +
    scale_color_manual(
      values = c(
        "95% CI Doesn't Cover Zero" = "#023E8A",  # Blue
        "95% CI Covers Zero" = "#BEBEBE",
        "Covariance Range" = "#676767"
      ),
      breaks = c("95% CI Doesn't Cover Zero", "Covariance Range") # Exclude "Covers Zero" from legend
    )
  
  # Add the covariance range error bars
  plot <- plot + geom_errorbar(
    aes(ymin = minsigma_absolute_minimum_covariance, ymax = maxsigma_absolute_maximum_covariance, color = "Covariance Range"),
    width = 0.4, size = 0.8
  )
  
  # Add the 95% confidence intervals
  plot <- plot + geom_errorbar(
    aes(ymin = ninetyfive_ci_lower, ymax = ninetyfive_ci_upper, color = highlight),
    width = 0.6, size = 1
  )
  
  # Add points for the p-values if provided
  if (any(combined_data$p_value_provided)) {
    plot <- plot + geom_point(
      data = combined_data %>% filter(p_value_provided == TRUE),
      aes(y = (ninetyfive_ci_lower + ninetyfive_ci_upper) / 2, size = -log10(p_value)),
      color = "black", shape = 21, fill = "white", stroke = 1
    ) +
      labs(size = expression("-log"[10]*"(p-value)")) +
      scale_size_continuous(range = c(2, 6))
  }
  
  # Customize the background based on the bg parameter
  if (bg == "transparent") {
    plot <- plot + 
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
  } else if (bg == "white") {
    plot <- plot + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
  } else {
    stop("bg parameter must be 'transparent' or 'white'")
  }
  
  # Add subtle gridlines
  plot <- plot + 
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray80", linetype = "dashed")
    )
  
  # Adjust y-axis breaks
  plot <- plot + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
  
  # Facet by Dataset to create side-by-side plots
  plot <- plot + facet_grid(. ~ Dataset, scales = "free_x", space = "free_x")
  
  # Adjust the plot size based on the number of comparisons
  scaling_factor <- 0.5
  dynamic_height <- length(unique(combined_data$comparison)) * scaling_factor
  min_height <- 10
  max_height <- 35
  final_height <- max(min_height, min(dynamic_height, max_height))
  
  # Save the plot if save is not NULL
  if (!is.null(save)) {
    # Check if directory exists
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    # Ensure the directory path ends with a slash
    if (!grepl("/$", dir_path)) {
      dir_path <- paste0(dir_path, "/")
    }
    # Set default filename if not provided
    if (is.null(filename)) {
      filename <- "forest_plot"
    }
    file_name <- paste0(dir_path, filename, ".", save)
    
    # Adjust width based on the number of data frames
    plot_width <- 8 * length(data_list)
    
    ggsave(file_name, plot, width = plot_width, height = final_height, dpi = 300, device = save, bg = bg, limitsize = FALSE)
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
#' @import GGally
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

    # Create bivariate plot grid
  bivariate_plot <- ggpairs(all_inner_results, columns = c("minsigma_correlation_relativetaxa1_scale", 
                                                           "maxsigma_correlation_relativetaxa1_scale", 
                                                           "minsigma_correlation_relativetaxa2_scale", 
                                                           "maxsigma_correlation_relativetaxa2_scale", 
                                                           "minsigma_scale_variance", 
                                                           "maxsigma_scale_variance", 
                                                           "minsigma_absolute_minimum_covariance", 
                                                           "maxsigma_absolute_maximum_covariance"), 
                            aes(color = comparison), 
                            title = "Bivariate Plot")
  
  # Print the bivariate plot grid
  print(bivariate_plot)
  
  # Save the plot if needed
  if (!is.null(filename)) {
    ggsave(filename = file.path(dir_path, filename), plot = bivariate_plot, width = 15, height = 15, units = "in", bg = bg)
  }
  
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

#' Diagnostics for Bootstrap Convergence Results
#'
#' This function generates diagnostics plots and statistics for bootstrap convergence results from the `estimate_covariance_convergence` function. It includes convergence plots, standard error and confidence intervals, effective sample size, Rhat (Gelman-Rubin Diagnostic), cumulative mean and variance plots, resampling diagnostics (Jackknife-after-Bootstrap, Bootstrap-after-Bootstrap), Monte Carlo standard error, and comparison of subsamples.
#'
#' @param combined_results Dataframe combining the results from all the incremental bootstrap samples. This dataframe should include columns such as \code{S} (number of bootstrap samples), \code{minsigma_absolute_minimum_covariance} (minimum covariance estimate), \code{maxsigma_absolute_maximum_covariance} (maximum covariance estimate), \code{ninetyfive_ci_lower} (lower bound of the 95% confidence interval), \code{ninetyfive_ci_upper} (upper bound of the 95% confidence interval), and \code{comparison} (taxa comparison identifier).
#' @param convergence_results List of results for each incrementally increased number of bootstrap samples.
#' @return A list containing various diagnostics plots and statistics:
#' \item{convergence_plot}{A plot showing the convergence of bootstrap estimates as a function of the number of bootstrap samples.}
#' \item{ci_plot}{A plot displaying the standard error and 95\% confidence intervals of the bootstrap estimates.}
#' \item{ess_values}{Effective sample size (ESS) for each combination of bootstrap sample size and taxa comparison.}
#' \item{rhat_values}{Rhat (Gelman-Rubin Diagnostic) values for each combination of bootstrap sample size and taxa comparison.}
#' \item{cumulative_mean_variance_plot}{Plots showing the cumulative mean and variance of bootstrap estimates as a function of the number of bootstrap samples.}
#' \item{jackknife_results}{Jackknife estimates of the bootstrap results.}
#' \item{bootstrap_resamples}{Bootstrap-after-Bootstrap estimates of the bootstrap results.}
#' \item{mcse_values}{Monte Carlo Standard Error (MCSE) values for each combination of bootstrap sample size and taxa comparison.}
#' \item{subsample_comparison_plot}{A plot comparing subsample estimates across different groups and bootstrap sample sizes.}
#' @details
#' This function provides a comprehensive set of diagnostics to assess the stability and convergence of bootstrap estimates obtained from the `estimate_covariance_convergence` function. The diagnostics include:
#' \enumerate{
#'   \item Convergence Plots: Visualize how the bootstrap estimates stabilize as the number of bootstrap samples increases.
#'   \item Standard Error and Confidence Intervals: Display the standard error and 95\% confidence intervals for the bootstrap estimates.
#'   \item Effective Sample Size (ESS): Calculate the ESS to determine the number of effectively independent samples in the bootstrap estimates.
#'   \item Rhat (Gelman-Rubin Diagnostic): Assess the convergence of the bootstrap chains using the Rhat statistic.
#'   \item Cumulative Mean and Variance Plots: Show how the cumulative mean and variance of the bootstrap estimates evolve with increasing bootstrap samples.
#'   \item Resampling Diagnostics: Perform Jackknife-after-Bootstrap and Bootstrap-after-Bootstrap resampling to further assess the stability of the bootstrap estimates.
#'   \item Monte Carlo Standard Error (MCSE): Calculate the MCSE to quantify the variability in the bootstrap estimates.
#'   \item Comparing Subsamples: Compare the bootstrap estimates across different subsample groups to assess consistency.
#' }
#' @examples
#' # Example usage:
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 10)
#' results <- estimate_covariance_convergence(Y)
#' diagnostics <- diagnose_bootstrap_convergence(results$combined_results, results$convergence_results)
#' @import ggplot2
#' @import coda
#' @import dplyr
#' @import gridExtra
#' @export
diagnose_bootstrap_convergence <- function(combined_results, convergence_results) {
  
  # Custom cumulative variance function
  cumvar <- function(x) {
    n <- length(x)
    cumsum((x - cumsum(x) / seq_along(x))^2) / seq_along(x)
  }
  
  # Plot convergence of bootstrap estimates
  plot_convergence <- function(bootstrap_results) {
    ggplot(bootstrap_results, aes(x = S, y = (minsigma_absolute_minimum_covariance + maxsigma_absolute_maximum_covariance) / 2, color = comparison)) +
      geom_line() +
      labs(title = "Convergence Plot of Bootstrap Estimates",
           x = "Bootstrap Sample Size",
           y = "Mean Estimate")
  }
  
  # Plot standard error and confidence intervals
  plot_standard_error_and_ci <- function(bootstrap_results) {
    ggplot(bootstrap_results, aes(x = S, y = (minsigma_absolute_minimum_covariance + maxsigma_absolute_maximum_covariance) / 2, color = comparison)) +
      geom_point() +
      geom_errorbar(aes(ymin = ninetyfive_ci_lower, ymax = ninetyfive_ci_upper), width = 0.2) +
      labs(title = "Bootstrap Estimates with 95% CI",
           x = "Bootstrap Sample Size",
           y = "Estimate")
  }
  
  # Calculate effective sample size
  calculate_effective_sample_size <- function(bootstrap_estimates) {
    ess <- function(x) {
      n <- length(x)
      acf_x <- acf(x, plot = FALSE)
      return(n / (1 + 2 * sum(acf_x$acf[-1])))
    }
    bootstrap_estimates %>%
      group_by(S, comparison) %>%
      summarise(
        ess_min = ess(minsigma_absolute_minimum_covariance),
        ess_max = ess(maxsigma_absolute_maximum_covariance)
      )
  }
  
  # Calculate Rhat (Gelman-Rubin Diagnostic)
  calculate_rhat <- function(bootstrap_estimates) {
    chains_min <- split(bootstrap_estimates$minsigma_absolute_minimum_covariance, bootstrap_estimates$S)
    chains_max <- split(bootstrap_estimates$maxsigma_absolute_maximum_covariance, bootstrap_estimates$S)
    mcmc_chains_min <- coda::mcmc.list(lapply(chains_min, mcmc))
    mcmc_chains_max <- coda::mcmc.list(lapply(chains_max, mcmc))
    list(
      rhat_min = coda::gelman.diag(mcmc_chains_min)$psrf,
      rhat_max = coda::gelman.diag(mcmc_chains_max)$psrf
    )
  }
  
  # Plot cumulative mean and variance
  plot_cumulative_mean_variance <- function(bootstrap_estimates) {
    cumulative_stats <- bootstrap_estimates %>%
      group_by(S, comparison) %>%
      arrange(S) %>%
      mutate(
        cumulative_mean_min = cummean(minsigma_absolute_minimum_covariance),
        cumulative_mean_max = cummean(maxsigma_absolute_maximum_covariance),
        cumulative_variance_min = cumvar(minsigma_absolute_minimum_covariance),
        cumulative_variance_max = cumvar(maxsigma_absolute_maximum_covariance)
      )
    
    mean_plot_min <- ggplot(cumulative_stats, aes(x = S, y = cumulative_mean_min, color = comparison)) +
      geom_line() +
      labs(title = "Cumulative Mean (Min)",
           x = "Bootstrap Sample Size",
           y = "Cumulative Mean")
    
    mean_plot_max <- ggplot(cumulative_stats, aes(x = S, y = cumulative_mean_max, color = comparison)) +
      geom_line() +
      labs(title = "Cumulative Mean (Max)",
           x = "Bootstrap Sample Size",
           y = "Cumulative Mean")
    
    variance_plot_min <- ggplot(cumulative_stats, aes(x = S, y = cumulative_variance_min, color = comparison)) +
      geom_line() +
      labs(title = "Cumulative Variance (Min)",
           x = "Bootstrap Sample Size",
           y = "Cumulative Variance")
    
    variance_plot_max <- ggplot(cumulative_stats, aes(x = S, y = cumulative_variance_max, color = comparison)) +
      geom_line() +
      labs(title = "Cumulative Variance (Max)",
           x = "Bootstrap Sample Size",
           y = "Cumulative Variance")
    
    gridExtra::grid.arrange(mean_plot_min, mean_plot_max, variance_plot_min, variance_plot_max, ncol = 2)
  }
  
  # Resampling diagnostics (Jackknife-after-Bootstrap, Bootstrap-after-Bootstrap)
  jackknife_after_bootstrap <- function(bootstrap_estimates) {
    n <- length(bootstrap_estimates$minsigma_absolute_minimum_covariance)
    list(
      jackknife_min = sapply(1:n, function(i) {
        mean(bootstrap_estimates$minsigma_absolute_minimum_covariance[-i])
      }),
      jackknife_max = sapply(1:n, function(i) {
        mean(bootstrap_estimates$maxsigma_absolute_maximum_covariance[-i])
      })
    )
  }
  
  bootstrap_after_bootstrap <- function(bootstrap_estimates, num_resamples = 1000) {
    list(
      bootstrap_min = replicate(num_resamples, {
        resample_indices <- sample(seq_along(bootstrap_estimates$minsigma_absolute_minimum_covariance), replace = TRUE)
        mean(bootstrap_estimates$minsigma_absolute_minimum_covariance[resample_indices])
      }),
      bootstrap_max = replicate(num_resamples, {
        resample_indices <- sample(seq_along(bootstrap_estimates$maxsigma_absolute_maximum_covariance), replace = TRUE)
        mean(bootstrap_estimates$maxsigma_absolute_maximum_covariance[resample_indices])
      })
    )
  }
  
  # Calculate Monte Carlo Standard Error
  calculate_mcse <- function(bootstrap_estimates) {
    mcse <- function(x) {
      sd(x) / sqrt(length(x))
    }
    bootstrap_estimates %>%
      group_by(S, comparison) %>%
      summarise(
        mcse_min = mcse(minsigma_absolute_minimum_covariance),
        mcse_max = mcse(maxsigma_absolute_maximum_covariance)
      )
  }
  
  # Compare subsamples
  compare_subsamples <- function(bootstrap_estimates, num_groups = 2) {
    subsample_comparison <- bootstrap_estimates %>%
      group_by(S, comparison) %>%
      mutate(group = ntile(row_number(), num_groups)) %>%
      group_by(group, add = TRUE) %>%
      summarise(
        mean_min = mean(minsigma_absolute_minimum_covariance),
        mean_max = mean(maxsigma_absolute_maximum_covariance)
      )
    
    ggplot(subsample_comparison, aes(x = group, y = (mean_min + mean_max) / 2, color = factor(S))) +
      geom_line() +
      geom_point() +
      labs(title = "Comparison of Subsamples",
           x = "Subsample Group",
           y = "Mean Estimate",
           color = "Sample Size")
  }
  
  # Generate all diagnostics
  convergence_plot <- plot_convergence(combined_results)
  ci_plot <- plot_standard_error_and_ci(combined_results)
  ess_values <- calculate_effective_sample_size(combined_results)
  rhat_values <- calculate_rhat(combined_results)
  cumulative_mean_variance_plot <- plot_cumulative_mean_variance(combined_results)
  jackknife_results <- jackknife_after_bootstrap(combined_results)
  bootstrap_resamples <- bootstrap_after_bootstrap(combined_results)
  mcse_values <- calculate_mcse(combined_results)
  subsample_comparison_plot <- compare_subsamples(combined_results)
  
  return(list(
    convergence_plot = convergence_plot,
    ci_plot = ci_plot,
    ess_values = ess_values,
    rhat_values = rhat_values,
    cumulative_mean_variance_plot = cumulative_mean_variance_plot,
    jackknife_results = jackknife_results,
    bootstrap_resamples = bootstrap_resamples,
    mcse_values = mcse_values,
    subsample_comparison_plot = subsample_comparison_plot
  ))
}


#' Plot and Save 3D Scatter Plots with Convex Hull for All Comparisons
#'
#' This function generates 3D scatter plots for `minsigma`, `maxsigma`,
#' and optionally highlights `SPSD < 0` points for all comparisons in the dataset.
#' Convex hulls are drawn around the points, and the user can choose to plot
#' all comparisons on the same plot or create individual plots.
#' The plots are saved as HTML files.
#'
#' @param data A data frame containing the optimization results. 
#' The data frame must contain columns named "comparison", 
#' "minsigma_correlation_relativetaxa1_scale", 
#' "minsigma_correlation_relativetaxa2_scale", 
#' "minsigma_scale_sd", 
#' "maxsigma_correlation_relativetaxa1_scale", 
#' "maxsigma_correlation_relativetaxa2_scale", and 
#' "maxsigma_scale_sd". Additionally, if the data frame contains an 
#' "SPSD" column, the plot will reflect this using color.
#' @param output_directory A string specifying the directory where the plots will be saved.
#' @param correlation_relativetaxa_scale_range A vector specifying the range of values for the x and y axes (default: c(-1, 1)).
#' @param scale_sd_range A vector specifying the range of values for the z-axis (default: c(0.49, 0.51)).
#' @param plot_all Logical. If TRUE, plot all comparisons on the same 3D scatter plot. 
#' If FALSE, generate separate plots for each comparison.
#' @param alpha_value A numeric value that controls the convex hull's tightness (default: 1).
#' @return This function saves the plots and returns no value.
#' @import plotly htmlwidgets alphashape3d
#' @export
plot_comparisons_3d_scatter <- function(data, output_directory = "plots/", 
                                        correlation_relativetaxa_scale_range = c(-1, 1), 
                                        scale_sd_range = c(0.49, 0.51),
                                        plot_all = TRUE,
                                        alpha_value = 1) {
  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Get unique comparisons from the data
  unique_comparisons <- unique(data$comparison)
  
  # Create a custom color palette excluding grey and red shades
  all_colors <- grDevices::colors()[grep('gr(a|e)y|red', grDevices::colors(), invert = TRUE)]
  
  # Sample N colors for each comparison
  color_palette <- sample(all_colors, min(length(unique_comparisons), length(all_colors)))

  # Initialize an empty plot if plotting all together
  if (plot_all) {
    plot <- plot_ly()
  }

  # Loop through each comparison
  for (i in seq_along(unique_comparisons)) {
    comparison_value <- unique_comparisons[i]
    # Filter the data for the specific comparison
    data_subset <- subset(data, comparison == comparison_value)
    
    # Marker shapes: circle for minsigma, square for maxsigma
    minsigma_shape <- 'circle'
    maxsigma_shape <- 'square'
    
    # Check if SPSD column exists and highlight points where SPSD < 0
    if ("SPSD" %in% names(data_subset)) {
      minsigma_colors <- ifelse(data_subset$SPSD < 0, "red", color_palette[i])
      maxsigma_colors <- ifelse(data_subset$SPSD < 0, "red", color_palette[i])
    } else {
      # Use comparison colors if SPSD is not available
      minsigma_colors <- color_palette[i]
      maxsigma_colors <- color_palette[i]
    }
    
    # Get the points for both minsigma and maxsigma
    minsigma_points <- data.frame(
      x = data_subset$minsigma_correlation_relativetaxa1_scale,
      y = data_subset$minsigma_correlation_relativetaxa2_scale,
      z = data_subset$minsigma_scale_sd
    )
    
    maxsigma_points <- data.frame(
      x = data_subset$maxsigma_correlation_relativetaxa1_scale,
      y = data_subset$maxsigma_correlation_relativetaxa2_scale,
      z = data_subset$maxsigma_scale_sd
    )
    
    # Combine minsigma and maxsigma points for creating the convex hull
    all_points <- rbind(minsigma_points, maxsigma_points)
    
    # Remove duplicate points
    all_points <- unique(all_points)
    
    # If there are fewer than 4 unique points, we can't create a 3D convex hull
    if (nrow(all_points) >= 4) {
      # Center the points by subtracting the mean (shifting the center of the points to the origin)
      all_points_centered <- scale(all_points, center = TRUE, scale = FALSE)
      
      # Try to create the convex hull using alphashape3d (alpha shape is like a convex hull)
      ashape <- tryCatch({
        ashape3d(all_points_centered, alpha = alpha_value)  # Adjust alpha for tighter or looser fit
      }, error = function(e) {
        message("Failed to create convex hull for comparison: ", comparison_value)
        return(NULL)
      })
      
      if (!is.null(ashape)) {
        # Extract the vertices and faces for the convex hull mesh
        vertices <- ashape$alpha3d$triang
        faces <- ashape$alpha3d$facets
      }
    }
    
    if (plot_all) {
      # Use the same color for minsigma and maxsigma for the same comparison
      minsigma_colors <- color_palette[i]
      maxsigma_colors <- color_palette[i]
      
      # Add minsigma scatter plot with 'circle' markers to the overall plot
      plot <- plot %>%
        add_markers(data = data_subset, 
                    x = ~minsigma_correlation_relativetaxa1_scale, 
                    y = ~minsigma_correlation_relativetaxa2_scale, 
                    z = ~minsigma_scale_sd, 
                    marker = list(symbol = minsigma_shape, color = minsigma_colors, size = 4),
                    name = paste('MinSigma', comparison_value))
      
      # Add maxsigma scatter plot with 'square' markers to the overall plot
      plot <- plot %>%
        add_markers(data = data_subset, 
                    x = ~maxsigma_correlation_relativetaxa1_scale, 
                    y = ~maxsigma_correlation_relativetaxa2_scale, 
                    z = ~maxsigma_scale_sd, 
                    marker = list(symbol = maxsigma_shape, color = maxsigma_colors, size = 4),
                    name = paste('MaxSigma', comparison_value))
      
      # Add convex hull if it was successfully created
      if (exists("ashape") && !is.null(ashape)) {
        plot <- plot %>%
          add_trace(type = 'mesh3d',
                    x = vertices[,1], y = vertices[,2], z = vertices[,3],
                    i = faces[,1] - 1, j = faces[,2] - 1, k = faces[,3] - 1, 
                    color = color_palette[i], opacity = 0.3, name = paste("Convex Hull", comparison_value))
      }
    } else {
      # For individual plots, use different colors for minsigma and maxsigma
      minsigma_colors <- sample(all_colors, 1)
      maxsigma_colors <- sample(all_colors, 1)
      
      # Generate individual plots for each comparison
      comparison_plot <- plot_ly() %>%
        
        # Add minsigma scatter plot with 'circle' markers
        add_markers(data = data_subset, 
                    x = ~minsigma_correlation_relativetaxa1_scale, 
                    y = ~minsigma_correlation_relativetaxa2_scale, 
                    z = ~minsigma_scale_sd, 
                    marker = list(symbol = minsigma_shape, color = minsigma_colors, size = 4),
                    name = paste('MinSigma', comparison_value)) %>%
        
        # Add maxsigma scatter plot with 'square' markers
        add_markers(data = data_subset, 
                    x = ~maxsigma_correlation_relativetaxa1_scale, 
                    y = ~maxsigma_correlation_relativetaxa2_scale, 
                    z = ~maxsigma_scale_sd, 
                    marker = list(symbol = maxsigma_shape, color = maxsigma_colors, size = 4),
                    name = paste('MaxSigma', comparison_value)) 
        
      # Add convex hull if it was successfully created
      if (exists("ashape") && !is.null(ashape)) {
        comparison_plot <- comparison_plot %>%
          add_trace(type = 'mesh3d',
                    x = vertices[,1], y = vertices[,2], z = vertices[,3],
                    i = faces[,1] - 1, j = faces[,2] - 1, k = faces[,3] - 1, 
                    color = minsigma_colors, opacity = 0.3, name = paste("Convex Hull", comparison_value))
      }
      
      # Set layout with custom axis ranges
      comparison_plot <- comparison_plot %>%
        layout(scene = list(xaxis = list(title = 'Correlation RelTaxa1 Scale', range = correlation_relativetaxa_scale_range),
                            yaxis = list(title = 'Correlation RelTaxa2 Scale', range = correlation_relativetaxa_scale_range),
                            zaxis = list(title = 'Scale SD', range = scale_sd_range)),
               title = paste("MinSigma and MaxSigma Scatter Plot for Comparison", comparison_value))
      
      # Define file path for saving the plot
      output_file <- file.path(output_directory, paste0("sigma_scatter_comparison_", comparison_value, ".html"))
      
      # Save the individual plot as an HTML file
      htmlwidgets::saveWidget(comparison_plot, file = output_file)
      
      # Optionally, print a message to confirm plot generation
      message("Generated and saved individual scatter plot for comparison: ", comparison_value)
    }
  }

  # If plotting all comparisons together, finalize and save the plot
  if (plot_all) {
    # Set layout with custom axis ranges
    plot <- plot %>%
      layout(scene = list(xaxis = list(title = 'Correlation RelTaxa1 Scale', range = correlation_relativetaxa_scale_range),
                          yaxis = list(title = 'Correlation RelTaxa2 Scale', range = correlation_relativetaxa_scale_range),
                          zaxis = list(title = 'Scale SD', range = scale_sd_range)),
             title = "MinSigma and MaxSigma Scatter Plot for All Comparisons")
    
    # Define file path for saving the combined plot
    output_file <- file.path(output_directory, "sigma_scatter_all_comparisons.html")
    
    # Save the combined plot as an HTML file
    htmlwidgets::saveWidget(plot, file = output_file)
    
    # Optionally, print a message to confirm plot generation
    message("Generated and saved combined scatter plot for all comparisons.")
  }
}

#' Plot and Save 3D Scatter Plots for rpars
#'
#' This function generates 3D scatter plots for `sigma`, `rho1`, `rho2`, and `scalestdevstep`
#' and optionally highlights `SPSD < 0` points for all comparisons in the dataset.
#' Convex hulls are drawn around the points, and the user can choose to plot
#' all comparisons on the same plot or create individual plots.
#' The plots are saved as HTML files.
#'
#' @param data A data frame containing the optimization results, including `sigma`, `rho1`, `rho2`, and `scalestdevstep`.
#' @param output_directory A string specifying the directory where the plots will be saved.
#' @param rho_scale_range A vector specifying the range of values for `rho1` and `rho2` axes (default: c(-1, 1)).
#' @param scalestdevstep_range A vector specifying the range of values for the `scalestdevstep` axis (default: c(0.49, 0.51)).
#' @param plot_all Logical. If TRUE, plot all comparisons on the same 3D scatter plot. 
#' If FALSE, generate separate plots for each comparison.
#' @param alpha_value A numeric value that controls the convex hull's tightness (default: 1).
#' @return This function saves the plots and returns no value.
#' @import plotly htmlwidgets alphashape3d
#' @export
plot_rpars_3d_scatter <- function(data, output_directory = "plots/", rho_scale_range = c(min(data$rho1, na.rm = TRUE), max(data$rho1, na.rm = TRUE)), scalestdevstep_range = c(min(data$scalestdevstep, na.rm = TRUE), max(data$scalestdevstep, na.rm = TRUE)), plot_all = TRUE, alpha_value = 1) {
  
  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Step 1: Group by 'comparison' and then group by the shared columns (rho1, rho2, scalestdevstep) within each comparison
  averaged_data <- data %>%
    group_by(comparison, rho1, rho2, scalestdevstep) %>%
    summarise(across(where(is.numeric), mean, .names = "mean_{.col}"), .groups = 'drop')

  # Get unique comparisons from the data
  unique_comparisons <- unique(averaged_data$comparison)
  
  # Initialize an empty plot if plotting all together
  if (plot_all) {
    plot <- plot_ly(type = 'scatter3d', mode = 'markers')
  }

  # Loop through each comparison
  for (i in seq_along(unique_comparisons)) {
    comparison_value <- unique_comparisons[i]
    
    # Filter the data for the specific comparison
    data_subset <- subset(averaged_data, comparison == comparison_value)

    # Separate the data into two parts: SPSD >= 0 and SPSD < 0
    data_spsd_pos <- subset(data_subset, mean_SPSD >= 0)
    data_spsd_neg <- subset(data_subset, mean_SPSD < 0)

    # Determine the size scale for sigma values
    size_min <- 1  # minimum marker size (same as SPSD < 0 markers)
    size_max <- 8  # maximum marker size (2x the size_min)
    sigma_range <- range(data_spsd_pos$mean_sigma, na.rm = TRUE)  # get the range of sigma values
    
    if (diff(sigma_range) == 0) {
      sigma_range[2] <- sigma_range[1] + 1  # prevent division by zero
    }

    # Map sigma to marker size for SPSD >= 0 points
    marker_size_spsd_pos <- size_min + (data_spsd_pos$mean_sigma - sigma_range[1]) / (sigma_range[2] - sigma_range[1]) * (size_max - size_min)

    # Plot points where SPSD < 0 (in red) with constant size
    if (nrow(data_spsd_neg) > 0) {
      plot <- plot %>%
        add_markers(data = data_spsd_neg, 
                    x = ~rho1, 
                    y = ~rho2, 
                    z = ~scalestdevstep, 
                    marker = list(symbol = 'circle', color = "red", size = size_min),
                    name = "SPSD < 0 (Red)",
                    hoverinfo = "text",
                    text = ~paste("rho1:", rho1, "<br>rho2:", rho2, "<br>scalestdevstep:", scalestdevstep, "<br>sigma:", mean_sigma, "<br>SPSD:", mean_SPSD))
    }

    # Plot points where SPSD >= 0 with a gradient color and varying size based on sigma values
    if (nrow(data_spsd_pos) > 0) {
      plot <- plot %>%
        add_markers(data = data_spsd_pos, 
                    x = ~rho1, 
                    y = ~rho2, 
                    z = ~scalestdevstep, 
                    marker = list(symbol = 'circle', 
                                  color = ~mean_sigma, 
                                  colorscale = 'Viridis', 
                                  colorbar = list(title = "Sigma", len = 0.4),
                                  size = marker_size_spsd_pos,
                                  sizemode = 'diameter',
                                  sizeref = 2 * max(marker_size_spsd_pos, na.rm = TRUE) / (size_max ^ 2)),
                    name = "SPSD >= 0 (Sigma Gradient)",
                    hoverinfo = "text",
                    text = ~paste("rho1:", rho1, "<br>rho2:", rho2, "<br>scalestdevstep:", scalestdevstep, "<br>sigma:", mean_sigma, "<br>SPSD:", mean_SPSD))
    }
  }

  # Finalize and add legend and color scale
  if (plot_all) {
    plot <- plot %>%
      layout(
        scene = list(
          xaxis = list(title = 'rho1', range = rho_scale_range),
          yaxis = list(title = 'rho2', range = rho_scale_range),
          zaxis = list(title = 'scalestdevstep', range = scalestdevstep_range)
        ),
        title = "Sigma Scatter Plot for All Comparisons",
        showlegend = TRUE
      )
    
    # Define file path for saving the combined plot
    output_file <- file.path(output_directory, "sigma_scatter_all_comparisons.html")
    
    # Save the combined plot as an HTML file
    htmlwidgets::saveWidget(plot, file = output_file)
    
    # Optionally, print a message to confirm plot generation
    message("Generated and saved combined scatter plot for all comparisons.")
  }
}


#' 3D Scatter Plot for Sigma and SPSD Values with Confidence Intervals
#'
#' This function generates a 3D scatter plot to visualize the relationship 
#' between `rho1`, `rho2`, `scalestdevstep`, and `sigma` values. It also 
#' highlights points where the SPSD condition is uncertain, based on confidence 
#' intervals (CI). Points where the SPSD condition is uncertain are displayed 
#' in red, while those where SPSD is certain are shown with a color gradient based 
#' on `sigma` values.
#'
#' @param data A data frame containing the variables `rho1`, `rho2`, `scalestdevstep`, 
#'        `SPSD`, `sigma`, and `comparison` for grouping. This dataset should 
#'        contain multiple rows with repeated `rho1`, `rho2`, `scalestdevstep`, 
#'        and `comparison` values to allow for averaging and confidence interval 
#'        calculation.
#' @param output_directory A character string specifying the directory where the 
#'        plot should be saved. Defaults to "plots/".
#' @param rho_scale_range A numeric vector of length 2 defining the range for the 
#'        `rho1` and `rho2` axes. Defaults to `c(-1, 1)`.
#' @param scalestdevstep_range A numeric vector of length 2 defining the range for 
#'        the `scalestdevstep` axis. Defaults to `c(0.49, 0.51)`.
#' @param plot_all A logical value indicating whether to generate a combined plot 
#'        of all comparisons in the dataset. Defaults to `TRUE`.
#' @param alpha_value A numeric value controlling the transparency of points where 
#'        the SPSD condition is uncertain. Defaults to `1` (fully opaque).
#'
#' @return The function returns a 3D scatter plot generated using `plotly`. It also 
#'         saves the plot as an HTML file in the specified `output_directory`.
#'         
#' @import plotly
#' @import dplyr
#' @import htmlwidgets
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' data <- data.frame(
#'   rho1 = runif(100, -1, 1),
#'   rho2 = runif(100, -1, 1),
#'   scalestdevstep = runif(100, 0.49, 0.51),
#'   SPSD = rnorm(100),
#'   sigma = runif(100, 0, 1),
#'   comparison = sample(1:5, 100, replace = TRUE)
#' )
#' 
#' # Generate the plot
#' plot_spsdparameter(data, output_directory = "plots/")
#' }
#' 
#' @export
plot_spsdparameter <- function(data, output_directory = "plots/", rho_scale_range = c(min(data$rho1), max(data$rho1)), scalestdevstep_range = c(min(data$scalestdevstep), max(data$scalestdevstep)), alpha_value = 1) {

  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Step 1: Group by 'comparison' and calculate both the mean and CI for each group
  summarised_data <- data %>%
    group_by(rho1, rho2, scalestdevstep) %>%
    summarise(
      mean_SPSD = mean(SPSD),
      ci_lower_SPSD = mean(SPSD) - qt(0.975, df = n() - 1) * sd(SPSD) / sqrt(n()),  # 95% CI lower bound
      ci_upper_SPSD = mean(SPSD) + qt(0.975, df = n() - 1) * sd(SPSD) / sqrt(n()),  # 95% CI upper bound
      mean_sigma = mean(sigma),
      .groups = 'drop'
    )
  
  plot <- plot_ly()

  # Separate based on CI bounds
  data_spsd_uncertain <- subset(summarized_data, ci_lower_SPSD < 0)
  data_spsd_certain <- subset(summarized_data, ci_lower_SPSD >= 0)

  if (nrow(data_spsd_uncertain) > 0) {
      plot <- plot %>%
        add_markers(data = data_spsd_uncertain, 
                    x = ~rho1, 
                    y = ~rho2, 
                    z = ~scalestdevstep, 
                    marker = list(symbol = 'circle', color = "rgba(255, 0, 0, 0.5)", size = 5),  # Red and semi-transparent
                    name = "SPSD Uncertain",
                    hoverinfo = "text",
                    text = ~paste("rho1:", rho1, "<br>rho2:", rho2, "<br>scalestdevstep:", scalestdevstep, "<br>sigma:", mean_sigma, "<br>SPSD CI:", ci_lower_SPSD, "-", ci_upper_SPSD))
  }

    # Plot certain SPSD points with gradient based on sigma values
  if (nrow(data_spsd_certain) > 0) {
      plot <- plot %>%
        add_markers(data = data_spsd_certain, 
                    x = ~rho1, 
                    y = ~rho2, 
                    z = ~scalestdevstep, 
                    marker = list(symbol = 'circle', 
                                  color = ~mean_sigma, 
                                  colorscale = 'Viridis', 
                                  colorbar = list(title = "Sigma", len = 0.4),  
                                  size = 8),
                    name = "SPSD Certain",
                    hoverinfo = "text",
                    text = ~paste("rho1:", rho1, "<br>rho2:", rho2, "<br>scalestdevstep:", scalestdevstep, "<br>sigma:", mean_sigma, "<br>SPSD CI:", ci_lower_SPSD, "-", ci_upper_SPSD))
  }

  plot <- plot %>%
      layout(
        scene = list(
          xaxis = list(title = 'rho1', range = rho_scale_range),
          yaxis = list(title = 'rho2', range = rho_scale_range),
          zaxis = list(title = 'scalestdevstep', range = scalestdevstep_range)
        ),
        title = "Sigma Scatter Plot with SPSD Uncertainty",
        showlegend = TRUE
      )
  output_file <- file.path(output_directory, "sigma_scatter_with_CI.html")
  htmlwidgets::saveWidget(plot, file = output_file)
  message("Generated and saved scatter plot with CI for SPSD uncertainty.")
}



#' Bivariate Grid Plot with SPSD Classification and Adjustable Sample Size
#'
#' This function generates a bivariate scatter plot matrix to visualize the relationships 
#' between `rho1`, `rho2`, `scalestdevstep`, and `sigma`, classifying points based on 
#' whether they are SPSD or not. The plot includes scatter plots, correlation coefficients, 
#' and density plots for each variable. Points are colored based on their SPSD category.
#'
#' @param data A data frame containing the variables `rho1`, `rho2`, `scalestdevstep`, 
#'        `sigma`, and `SPSD`. The dataset should contain multiple observations of these 
#'        variables to generate meaningful bivariate plots.
#' @param output_directory A character string specifying the directory where the plot 
#'        should be saved. Defaults to "plots/".
#' @param sample_size A numeric value specifying the number of points to sample from the 
#'        dataset for plotting. Defaults to `3000`. If the dataset has fewer than 
#'        `sample_size` rows, the entire dataset will be used.
#'
#' @return The function returns a bivariate scatter plot matrix generated using `ggpairs`. 
#'         It also saves the plot as a PNG file in the specified `output_directory`.
#'         
#' @import dplyr
#' @import GGally
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' data <- data.frame(
#'   rho1 = runif(1000, -1, 1),
#'   rho2 = runif(1000, -1, 1),
#'   scalestdevstep = runif(1000, 0.49, 0.51),
#'   sigma = runif(1000, 0, 1),
#'   SPSD = rnorm(1000)
#' )
#' 
#' # Generate the bivariate plot, using the default sample size of 3000
#' plot_bivariate_grid(data, output_directory = "plots/")
#' }
#' 
#' @export
plot_bivariate_grid <- function(data, output_directory = "plots/", sample_size = 3000) {

  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Add a column classifying SPSD values
  plot_data <- data %>%
    dplyr::mutate(SPSD_category = ifelse(SPSD >= 0, "SPSD >= 0", "SPSD < 0")) %>%
    dplyr::select(rho1, rho2, scalestdevstep, sigma, SPSD_category)

  # Sample the data based on the specified sample_size
  if (sample_size < nrow(plot_data)) {
    plot_data_sample <- plot_data %>% sample_n(sample_size)
  } else {
    plot_data_sample <- plot_data
  }

  # Create custom color palette: one for SPSD >= 0 and one for SPSD < 0
  color_palette <- c("SPSD >= 0" = "blue", "SPSD < 0" = "red")

  # Custom function for scatter plots in lower triangle, colored by SPSD category
  lower_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(aes(color = SPSD_category), alpha = 0.7, size = 1.5) +
      scale_color_manual(values = color_palette) +
      theme(legend.position = "bottom")
  }

  # Create a grid of bivariate plots using ggpairs
  pairwise_plot <- ggpairs(
    plot_data_sample,
    mapping = aes(color = SPSD_category),  # Map color to SPSD category
    title = "Bivariate Scatter Plot Matrix",
    upper = list(continuous = wrap("cor", size = 4)),  # Add correlation in upper triangle
    lower = list(continuous = lower_fn),  # Custom lower function to color points by SPSD
    diag = list(continuous = "densityDiag")  # Add density plots on diagonal
  )

  # Display the plot
  print(pairwise_plot)

  # Save the plot as a PNG file
  output_file <- file.path(output_directory, "bivariate_plot_matrix.png")
  ggsave(output_file, pairwise_plot, width = 8, height = 8, dpi = 300)

  # Optionally, print a message to confirm plot generation
  message("Generated and saved bivariate scatter plot matrix.")
}

#' Plot Proportion of Intervals That Do Not Cover Zero (Pre-calculated)
#'
#' This function plots the proportion of intervals that do not cover zero across unique comparisons, 
#' assuming the calculations have already been done.
#'
#' @param data A dataframe that contains the proportion of intervals that do not cover zero.
#' @param comparison_col A string representing the name of the column that identifies the unique comparisons (default: "comparison").
#' @param proportion_col A string representing the name of the column that contains the calculated proportions (default: "proportion_intervals_dontcoverzero").
#'
#' @return A bar plot showing the proportion of intervals that do not cover zero for each unique comparison.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Assuming you have a dataframe called `results_df` with columns for comparison and proportion_intervals_dontcoverzero
#' # proportiondontcoverzerobars(results_df)
proportiondontcoverzerobars <- function(data, comparison_col = "comparison", proportion_col = "proportion_intervals_dontcoverzero", 
                                       ci_upper = "ninetyfive_ci_upper", ci_lower = "ninetyfive_ci_lower", outputdirectory = "./plots/") {
  
  # Check if the directory exists
  if (!dir.exists(outputdirectory)) {
    # Create the directory
    dir.create(outputdirectory)
    if (dir.exists(outputdirectory)) {
      cat("Directory created successfully!\n")
    } else {
      cat("Failed to create directory.\n")
    }
  } else {
    cat("Directory already present.\n")
  }
  
  # Create a new column to flag whether the 95% CI does not cover zero
  data$does_not_cover_zero <- ifelse(
    (data[[ci_upper]] > 0 & data[[ci_lower]] > 0) | (data[[ci_upper]] < 0 & data[[ci_lower]] < 0),
    "Does Not Cover Zero", "Covers Zero"
  )
  
  # Reorder the factor levels of the comparison column based on the proportion column
  data[[comparison_col]] <- factor(data[[comparison_col]], levels = data[[comparison_col]][order(data[[proportion_col]], decreasing = TRUE)])
  
  # Create the bar plot
  p <- ggplot(data, aes_string(x = comparison_col, y = proportion_col, fill = "does_not_cover_zero")) +
    geom_bar(stat = "identity", colour = "black", size = 0.5) +
    scale_fill_manual(values = c("Does Not Cover Zero" = "green", "Covers Zero" = "magenta")) +
    labs(
      title = "Proportion of Range Intervals That Do Not Cover Zero",
      x = "Comparison",
      y = "Proportion",
      fill = "95% CI Status"
    ) +
    theme_classic() +
    theme(
      # Increase font sizes for text elements
      text = element_text(size = 12, family = "Arial"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "top"
    ) +
    # Add horizontal dashed lines at specified y-values
    geom_hline(yintercept = c(0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975), linetype = "dashed", color = "grey50") +
    # Set y-axis limits
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.05), expand = c(0, 0)) +
    
    # Add text labels for values < 0.05 using geom_text
    geom_text(data = subset(data, data[[proportion_col]] == 0), 
              aes_string(label = proportion_col), 
              vjust = -0.5, color = "black", size = 3.5)  # Adjust position above bars
  
  # Save the plot in high resolution suitable for publications
  ggsave(paste0(outputdirectory, "proportion_coverage_plot.jpg"), plot = p, width = 15, height = 5, dpi = 300, units = "in")
  
  # Return the plot object
  return(p)
}

#' Plot Ridge Plot with Boxplots for RhoLower and RhoUpper per Taxa
#'
#' This function creates a professional-looking ridge plot showing the distribution of RhoLower and RhoUpper 
#' for each taxa, stratified by samples, and overlays boxplots for each ridge. It saves the plot in a "plots" directory, 
#' creating the directory if it does not exist, and also prints the plot to the screen.
#'
#' @param rhobounds A 3D array of dimensions [D x 2 x S], where D is the number of taxa,
#'   the second dimension contains RhoLower and RhoUpper values, and S is the number 
#'   of samples.
#' @param D Integer representing the number of taxa.
#' @param S Integer representing the number of samples.
#'
#' @return A ggplot object containing the ridge plot with boxplots.
#' 
#' @importFrom ggplot2 ggplot aes labs theme_classic theme element_text element_rect
#' @importFrom ggplot2 scale_fill_manual scale_color_manual ggsave after_stat geom_boxplot
#' @importFrom ggplot2 scale_y_discrete scale_x_continuous
#' @importFrom ggridges geom_density_ridges
#' @importFrom reshape2 melt
#' @keywords internal
plot_rho_ridges <- function(rhobounds, D, S) {
  
  # Reshape the rhobounds array into a long-format data frame
  melted_data <- reshape2::melt(rhobounds)
  colnames(melted_data) <- c("Taxa", "Bound", "Sample", "Rho")
  
  # Map the 'Bound' variable to 'RhoLower' and 'RhoUpper'
  melted_data$Bound <- factor(melted_data$Bound, levels = c(1, 2), labels = c("RhoLower", "RhoUpper"))
  
  # Reorder Taxa for plotting
  taxa_levels <- rev(unique(melted_data$Taxa))
  melted_data$Taxa <- factor(melted_data$Taxa, levels = taxa_levels)
  
  # Convert Taxa to numeric for adjusted positions
  melted_data$Taxa_numeric <- as.numeric(melted_data$Taxa)
  
  # Create adjusted positions for boxplots
  melted_data$Taxa_position <- melted_data$Taxa_numeric + ifelse(melted_data$Bound == "RhoLower", -0.2, 0.2)
  
  # Colors for bounds
  bound_colors <- c("RhoLower" = "#0072B2", "RhoUpper" = "#D55E00")
  
  # Create the ridge plot with boxplots
  plot <- ggplot2::ggplot(melted_data, ggplot2::aes(x = Rho, y = Taxa, fill = Bound)) +
    # Density ridges with increased transparency and more spacing
    ggridges::geom_density_ridges(
      ggplot2::aes(color = Bound),
      alpha = 0.4, scale = 0.7, rel_min_height = 0.01
    ) +
    # Overlay skinnier boxplots with adjusted positions
    ggplot2::geom_boxplot(
      ggplot2::aes(y = Taxa_position, group = interaction(Taxa, Bound)),
      width = 0.15, outlier.shape = NA, alpha = 0.6, color = "black"
    ) +
    # Labels and theme
    ggplot2::labs(
      title = "Distribution of RhoLower and RhoUpper for Each Taxa",
      x = "Rho Value",
      y = "Taxa",
      fill = "Rho Bound",
      color = "Rho Bound"
    ) +
    # Color palettes
    ggplot2::scale_fill_manual(values = bound_colors) +
    ggplot2::scale_color_manual(values = bound_colors) +
    # Adjust y-axis labels to display taxa names
    ggplot2::scale_y_discrete(
      expand = c(0.1, 0)
    ) +
    # Limit x-axis between -1 and 1 and set ticks every 0.1
    ggplot2::scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.1)
    ) +
    # Professional theme
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  
  # Print plot to the screen
  print(plot)
  
  # Check if "plots" directory exists, if not, create it
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  
  # Save the plot to the "plots" directory with white background
  ggplot2::ggsave(
    filename = "plots/correlation_ridge_plot.png",
    plot = plot,
    width = 12, height = 8, dpi = 300, bg = "white"
  )
  
  return(plot)
}







#' Plot Histogram with Density Overlay for Scale Standard Deviation Bounds
#'
#' This function creates a professional-looking histogram for the scale standard deviation 
#' (`scalestdev`) for each sample, overlaid with density plots. The histogram 
#' and density lines are colored differently for the lower and upper bounds.
#' Vertical lines indicating the means of each bound are added to the plot.
#' It saves the plot in a "plots" directory, creating the directory if it does not exist,
#' and also prints the plot to the screen.
#'
#' @param scalestdev A matrix of dimensions [S x 2], where S is the number of samples.
#'   The first column contains the lower bound of the standard deviation, and the 
#'   second column contains the upper bound.
#' @param S Integer representing the number of samples.
#'
#' @return A ggplot object containing the histogram with density overlays.
#' 
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs theme_classic theme element_line
#' @importFrom ggplot2 element_text scale_fill_manual scale_color_manual ggsave after_stat geom_vline
#' @importFrom ggplot2 annotate
#' @keywords internal
plot_scalestdev_histogram <- function(scalestdev, S) {
  
  # Prepare the data
  plot_data <- data.frame(
    ScaleSD = c(scalestdev[, 1], scalestdev[, 2]),
    Bound = factor(rep(c("Lower", "Upper"), each = S), levels = c("Lower", "Upper"))
  )
  
  # Calculate means for each bound
  mean_values <- aggregate(ScaleSD ~ Bound, data = plot_data, FUN = mean)
  
  # Create the plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ScaleSD, fill = Bound, color = Bound)) +
    # Histogram
    ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), 
                            bins = 30, position = "identity", alpha = 0.6, 
                            show.legend = TRUE, color = "black") +
    # Density plot with increased transparency
    ggplot2::geom_density(size = 1, adjust = 1.5, alpha = 0.2) +
    # Vertical lines at means
    ggplot2::geom_vline(data = mean_values, ggplot2::aes(xintercept = ScaleSD, color = Bound),
                        linetype = "dashed", size = 1) +
    # Annotate mean values
    ggplot2::annotate("text", x = mean_values$ScaleSD, y = Inf, label = paste0("Mean = ", round(mean_values$ScaleSD, 4)),
                      color = c("#1b9e77", "#d95f02"), angle = 90, vjust = -0.5, hjust = 1.1, size = 5) +
    # Labels and theme
    ggplot2::labs(title = "Histogram of Scale Standard Deviation Bounds",
                  x = "Scale Standard Deviation",
                  y = "Density",
                  fill = "Bound",
                  color = "Bound") +
    # Color palette
    ggplot2::scale_fill_manual(values = c("Lower" = "#0072B2", "Upper" = "#D55E00")) +
    ggplot2::scale_color_manual(values = c("Lower" = "#0072B2", "Upper" = "#D55E00")) +
    # Professional theme
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", size = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Print plot to the screen
  print(plot)
  
  # Check if "plots" directory exists, if not, create it
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  
  # Save the plot to the "plots" directory with white background
  ggplot2::ggsave(filename = "plots/scalestdev_histogram_plot.png", plot = plot, 
                  width = 10, height = 8, dpi = 300, bg = "white")
  
  return(plot)
}

#' Plot Posterior Density for All Taxa
#'
#' This internal function generates a density plot to visualize the posterior distributions 
#' of relative abundances for all taxa.
#' The plot is saved as a high-resolution PNG file.
#'
#' @param rWparaoriginal A 3D array containing posterior samples of relative abundances. 
#'                       Dimensions should be [taxa, samples, iterations].
#' @param file_name A string representing the name of the file to save the plot.
#'                  Default is "taxa_posterior_density.png".
#' @param width The width of the saved plot in inches. Default is 8.
#' @param height The height of the saved plot in inches. Default is 6.
#' @param dpi The resolution of the saved plot in dots per inch (dpi). Default is 300.
#'
#' @keywords internal
plot_posterior_density <- function(rWparaoriginal, save = TRUE, file_name = "taxa_posterior_density.png", width = 8, height = 8, dpi = 300) {
  
  # Get the number of taxa (D) and iterations (S)
  D <- dim(rWparaoriginal)[1]
  S <- dim(rWparaoriginal)[3]
  
  # Create a data frame from the array for density plotting
  plot_data <- data.frame(
    Iteration = rep(1:S, times = D),
    Taxa = rep(paste0("Taxa", 1:D), each = S),
    Value = as.vector(rWparaoriginal[1:D, 1, ])  # Select the first sample to plot posterior distributions
  )
  
  # Generate density plot
  p <- ggplot(plot_data, aes(x = Value, fill = Taxa, color = Taxa)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(
      title = "Posterior Density of Relative Abundance for All Taxa",
      x = "Relative Abundance (Log)",
      y = "Density",
      fill = "Taxa"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.box.margin = margin(0, 20, 0, 0),
      plot.margin = margin(5.5, 40, 5.5, 5.5)
    )
  
  if (save == TRUE){	
  # Save the plot with high resolution
  ggsave(filename = file_name, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  } 
  return(p)
}

#' Plot Posterior Boxplot for All Taxa
#'
#' This internal function generates a boxplot to visualize the posterior samples 
#' of relative abundances for all taxa.
#' The plot is saved as a high-resolution PNG file.
#'
#' @param rWparaoriginal A 3D array containing posterior samples of relative abundances. 
#'                       Dimensions should be [taxa, samples, iterations].
#' @param file_name A string representing the name of the file to save the plot.
#'                  Default is "taxa_posterior_boxplot.png".
#' @param width The width of the saved plot in inches. Default is 8.
#' @param height The height of the saved plot in inches. Default is 6.
#' @param dpi The resolution of the saved plot in dots per inch (dpi). Default is 300.
#'
#' @keywords internal
plot_posterior_boxplot <- function(rWparaoriginal, save = TRUE, file_name = "taxa_posterior_boxplot.png", width = 8, height = 6, dpi = 300) {
  
  # Get the number of taxa (D) and iterations (S)
  D <- dim(rWparaoriginal)[1]
  S <- dim(rWparaoriginal)[3]
  
  # Create a data frame for boxplot
  plot_data <- data.frame(
    Iteration = rep(1:S, times = D),
    Taxa = rep(paste0("Taxa", 1:D), each = S),
    Value = as.vector(rWparaoriginal[1:D, 1, ])  # Select the first sample to plot posterior distributions
  )
  
  # Generate the box plot
  p <- ggplot(plot_data, aes(x = Taxa, y = Value, fill = Taxa)) +
    geom_boxplot(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Posterior Boxplot of Relative Abundance for All Taxa",
      x = " ",
      y = "Relative Abundance (Log)",
      fill = "Taxa"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "none",  # Boxplots typically don't need legends
      plot.margin = margin(5.5, 40, 5.5, 5.5)
    )
  
  if (save == TRUE){	
  # Save the plot with high resolution
  ggsave(filename = file_name, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  } 
  return(p)
}

#' Plot Posterior Violin Plot for All Taxa
#'
#' This internal function generates a violin plot to visualize the posterior samples 
#' of relative abundances for all taxa.
#' The plot is saved as a high-resolution PNG file.
#'
#' @param rWparaoriginal A 3D array containing posterior samples of relative abundances. 
#'                       Dimensions should be [taxa, samples, iterations].
#' @param file_name A string representing the name of the file to save the plot.
#'                  Default is "taxa_posterior_violin.png".
#' @param width The width of the saved plot in inches. Default is 8.
#' @param height The height of the saved plot in inches. Default is 6.
#' @param dpi The resolution of the saved plot in dots per inch (dpi). Default is 300.
#'
#' @keywords internal
plot_posterior_violin <- function(rWparaoriginal, save = TRUE, file_name = "taxa_posterior_violin.png", width = 8, height = 8, dpi = 300) {
  
  # Get the number of taxa (D) and iterations (S)
  D <- dim(rWparaoriginal)[1]
  S <- dim(rWparaoriginal)[3]
  
  # Create a data frame for violin plot
  plot_data <- data.frame(
    Iteration = rep(1:S, times = D),
    Taxa = rep(paste0("Taxa", 1:D), each = S),
    Value = as.vector(rWparaoriginal[1:D, 1, ])  # Select the first sample to plot posterior distributions
  )
  
  # Generate the violin plot
  p <- ggplot(plot_data, aes(x = Taxa, y = Value, fill = Taxa)) +
    geom_violin(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Posterior Violin Plot of Relative Abundance for All Taxa",
      x = " ",
      y = "Relative Abundance (Log)",
      fill = "Taxa"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "none",
      plot.margin = margin(5.5, 40, 5.5, 5.5)
    )
  
  if (save == TRUE){	
  # Save the plot with high resolution
  ggsave(filename = file_name, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  } 
  return(p)
}

plot_posterior_samples <- function(rWparaoriginal, save = TRUE, file_name = "taxa_posterior_samples.png", width = 8, height = 8, dpi = 300) {
  
  # Get the number of taxa (D)
  D <- dim(rWparaoriginal)[1]
  
  # Create a data frame from the array for easier plotting
  plot_data <- data.frame(
    Sample = rep(1:dim(rWparaoriginal)[2], times = D),
    Taxa = rep(paste0("Taxa", 1:D), each = dim(rWparaoriginal)[2]),
    Value = as.vector(rWparaoriginal[1:D, , 1])  # You can change this to select the appropriate slice of data
  )
  
  # Generate the plot with the legend outside
  p <- ggplot(plot_data, aes(x = Sample, y = Value, color = Taxa)) +
    geom_line(size = 1) +
    theme_minimal() +
    labs(
      title = "Variation Across Taxa Posterior Samples",
      x = "Sample",
      y = "Relative Abundance (Log)",
      color = "Taxa"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right",          # Move legend to the right
      legend.box = "vertical",            # Arrange legend vertically
      legend.box.margin = margin(0, 20, 0, 0), # Add space between the plot and the legend
      plot.margin = margin(5.5, 40, 5.5, 5.5)  # Ensure enough space for the legend outside the plot
    )

  if (save == TRUE){	
  # Save the plot with high resolution
  ggsave(filename = file_name, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  } 
	
  return(p)
}

#' Combine All Posterior Plots into a Grid
#'
#' This internal function generates and combines the density plot, boxplot, violin plot,
#' and line plot of posterior samples for all taxa into a single grid layout.
#' The combined plot is saved as a high-resolution PNG file.
#'
#' @param rWparaoriginal A 3D array containing posterior samples of relative abundances. 
#'                       Dimensions should be [taxa, samples, iterations].
#' @param file_name A string representing the name of the file to save the combined plot.
#'                  Default is "combined_posterior_plots.png".
#' @param width The width of the saved plot in inches. Default is 16.
#' @param height The height of the saved plot in inches. Default is 12.
#' @param dpi The resolution of the saved plot in dots per inch (dpi). Default is 300.
#'
#' @keywords internal
#' @import ggplot2
#' @import gridExtra
combine_posterior_plots <- function(rWparaoriginal, file_name = "combined_posterior_plots.png", width = 16, height = 16, dpi = 300) {
  
  # Generate individual plots
  p1 <- plot_posterior_density(rWparaoriginal, save=FALSE)   # Density plot
  p2 <- plot_posterior_boxplot(rWparaoriginal, save=FALSE)   # Box plot
  p3 <- plot_posterior_violin(rWparaoriginal, save=FALSE)    # Violin plot
  p4 <- plot_posterior_samples(rWparaoriginal, save=FALSE)   # Line plot

    # Convert ggplot objects to grobs
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(p3)
  g4 <- ggplotGrob(p4)
	
  # Arrange all plots into a 2x2 grid layout
  combined_plot <- grid.arrange(g1, g2, g3, g4, ncol = 2)
  
  # Save the combined plot as a PNG
  ggsave(filename = file_name, plot = combined_plot, width = width, height = height, dpi = dpi, bg = "white")
  
  return(combined_plot)
}

#' Plot Rho Correlations for Taxa vs Flow Data, Covariance Heatmap, and Print Scale SD
#'
#' This function calculates the true rho correlations between taxa and flow data, 
#' generates a plot of these correlations, a pheatmap of covariances, and prints the standard deviation 
#' of the log-transformed flow data. The correlation plot and covariance heatmap are displayed side by side.
#'
#' @param rdat Data frame. The resampled abundance data from `simulate_prepost` with taxa counts.
#' @param flow_data Data frame. The flow cytometry measurements from `simulate_prepost` with flow values.
#' @param dat Data frame. The true abundance data from `simulate_prepost` to calculate covariances.
#' @param file_path Character. File path to save the grid of plots as a high-quality PNG image. Default is `"rho_and_covariances_plot.png"`.
#'
#' @return None. The function produces a correlation plot, a covariance heatmap, and prints the standard deviation.
#'
#' @importFrom stats cor sd cov
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot axis grid abline text
#' @import pheatmap
#' @import gridExtra
#' @export
#'
#' @examples
#' \dontrun{
#'   # Generate correlation plot and covariance heatmap
#'   plot_true_abundances(rdat, flow_data, dat, file_path = "true_abundances_plot.png")
#' }
plot_true_abundances <- function(flow_data, dat, file_path = "true_abundances_plot.png") {
  
  # Transpose rdat to have taxa as rows
  taxa_data <- t(dat[,-1])
  
  # Remove sample names from flow_data
  flow_data_numeric <- flow_data[,-1]
  
  # Calculate true rho correlations between each taxon and flow cytometry data
  truerhocorrelation <- numeric(nrow(taxa_data))
  for (i in 1:nrow(taxa_data)) {
    truerhocorrelation[i] <- cor(log(taxa_data[i,]), log(flow_data_numeric))
  }
  
  # Calculate standard deviation (scale) for log-transformed flow_data
  scale_sd_flow_data <- sd(log(flow_data_numeric))
  # Print the calculated scale SD
  print(paste("Scale standard deviation (SD) of flow data: ", scale_sd_flow_data))
  
  # Extract the taxa columns (excluding the first column which is 'Condition')
  truecorrelations <- cor(t(taxa_data))  # Calculate correlation matrix
  truecovariances <- cov(t(taxa_data))  # Calculate covariance matrix
	
  # Create forest plot data for covariances
  forest_plot_data <- data.frame(
    comparison = character(),
    ninetyfive_ci_lower = numeric(),
    ninetyfive_ci_upper = numeric(),
    minsigma_absolute_minimum_covariance = numeric(),
    maxsigma_absolute_maximum_covariance = numeric(),
    p_value = numeric()  # Placeholder
  )

  # Loop through the upper triangle of the covariance matrix to extract comparisons
  for (i in 1:(ncol(truecovariances) - 1)) {
    for (j in (i + 1):ncol(truecovariances)) {
      # Extract the covariance value
      cov_value <- truecovariances[i, j]
      
      # Create the comparison name (e.g., "Taxa1:Taxa2")
      comparison_name <- paste0(colnames(truecovariances)[i], ":", colnames(truecovariances)[j])
      
      # Append this comparison to the data frame
      forest_plot_data <- rbind(forest_plot_data, data.frame(
        comparison = comparison_name,
        ninetyfive_ci_lower = cov_value,  # True value for lower CI
        ninetyfive_ci_upper = cov_value,  # True value for upper CI
        minsigma_absolute_minimum_covariance = cov_value,  # Min sigma is the true value
        maxsigma_absolute_maximum_covariance = cov_value,  # Max sigma is the true value
        p_value = NA  # Placeholder, can be removed if not needed
      ))
    }
  }

  # Create the forest plot using ggplot
  taxa_labels <- rownames(taxa_data)
  forest_plot <- ggplot(data = data.frame(taxa = taxa_labels, rho = truerhocorrelation), aes(x = rho, y = taxa)) +
	  geom_point(color = "steelblue", size = 4) +
	  geom_segment(aes(x = 0, xend = rho, y = taxa, yend = taxa), color = "steelblue", size = 1.2) +
	  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
	  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +  # Set x-axis limits and breaks
	  labs(x = "Correlation") +  # Set x-axis label
	  theme_minimal() +
	  theme(
	    axis.title.y = element_blank(),
	    axis.text.y = element_text(size = 12),
	    plot.title = element_text(hjust = 0.5, size = 15)
	  ) +
	  ggtitle("True Correlations of Taxa and Scale")
  
  # Create the correlation heatmap using pheatmap
  correlation_heatmap <- pheatmap::pheatmap(
	  truecorrelations, 
	  main = "Correlation Matrix of Taxa", 
	  color = colorRampPalette(c("blue", "white", "red"))(100),
	  breaks = seq(-1, 1, length.out = 101),  # Ensures the color range goes from -1 to 1
	  border_color = NA, 
	  silent = TRUE
	)
  # Convert to a grob
  forest_plot_grob <- ggplotGrob(forest_plot)
  correlation_grob <- grid::grid.grabExpr(grid::grid.draw(correlation_heatmap$gtable))

  # Save the combined plot as a PNG file
  png(filename = file_path, width = 15, height = 8, units = "in", res = 300)
  
  # Arrange the forest plot and correlation heatmap side by side
  gridExtra::grid.arrange(forest_plot_grob, correlation_grob, ncol = 2)
  
  # Close the device to save the file
  dev.off()

  return(forest_plot_data)
}

#' Generate a professional network plot from confidence interval data
#'
#' This function generates a professional network plot from taxa labels and confidence interval data. 
#' It provides an option to save the plot (default in PNG format) and returns the plot object for further 
#' customization. The network plot adjusts edge color based on the sign of the covariance inferred from the 
#' confidence intervals, edge width based on the CI range (to reflect magnitude), and transparency reflects 
#' uncertainty. Edges are grey if the confidence interval covers zero.
#'
#' @param results A dataframe containing the confidence interval information. The dataframe should include columns:
#'   - `taxa1`: The first taxa in the pair.
#'   - `taxa2`: The second taxa in the pair.
#'   - `ninetyfive_ci_lower`: Lower bound of the 95% CI for covariance.
#'   - `ninetyfive_ci_upper`: Upper bound of the 95% CI for covariance.
#' @param file_name Character; the name of the file to save the plot. Default is 'network_plot.png'.
#' @param save_plot Logical; if TRUE, saves the plot to file. Default is TRUE.
#' 
#' @return The network plot object for further customization.
#'
#' @import igraph
#' @import qgraph
#'
#' @examples
#' prism.network(results)

prism.network <- function(results, file_name = "network_plot.png", save_plot = TRUE) {
  # Filter rows with non-NA taxa
  filtered_results <- results[!is.na(results$taxa1) & !is.na(results$taxa2), ]
  
  # Initialize an adjacency matrix for the network
  taxa <- unique(c(filtered_results$taxa1, filtered_results$taxa2))
  adj_matrix <- matrix(0, nrow = length(taxa), ncol = length(taxa), 
                       dimnames = list(taxa, taxa))
  
  # Set edge color, width, and transparency based on the 95% CI for covariances
  edge_colors <- c()
  edge_widths <- c()
  
  for (i in 1:nrow(filtered_results)) {
    taxa1 <- filtered_results$taxa1[i]
    taxa2 <- filtered_results$taxa2[i]
    
    # Extract confidence interval
    ci_lower_cov <- filtered_results$ninetyfive_ci_lower[i]
    ci_upper_cov <- filtered_results$ninetyfive_ci_upper[i]
    
    # Determine sign based on CI
    if (ci_lower_cov <= 0 & ci_upper_cov >= 0) {
      # Confidence interval covers zero, so edge should be grey
      edge_color <- "grey"
    } else {
      # Positive or negative covariance based on the CI bounds
      edge_color <- ifelse(ci_lower_cov > 0, "blue", "red")
    }
    
    # Edge width based on CI range (larger range -> thinner line)
    ci_range_cov <- ci_upper_cov - ci_lower_cov
    edge_width <- 1 / ci_range_cov  # Inverse of CI range for width
    
    # Set edge colors and widths
    edge_colors <- c(edge_colors, edge_color)
    edge_widths <- c(edge_widths, edge_width)
    
    # Update adjacency matrix with the inferred covariance strength
    adj_matrix[taxa1, taxa2] <- edge_width  # Magnitude based on CI width
    adj_matrix[taxa2, taxa1] <- edge_width  # Symmetric matrix
  }
  
  # Generate the network plot using qgraph
  plot_object <- qgraph(adj_matrix, 
                        layout = "spring",            # Spring layout for professional appearance
                        labels = TRUE,                # Add taxa labels
                        edge.color = edge_colors,     # Set edge colors based on CI sign
                        edge.width = edge_widths,     # Set edge widths based on CI range
                        posCol = c("#009900", "darkgreen"),  # Positive covariance color
                        negCol = c("#BF0000", "red"),        # Negative covariance color
                        unCol = "#808080",                   # Grey for edges where CI covers zero
                        trans = TRUE,                        # Enable transparency based on edge weight
                        fade = TRUE,                         # Enable fading based on edge weight
                        esize = 15 * exp(-length(taxa) / 90) + 1,  # Scalar for edge size
                        vsize = 8,                    # Set larger node size for clarity
                        borders = FALSE,              # Remove node borders for cleaner presentation
                        label.cex = 1.5)              # Increase label font size for clarity
  
  # Save the plot if save_plot is TRUE
  if (save_plot) {
    png(file_name, width = 1200, height = 1200, res = 300)  # High resolution PNG
    print(plot_object)  # Save the plot to file
    dev.off()           # Close the PNG device
  }
  
  # Return the plot object for further customization
  return(plot_object)
}

# Example usage with the 'results' dataframe
# prism.network(results)

