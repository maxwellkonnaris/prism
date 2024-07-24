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
