# scripts/generate_plots.R

library(PRISM)

# Load results
covarianceresults <- readRDS('results/covariance_comparisons.rds')

# Generate forest plot
PRISM::prism.forestplot(
  covarianceresults, 
  save = "png", 
  filename = "plots/forest_plot", 
  color_y_axis_by_ci = TRUE
)

# Optionally, generate circle network plot
PRISM::prism.circlenetwork(
  covarianceresults, 
  metric = "covariance", 
  combine_plots = TRUE, 
  filename = "plots/circle_network"
)
