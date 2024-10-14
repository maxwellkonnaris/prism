# scripts/simulate_data.R

# Load necessary libraries
if (!requireNamespace("PRISM", quietly = TRUE)) {
  devtools::install_github("maxwellkonnaris/prism", force=TRUE)
}
library(PRISM)

# ... (rest of your code from simulate_data.R)

# Save the abundances object and others to RData files
save(abundances, file = "data/abundances.RData")
save(externalscalemeasurements, file = "data/externalscalemeasurements.RData")
save(forest_plot_data, file = "data/forest_plot_data.RData")
