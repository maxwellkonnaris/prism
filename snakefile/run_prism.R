# scripts/run_prism.R

library(PRISM)

# Load data
load('data/abundances.RData')
load('data/externalscalemeasurements.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Run PRISM method
prism_results <- PRISM::prism.covariance(
  Y = Y, 
  S = 4000, 
  uncertaintydistribution = "multinomiallognormal",
  externalscalemeasurements = externalscalemeasurements, 
  algorithm = "GRID_SEARCH", 
  outputdirectory = "results/", 
  prefix = "simulated_5000_", 
  pvalue = FALSE, 
  logfile = "results/log_prismcovariance.txt"
)

# Save results
saveRDS(prism_results$final_results, file = "results/prism_results.rds")
