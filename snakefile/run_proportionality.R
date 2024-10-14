# scripts/run_proportionality.R

library(propr)
library(dplyr)
library(tidyr)
library(parallel)

# Load data
load('data/abundances.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Normalize counts
normalize_counts <- function(count_matrix) {
  return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
}
normalized_data <- normalize_counts(as.matrix(Y))

# Run Proportionality analysis
proportionality_results <- list()
metrics <- c("rho", "phi", "phs")
for (metric in metrics) {
  prop_obj <- propr(
    t(normalized_data), metric = metric, 
    ivar = "clr", alpha = 0, p = 100
  )
  pr_obj <- updateCutoffs(
    prop_obj, tails = ifelse(metric == "rho", "both", "right"), 
    ncores = detectCores() - 1
  )
  proportionality_results[[metric]] <- pr_obj
}

# Save results
save(proportionality_results, file = "results/proportionality_results.RData")
