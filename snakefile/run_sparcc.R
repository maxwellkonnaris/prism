# scripts/run_sparcc.R

library(PRISM)

# Load data
load('data/abundances.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Run SparCC
sparcc_result <- PRISM::SparCC_count(x = t(Y))

# Prepare results as in your method comparison function
# ... (process sparcc_result to match the expected format)

# Save results
saveRDS(sparccresults, file = "results/sparcc_results.rds")
