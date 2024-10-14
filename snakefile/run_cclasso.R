# scripts/run_cclasso.R

library(PRISM)

# Load data
load('data/abundances.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Run CClasso
cclasso_result <- PRISM::cclasso(
  t(Y), counts = TRUE, pseudo = 0.5, k_cv = 10, 
  lam_int = c(1e-4, 10), k_max = 100, n_boot = 100
)

# Prepare results as in your method comparison function
# ... (process cclasso_result to match the expected format)

# Save results
saveRDS(cclassoresults, file = "results/cclasso_results.rds")
