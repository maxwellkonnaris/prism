# scripts/run_banocc.R

library(banocc)
library(phyloseq)
library(rstan)
library(parallel)

# Load data
load('data/abundances.RData')

Y = abundances$Y
Y = t(Y) # Transpose Y

# Prepare data for BAnOCC
normalize_counts <- function(count_matrix) {
  return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
}
normalized_data <- normalize_counts(as.matrix(Y))
ps <- phyloseq::phyloseq(otu_table(normalized_data, taxa_are_rows = TRUE))
ps <- t(ps)

# ... (include your BAnOCC code here, adapted to work as a standalone script)

# Save results
saveRDS(banoccresults, file = "results/banocc_sensitivity_results.rds")
