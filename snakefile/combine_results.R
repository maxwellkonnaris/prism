# scripts/combine_results.R

library(dplyr)
library(tidyr)
library(PRISM)

# Load results
load('data/forest_plot_data.RData')
prism_results <- readRDS('results/prism_results.rds')
banocc_results <- readRDS('results/banocc_sensitivity_results.rds')
spieceasi_results <- readRDS('results/spieceasi_results.rds')
sparcc_results <- readRDS('results/sparcc_results.rds')
cclasso_results <- readRDS('results/cclasso_results.rds')
load('results/proportionality_results.RData')

# Combine results into a list
covarianceresults <- list(
  PRISM = prism_results,
  BAnOCC = banocc_results,
  SpiecEasi = spieceasi_results,
  SparCC = sparcc_results,
  CClasso = cclasso_results,
  Propr = proportionality_results,
  TrueAbundances = forest_plot_data
)

# Save combined results
saveRDS(covarianceresults, file = "results/covariance_comparisons.rds")
