#' Estimate Covariance from Sequence Count Data with CLR/ILR Transformations
#'
#' This function takes in a sequence count dataset and runs various covariance estimation functions
#' commonly used in bioinformatics, including Banocc, SparCC, SPIEC-EASI, CCLasso, and Proportionality.
#' You can also choose to apply CLR or ILR transformations to close the data.
#'
#' @param count_data A matrix or data frame of sequence counts, where rows represent taxa and columns represent samples.
#' @param normalize Logical, indicating whether the count data should be normalized by column sums. Defaults to TRUE.
#' @param transformation Character, specifying the type of transformation to apply: "CLR" (Centered Log-Ratio),
#'        "ILR" (Isometric Log-Ratio), or "none". Defaults to "none".
#' @return A list containing the estimated covariance matrices from each method: Banocc, SparCC, SPIEC-EASI, CCLasso, and Proportionality.
#' @examples
#' \dontrun{
#' # Load example data (replace with your actual data)
#' count_data <- read.csv("your_sequence_data.csv", row.names = 1)
#'
#' # Run covariance estimation methods with CLR transformation
#' result <- prism.covariance_comparison(count_data, transformation = "CLR")
#'
#' # Access the covariance matrices
#' result$cov_matrix_banocc
#' result$cov_matrix_sparcc
#' result$cov_matrix_spiec_easi
#' result$cov_matrix_cclasso
#' result$cov_matrix_proportionality
#' }
#' @import BANOVA
#' @import SpiecEasi
#' @import propr
#' @import phyloseq
#' @import Matrix
#' @import igraph
#' @import compositions
#' @export
prism.covariance_comparison <- function(count_data, normalize = TRUE, transformation = "none") {
  
  # Normalization function (Common pre-processing step)
  normalize_counts <- function(count_matrix) {
    return(sweep(count_matrix, 2, colSums(count_matrix), "/"))
  }
  
  # Apply normalization if needed
  if (normalize) {
    normalized_data <- normalize_counts(as.matrix(count_data))
  } else {
    normalized_data <- as.matrix(count_data)
  }
  
  # Apply chosen transformation
  if (transformation == "CLR") {
    transformed_data <- clr(normalized_data)  # Centered Log-Ratio transformation
  } else if (transformation == "ILR") {
    transformed_data <- ilr(normalized_data)  # Isometric Log-Ratio transformation
  } else {
    transformed_data <- normalized_data  # No transformation
  }

  ### Banocc ###
  # Convert to phyloseq object for Banocc
  ps <- phyloseq(otu_table(transformed_data, taxa_are_rows = TRUE))
  # Define the prior information (you can adjust these as needed)
  priors <- list(
    Sigma = diag(ncol(ps)),  # The covariance matrix prior
    alpha = 0.5,  # Dirichlet prior on proportions
    mu = rep(0, ncol(ps))  # Prior mean for log-ratios
  )
  
  # Run Banocc model
  banocc_results <- banocc::run_banocc(
    ps,
    priors = priors,
    num_iterations = 10000,  # Number of iterations
    burnin = 5000  # Number of burn-in iterations
  )
  cov_matrix_banocc <- banocc_results$posterior$Sigma

  ### SparCC ###
  # Run SparCC using the prism_SparCC_count function
  sparcc_result <- SparCC_count(x = transformed_data)
  
  # Extract the covariance matrix from the result
  cov_matrix_sparcc <- sparcc_result$cov.w

  ### SPIEC-EASI ###
  spiec_easi_result <- spiec.easi(transformed_data, method = "mb", lambda.min.ratio = 1e-2, nlambda = 20)
  cov_matrix_spiec_easi <- cov2cor(spiec_easi_result$refit$cov)

  ### CCLasso ###
  cclasso_result <- cclasso(transformed_data)
  cov_matrix_cclasso <- cclasso_result$Theta

  ### Proportionality ###
  prop_result <- propr(transformed_data, metric = "rho")
  cov_matrix_proportionality <- prop_result@matrix

  # Return results as a list
  return(list(
    cov_matrix_banocc = cov_matrix_banocc,
    cov_matrix_sparcc = cov_matrix_sparcc,
    cov_matrix_spiec_easi = cov_matrix_spiec_easi,
    cov_matrix_cclasso = cov_matrix_cclasso,
    cov_matrix_proportionality = cov_matrix_proportionality
  ))
}
