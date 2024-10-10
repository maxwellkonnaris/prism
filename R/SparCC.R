
#' Title: SparCC for Count Data
#' 
#' Description: This function estimates covariance and correlation for count data using the SparCC algorithm.
#' 
#' @param x A nxp count data matrix where rows are samples and columns are variables.
#' @param imax Resampling times from the posterior distribution (default: 20).
#' @param kmax Maximum iteration steps for SparCC (default: 10).
#' @param alpha Threshold for strong correlation (default: 0.1).
#' @param Vmin Minimal variance if negative variance appears (default: 1e-4).
#' 
#' @return A list with two elements: `cov.w` for covariance estimation and `cor.w` for correlation estimation.
#' 
#' @examples
#' result <- SparCC_count(x = matrix_data)
#' 
#' @import gtools
#' @export
# Adapted from SparCC.R by Fang Huaying (Peking University)
# Original code: https://github.com/huayingfang/CCLasso
# This adaptation is for use within the PRISM package
################################################################################
# File: SparCC.R
# Aim : SparCC 
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 11/12/2014
#-------------------------------------------------------------------------------
# SparCC for counts known
#   function: SparCC_count
#   input:
#          x ------ nxp count data matrix, row is sample, col is variable
#       imax ------ resampling times from posterior distribution. default 20
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov_w ------ covariance estimation
#      cor_w ------ correlation estimation
#
require(gtools);
SparCC_count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # dimension for w (latent variables)
  p <- ncol(x);
  n <- nrow(x);
  # posterior distribution (alpha)
  x <- x + 1;
  # store generate data
  y <- matrix(0, n, p);
  # store covariance/correlation matrix
  cov_w <- cor_w <- matrix(0, p, p);
  indLow <- lower.tri(cov_w, diag = T);
  # store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x) 
      gtools::rdirichlet(n = 1, alpha = x)));
    # estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # store variance/correlation only low triangle 
    covs[, i] <- cov_cor$cov_w[indLow];
    cors[, i] <- cov_cor$cor_w[indLow];
  }
  # calculate median for several posterior samples
  cov_w[indLow] <- apply(covs, 1, median); 
  cor_w[indLow] <- apply(cors, 1, median);
  #
  cov_w <- cov_w + t(cov_w);
  diag(cov_w) <- diag(cov_w) / 2;
  cor_w <- cor_w + t(cor_w);
  diag(cor_w) <- 1;
  # 
  # Maxwell Edits to Sparcc below -- to extract uncertainty estimates for the model above:
  # Calculate 95% CI (2.5% and 97.5% quantiles)
  cor_ci_lower <- matrix(0, p, p);
  cor_ci_upper <- matrix(0, p, p);
  cov_ci_lower <- matrix(0, p, p);
  cov_ci_upper <- matrix(0, p, p);
  #
  # Correlation confidence intervals
  cor_ci_lower[indLow] <- apply(cors, 1, quantile, probs = 0.025);
  cor_ci_upper[indLow] <- apply(cors, 1, quantile, probs = 0.975);
  cor_ci_lower <- cor_ci_lower + t(cor_ci_lower);
  cor_ci_upper <- cor_ci_upper + t(cor_ci_upper);
  diag(cor_ci_lower) <- diag(cor_ci_upper) <- 1;
  #
  # Covariance confidence intervals
  cov_ci_lower[indLow] <- apply(covs, 1, quantile, probs = 0.025);
  cov_ci_upper[indLow] <- apply(covs, 1, quantile, probs = 0.975);
  cov_ci_lower <- cov_ci_lower + t(cov_ci_lower);
  cov_ci_upper <- cov_ci_upper + t(cov_ci_upper);
  diag(cov_ci_lower) <- diag(cov_ci_upper) <- diag(cov_w) / 2;
  #
  return(list(
    cov_w = cov_w, cor_w = cor_w, 
    cor_ci_lower = cor_ci_lower, cor_ci_upper = cor_ci_upper, 
    cov_ci_lower = cov_ci_lower, cov_ci_upper = cov_ci_upper
  ));
}
#-------------------------------------------------------------------------------
# SparCC for fractions known
#   function: SparCC.frac
#   input:
#          x ------ nxp fraction data matrix, row is sample, col is variable
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov_w ------ covariance estimation
#      cor_w ------ correlation estimation
SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # Log transformation
  x <- log(x);
  p <- ncol(x);
  # T0 = var(log(xi/xj)) variation matrix
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC  
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  #cor_w <- (outer(var.w, var.w, "+") - T0 ) / 
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  Is <- sqrt(1/var.w);
  cor_w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor_w[cor_w <= - 1] <- - 1; 
  cor_w[cor_w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1; 
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;  
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor_w <- cor_w;
    # Remove diagonal
    diag(curr_cor_w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor_w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor_w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor_w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor_w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;    
      Is <- sqrt(1/var.w);
      cor_w <- (var.w + rep(var.w, each = p) - T0) * 
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor_w[cor_w <= - 1] <- - 1;
      cor_w[cor_w >= 1] <- 1;
    }
    else {
      break;
    }
    # 
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov_w <- cor_w * Is * rep(Is, each = p);
  #
  return(list(cov_w = cov_w, cor_w = cor_w));
}
#-------------------------------------------------------------------------------
