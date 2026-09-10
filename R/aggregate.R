#' Extract one draw-array entry per pair as a pairs x draws matrix
#' @keywords internal
.pair_draws <- function(arr, i, j) {
  S <- dim(arr)[3]
  t(vapply(seq_along(i), function(k) arr[i[k], j[k], ], numeric(S)))
}

#' Off-diagonal (i < j) pair indices for a D x D matrix, column-major order
#' @keywords internal
.offdiag_pairs <- function(D) {
  idx <- which(upper.tri(matrix(0, D, D)), arr.ind = TRUE)
  idx[order(idx[, 2], idx[, 1]), , drop = FALSE]
}

#' Upper-triangle-including-diagonal (i <= j) pair indices, column-major order
#' @keywords internal
.upper_pairs <- function(D) {
  idx <- which(upper.tri(matrix(0, D, D), diag = TRUE), arr.ind = TRUE)
  idx[order(idx[, 2], idx[, 1]), , drop = FALSE]
}

#' Aggregate accepted bootstrap draws into CI matrices and pairwise tables
#'
#' `ci_lower`/`ci_upper` are the empirical `ci_alpha/2` and
#' `1 - ci_alpha/2` quantiles of the accepted per-draw lower/upper bound
#' arrays -- for a single deterministic draw (`S_eff == 1`), the quantile
#' of one value is that value, so this reduces automatically to the
#' identification interval with no special case needed. `p_value` is the
#' one genuine special case: with only one draw there is no bootstrap
#' distribution to test against the null region, so it is reported as
#' `NA` rather than a degenerate one-sample value.
#'
#' @param draws Output of `.prism_bootstrap()`: `list(lower, upper, rel)`,
#'   each a D x D x S_eff array over accepted draws.
#' @param feature_names Character vector of length D, or `NULL`.
#' @param ci_alpha Tail probability for CI endpoints.
#' @param delta Half-width of the pairwise null region used for p-values.
#' @param p_adjust_method Passed to `stats::p.adjust()`.
#' @return `list(ci_lower, ci_upper, pairwise, pairwise_rel)`.
#' @keywords internal
.prism_aggregate <- function(draws, feature_names = NULL, ci_alpha = 0.05,
                              delta = 0.1, p_adjust_method = "BH") {
  D <- dim(draws$lower)[1]
  S_eff <- dim(draws$lower)[3]
  ci_alpha <- .assert_scalar_finite(ci_alpha, "ci_alpha", lower = 0, upper = 1, upper_inclusive = FALSE)
  delta <- .assert_scalar_finite(delta, "delta", lower = 0)
  if (!p_adjust_method %in% stats::p.adjust.methods) {
    stop("p_adjust_method must be one of: ", paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)
  }
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(D))

  ci_lower <- apply(draws$lower, c(1, 2), stats::quantile, probs = ci_alpha / 2, names = FALSE)
  ci_upper <- apply(draws$upper, c(1, 2), stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)
  dimnames(ci_lower) <- dimnames(ci_upper) <- list(feature_names, feature_names)

  off <- .offdiag_pairs(D)
  i <- off[, 1]
  j <- off[, 2]

  if (S_eff == 1L) {
    p_value <- rep(NA_real_, length(i))
  } else {
    lower_draws <- .pair_draws(draws$lower, i, j)
    upper_draws <- .pair_draws(draws$upper, i, j)
    pos_count <- rowSums(lower_draws > delta)
    neg_count <- rowSums(upper_draws < -delta)
    p_plus <- (1 + pos_count) / (1 + S_eff)
    p_minus <- (1 + neg_count) / (1 + S_eff)
    p_value <- pmin(1, 2 * pmin(p_plus, p_minus))
  }
  q_value <- stats::p.adjust(p_value, method = p_adjust_method)

  pairwise <- data.frame(
    i = i, j = j,
    feature_i = feature_names[i], feature_j = feature_names[j],
    ci_lower = ci_lower[cbind(i, j)], ci_upper = ci_upper[cbind(i, j)],
    p_value = p_value, q_value = q_value,
    stringsAsFactors = FALSE
  )

  full <- .upper_pairs(D)
  fi <- full[, 1]
  fj <- full[, 2]
  rel_hat_mat <- apply(draws$rel, c(1, 2), mean)
  rel_ci_lower_mat <- apply(draws$rel, c(1, 2), stats::quantile, probs = ci_alpha / 2, names = FALSE)
  rel_ci_upper_mat <- apply(draws$rel, c(1, 2), stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)

  pairwise_rel <- data.frame(
    i = fi, j = fj,
    feature_i = feature_names[fi], feature_j = feature_names[fj],
    rel_hat = rel_hat_mat[cbind(fi, fj)],
    rel_ci_lower = rel_ci_lower_mat[cbind(fi, fj)],
    rel_ci_upper = rel_ci_upper_mat[cbind(fi, fj)],
    stringsAsFactors = FALSE
  )
  pairwise_rel$rel_ci_width <- pairwise_rel$rel_ci_upper - pairwise_rel$rel_ci_lower

  list(ci_lower = ci_lower, ci_upper = ci_upper, pairwise = pairwise, pairwise_rel = pairwise_rel)
}
