#' Extract one draw-array entry per pair as a pairs x draws matrix
#'
#' `vapply()` simplifies its result to a plain vector rather than a
#' matrix when `S == 1` (each call returns one value), so `t()` alone
#' would transpose the wrong way in that case; reshaping explicitly via
#' `matrix(..., byrow = TRUE)` handles `S == 1` and `S > 1` uniformly.
#' @keywords internal
.pair_draws <- function(arr, i, j) {
  S <- dim(arr)[3]
  matrix(vapply(seq_along(i), function(k) arr[i[k], j[k], ], numeric(S)), nrow = length(i), ncol = S, byrow = TRUE)
}

#' Upper-triangle-including-diagonal (i <= j) pair indices, column-major order
#' @keywords internal
.upper_pairs <- function(D) {
  idx <- which(upper.tri(matrix(0, D, D), diag = TRUE), arr.ind = TRUE)
  idx[order(idx[, 2], idx[, 1]), , drop = FALSE]
}

#' Per-pair summary statistics for every i <= j pair, from in-memory draws
#'
#' Returns one row per pair in `full` (in the same order): the
#' `ci_alpha/2`/`1 - ci_alpha/2` quantiles of the lower/upper bound draws,
#' the mean and quantiles of the relative-covariance draws, and counts of
#' draws whose lower bound exceeds `delta` / upper bound is below
#' `-delta` (the ingredients for a p-value; meaningless for the diagonal
#' but cheap to compute for every pair uniformly).
#' @keywords internal
.memory_pair_stats <- function(draws, full, ci_alpha, delta) {
  i <- full[, 1]
  j <- full[, 2]
  lower_draws <- .pair_draws(draws$lower, i, j)
  upper_draws <- .pair_draws(draws$upper, i, j)
  rel_draws <- .pair_draws(draws$rel, i, j)
  list(
    ci_lower = apply(lower_draws, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE),
    ci_upper = apply(upper_draws, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE),
    rel_hat = rowMeans(rel_draws),
    rel_ci_lower = apply(rel_draws, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE),
    rel_ci_upper = apply(rel_draws, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE),
    pos_count = rowSums(lower_draws > delta),
    neg_count = rowSums(upper_draws < -delta)
  )
}

#' Per-pair summary statistics for every i <= j pair, from streamed draws
#'
#' Same contract as `.memory_pair_stats()`, but reads pair-blocks from the
#' temporary binary files written by `.prism_bootstrap()` instead of
#' indexing an in-memory array, bounding peak memory to
#' `block_size x S_eff` doubles regardless of `nrow(full)`.
#' @keywords internal
.stream_pair_stats <- function(draws, full, ci_alpha, delta, block_size = 20000L) {
  n_pairs <- nrow(full)
  S <- draws$S_eff
  out <- lapply(
    c("ci_lower", "ci_upper", "rel_hat", "rel_ci_lower", "rel_ci_upper", "pos_count", "neg_count"),
    function(nm) numeric(n_pairs)
  )
  names(out) <- c("ci_lower", "ci_upper", "rel_hat", "rel_ci_lower", "rel_ci_upper", "pos_count", "neg_count")

  starts <- seq.int(1L, n_pairs, by = block_size)
  for (start in starts) {
    end <- min(start + block_size - 1L, n_pairs)
    idx <- start:end
    lower_block <- .stream_read_block(draws$files$lower, n_pairs, S, start, end)
    upper_block <- .stream_read_block(draws$files$upper, n_pairs, S, start, end)
    rel_block <- .stream_read_block(draws$files$rel, n_pairs, S, start, end)
    out$ci_lower[idx] <- apply(lower_block, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE)
    out$ci_upper[idx] <- apply(upper_block, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)
    out$rel_hat[idx] <- rowMeans(rel_block)
    out$rel_ci_lower[idx] <- apply(rel_block, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE)
    out$rel_ci_upper[idx] <- apply(rel_block, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)
    out$pos_count[idx] <- rowSums(lower_block > delta)
    out$neg_count[idx] <- rowSums(upper_block < -delta)
  }
  out
}

#' Assemble CI matrices and pairwise tables from per-pair summary statistics
#' @keywords internal
.assemble_prism_aggregate <- function(full, stats_full, feature_names, D, S_eff, p_adjust_method) {
  fi <- full[, 1]
  fj <- full[, 2]

  ci_lower_mat <- matrix(0, D, D)
  ci_upper_mat <- matrix(0, D, D)
  ci_lower_mat[cbind(fi, fj)] <- stats_full$ci_lower
  ci_lower_mat[cbind(fj, fi)] <- stats_full$ci_lower
  ci_upper_mat[cbind(fi, fj)] <- stats_full$ci_upper
  ci_upper_mat[cbind(fj, fi)] <- stats_full$ci_upper
  dimnames(ci_lower_mat) <- dimnames(ci_upper_mat) <- list(feature_names, feature_names)

  is_offdiag <- fi != fj
  i <- fi[is_offdiag]
  j <- fj[is_offdiag]
  if (S_eff == 1L) {
    p_value <- rep(NA_real_, length(i))
  } else {
    pos_count <- stats_full$pos_count[is_offdiag]
    neg_count <- stats_full$neg_count[is_offdiag]
    p_plus <- (1 + pos_count) / (1 + S_eff)
    p_minus <- (1 + neg_count) / (1 + S_eff)
    p_value <- pmin(1, 2 * pmin(p_plus, p_minus))
  }
  q_value <- stats::p.adjust(p_value, method = p_adjust_method)

  pairwise <- data.frame(
    i = i, j = j,
    feature_i = feature_names[i], feature_j = feature_names[j],
    ci_lower = stats_full$ci_lower[is_offdiag], ci_upper = stats_full$ci_upper[is_offdiag],
    p_value = p_value, q_value = q_value,
    stringsAsFactors = FALSE
  )

  pairwise_rel <- data.frame(
    i = fi, j = fj,
    feature_i = feature_names[fi], feature_j = feature_names[fj],
    rel_hat = stats_full$rel_hat,
    rel_ci_lower = stats_full$rel_ci_lower,
    rel_ci_upper = stats_full$rel_ci_upper,
    stringsAsFactors = FALSE
  )
  pairwise_rel$rel_ci_width <- pairwise_rel$rel_ci_upper - pairwise_rel$rel_ci_lower

  list(ci_lower = ci_lower_mat, ci_upper = ci_upper_mat, pairwise = pairwise, pairwise_rel = pairwise_rel)
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
#' @param draws Output of `.prism_bootstrap()`.
#' @param feature_names Character vector of length D, or `NULL`.
#' @param ci_alpha Tail probability for CI endpoints.
#' @param delta Half-width of the pairwise null region used for p-values.
#' @param p_adjust_method Passed to `stats::p.adjust()`.
#' @return `list(ci_lower, ci_upper, pairwise, pairwise_rel)`.
#' @keywords internal
.prism_aggregate <- function(draws, feature_names = NULL, ci_alpha = 0.05,
                              delta = 0, p_adjust_method = "BH") {
  D <- draws$D
  S_eff <- draws$S_eff
  ci_alpha <- .assert_scalar_finite(ci_alpha, "ci_alpha", lower = 0, upper = 1, upper_inclusive = FALSE)
  delta <- .assert_scalar_finite(delta, "delta", lower = 0)
  if (!p_adjust_method %in% stats::p.adjust.methods) {
    stop("p_adjust_method must be one of: ", paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)
  }
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(D))

  full <- .upper_pairs(D)
  stats_full <- if (identical(draws$storage, "stream")) {
    .stream_pair_stats(draws, full, ci_alpha, delta)
  } else {
    .memory_pair_stats(draws, full, ci_alpha, delta)
  }
  .assemble_prism_aggregate(full, stats_full, feature_names, D, S_eff, p_adjust_method)
}
