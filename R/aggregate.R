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
#' quantiles of the relative-covariance draws, bootstrap sign counts, and
#' counts used for a two-sided bootstrap p-value against `[-delta, delta]`.
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
    rel_ci_lower = apply(rel_draws, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE),
    rel_ci_upper = apply(rel_draws, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE),
    positive_count = rowSums(lower_draws > 0),
    negative_count = rowSums(upper_draws < 0),
    covers_zero_count = rowSums(lower_draws <= 0 & upper_draws >= 0),
    positive_null_count = rowSums(lower_draws > delta),
    negative_null_count = rowSums(upper_draws < -delta)
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
    c(
      "ci_lower", "ci_upper", "rel_ci_lower", "rel_ci_upper",
      "positive_count", "negative_count", "covers_zero_count",
      "positive_null_count", "negative_null_count"
    ),
    function(nm) numeric(n_pairs)
  )
  names(out) <- c(
    "ci_lower", "ci_upper", "rel_ci_lower", "rel_ci_upper",
    "positive_count", "negative_count", "covers_zero_count",
    "positive_null_count", "negative_null_count"
  )

  starts <- seq.int(1L, n_pairs, by = block_size)
  for (start in starts) {
    end <- min(start + block_size - 1L, n_pairs)
    idx <- start:end
    lower_block <- .stream_read_block(draws$files$lower, n_pairs, S, start, end)
    upper_block <- .stream_read_block(draws$files$upper, n_pairs, S, start, end)
    rel_block <- .stream_read_block(draws$files$rel, n_pairs, S, start, end)
    out$ci_lower[idx] <- apply(lower_block, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE)
    out$ci_upper[idx] <- apply(upper_block, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)
    out$rel_ci_lower[idx] <- apply(rel_block, 1, stats::quantile, probs = ci_alpha / 2, names = FALSE)
    out$rel_ci_upper[idx] <- apply(rel_block, 1, stats::quantile, probs = 1 - ci_alpha / 2, names = FALSE)
    out$positive_count[idx] <- rowSums(lower_block > 0)
    out$negative_count[idx] <- rowSums(upper_block < 0)
    out$covers_zero_count[idx] <- rowSums(lower_block <= 0 & upper_block >= 0)
    out$positive_null_count[idx] <- rowSums(lower_block > delta)
    out$negative_null_count[idx] <- rowSums(upper_block < -delta)
  }
  out
}

#' Assemble CI matrices and pairwise tables from per-pair summary statistics
#' @keywords internal
.assemble_prism_aggregate <- function(full, stats_full, feature_names, D, S_eff,
                                      bootstrap_active, p_adjust_method) {
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
  positive_count <- stats_full$positive_count[is_offdiag]
  negative_count <- stats_full$negative_count[is_offdiag]
  covers_zero_count <- stats_full$covers_zero_count[is_offdiag]
  if (!bootstrap_active) {
    p_value <- rep(NA_real_, length(i))
    bootstrap_positive_count <- bootstrap_negative_count <- bootstrap_covers_zero_count <- rep(NA_integer_, length(i))
    bootstrap_positive_frequency <- bootstrap_negative_frequency <- bootstrap_covers_zero_frequency <- rep(NA_real_, length(i))
  } else {
    positive_null_count <- stats_full$positive_null_count[is_offdiag]
    negative_null_count <- stats_full$negative_null_count[is_offdiag]
    p_positive <- (1 + S_eff - positive_null_count) / (1 + S_eff)
    p_negative <- (1 + S_eff - negative_null_count) / (1 + S_eff)
    p_value <- pmin(1, 2 * pmin(p_positive, p_negative))
    bootstrap_positive_count <- as.integer(positive_count)
    bootstrap_negative_count <- as.integer(negative_count)
    bootstrap_covers_zero_count <- as.integer(covers_zero_count)
    bootstrap_positive_frequency <- positive_count / S_eff
    bootstrap_negative_frequency <- negative_count / S_eff
    bootstrap_covers_zero_frequency <- covers_zero_count / S_eff
  }
  q_value <- stats::p.adjust(p_value, method = p_adjust_method)

  pairwise <- data.frame(
    i = i, j = j,
    feature_i = feature_names[i], feature_j = feature_names[j],
    ci_lower = stats_full$ci_lower[is_offdiag], ci_upper = stats_full$ci_upper[is_offdiag],
    covers_zero = stats_full$ci_lower[is_offdiag] <= 0 & stats_full$ci_upper[is_offdiag] >= 0,
    bootstrap_positive_count = bootstrap_positive_count,
    bootstrap_negative_count = bootstrap_negative_count,
    bootstrap_covers_zero_count = bootstrap_covers_zero_count,
    bootstrap_positive_frequency = bootstrap_positive_frequency,
    bootstrap_negative_frequency = bootstrap_negative_frequency,
    bootstrap_covers_zero_frequency = bootstrap_covers_zero_frequency,
    p_value = p_value, q_value = q_value,
    stringsAsFactors = FALSE
  )

  pairwise_rel <- data.frame(
    i = fi, j = fj,
    feature_i = feature_names[fi], feature_j = feature_names[fj],
    rel_ci_lower = stats_full$rel_ci_lower,
    rel_ci_upper = stats_full$rel_ci_upper,
    stringsAsFactors = FALSE
  )
  pairwise_rel$rel_ci_width <- pairwise_rel$rel_ci_upper - pairwise_rel$rel_ci_lower

  list(ci_lower = ci_lower_mat, ci_upper = ci_upper_mat, pairwise = pairwise, pairwise_rel = pairwise_rel)
}

#' Aggregate bootstrap draws into CI matrices and pairwise tables
#'
#' `ci_lower`/`ci_upper` are the empirical `ci_alpha/2` and
#' `1 - ci_alpha/2` quantiles of the per-draw lower/upper bound
#' arrays -- for a single deterministic draw (`S_eff == 1`), the quantile
#' of one value is that value, so this reduces automatically to the
#' identification interval with no special case needed. `p_value` is the
#' reported only when paired subject bootstrapping is active. Without subject
#' bootstrapping, `p_value` and `q_value` are `NA`, while `covers_zero` reports
#' whether the aggregated identification interval includes zero.
#'
#' @param draws Output of `.prism_bootstrap()`.
#' @param feature_names Character vector of length D, or `NULL`.
#' @param ci_alpha Tail probability for CI endpoints.
#' @param delta Half-width of the pairwise null region used for p-values.
#' @param bootstrap_active Whether paired subject bootstrapping was performed.
#' @param p_adjust_method Passed to `stats::p.adjust()`.
#' @return `list(ci_lower, ci_upper, pairwise, pairwise_rel)`.
#' @keywords internal
.prism_aggregate <- function(draws, feature_names = NULL, ci_alpha = 0.05,
                              delta = 0, bootstrap_active = TRUE,
                              p_adjust_method = "BH") {
  D <- draws$D
  S_eff <- draws$S_eff
  ci_alpha <- .assert_scalar_finite(ci_alpha, "ci_alpha", lower = 0, upper = 1, upper_inclusive = FALSE)
  delta <- .assert_scalar_finite(delta, "delta", lower = 0)
  if (!is.logical(bootstrap_active) || length(bootstrap_active) != 1L || is.na(bootstrap_active)) {
    stop("bootstrap_active must be TRUE or FALSE.", call. = FALSE)
  }
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
  .assemble_prism_aggregate(
    full, stats_full, feature_names, D, S_eff,
    bootstrap_active, p_adjust_method
  )
}
