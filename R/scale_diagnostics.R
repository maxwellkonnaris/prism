#' Diagnose a log-scale total vector for common failure modes
#'
#' Runs a fixed battery of checks on a `scale_log`-like input before it
#' reaches [.estimate_scale_log_bounds()] / [cov_bounds()]: a log-base
#' consistency check, a bimodality check (unit mixing), robust-vs-classical
#' spread, a physical fold-range check, zero/non-positive handling on the raw
#' scale, and (if `group` is supplied) a one-way within/between variance
#' decomposition. Each check is a diagnostic flag for a human to look at, not
#' a pass/fail gate -- none of them alter `log_values`.
#'
#' @param log_values Numeric vector, already on the log scale (natural log
#'   unless the values are known to be on a different base -- this function
#'   does not assume a base for `log_values` itself, only checks whether one
#'   is plausible).
#' @param raw_values Optional numeric vector of pre-log values, same
#'   length/order as `log_values`. Enables the log-base-consistency and
#'   zero/LOD checks; without it those two are skipped (`NULL` in the
#'   result) and only `log_values`-only checks run.
#' @param group Optional grouping factor/vector (e.g. cohort, study, batch),
#'   same length as `log_values`. Enables the within/between decomposition;
#'   without it that check is skipped (`NULL` in the result).
#' @param max_fold_range Fold-range on the natural scale above which the
#'   range check is flagged (default `1e4`). Biological total-load variation
#'   across subjects within one study is rarely more than 2-3 orders of
#'   magnitude; this default is deliberately generous so it only fires on
#'   genuinely implausible spans.
#' @param bimodality_threshold Sarle's bimodality coefficient above which the
#'   unit-mixing check is flagged (default `0.555`, the standard rule-of-thumb
#'   value for a uniform distribution; higher indicates more bimodal).
#' @param sd_mad_ratio_threshold `sd / (1.4826 * mad)` ratio above which the
#'   robust-vs-classical check is flagged (default `1.3`).
#'
#' @return A `prism_scale_diagnostics` list with one named element per check
#'   (`log_base`, `bimodality`, `robust_sd`, `fold_range`, `zero_lod`,
#'   `within_between`), each a list with at least `flag` (logical) and
#'   `detail` (a one-line character description), plus check-specific
#'   numeric fields. Elements for checks that could not run (missing
#'   `raw_values` or `group`) are `NULL`.
#'
#' @details
#' **Log-base consistency** only detects a *base* mismatch between
#' `log_values` and `log(raw_values)` (e.g. log10 stored where natural log
#' was assumed, as found in a real pipeline where `logCopyNumber -
#' log10(CopyNumber)` was a near-constant offset). It cannot detect
#' `log_values` being a transformation of the *wrong physical quantity*
#' (e.g. `log10(Ct)` mislabeled as `log10(total abundance)` -- Ct is
#' inversely related to abundance, not a rescaling of it); that failure mode
#' has no numeric signature to check against `raw_values` alone and requires
#' reading the column's actual provenance.
#'
#' **Within/between decomposition** uses the classical method-of-moments
#' one-way random-effects ANOVA estimator (Henderson's Method 1), not a
#' mixed-model fit, to avoid adding a package dependency for a diagnostic:
#' `sigma2_between = max(0, (MSB - MSE) / n0)`, `sigma2_within = MSE`, with
#' `n0` the standard unbalanced-design correction factor.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 10, 0.5), rnorm(50, 10, 0.5))
#' x[1] <- 40  # inject one outlier
#' prism_scale_diagnostics(x)
#'
#' @export
prism_scale_diagnostics <- function(log_values, raw_values = NULL, group = NULL,
                                     max_fold_range = 1e4, bimodality_threshold = 0.555,
                                     sd_mad_ratio_threshold = 1.3) {
  log_values <- as.numeric(log_values)
  n <- length(log_values)
  if (sum(is.finite(log_values)) < 2L) {
    stop("log_values must be numeric with at least 2 finite entries.", call. = FALSE)
  }
  if (!is.null(raw_values) && length(raw_values) != n) {
    stop("raw_values must be the same length as log_values.", call. = FALSE)
  }
  if (!is.null(group) && length(group) != n) {
    stop("group must be the same length as log_values.", call. = FALSE)
  }
  max_fold_range <- .assert_scalar_finite(max_fold_range, "max_fold_range", lower = 0, lower_inclusive = FALSE)
  bimodality_threshold <- .assert_scalar_finite(bimodality_threshold, "bimodality_threshold", lower = 0)
  sd_mad_ratio_threshold <- .assert_scalar_finite(sd_mad_ratio_threshold, "sd_mad_ratio_threshold", lower = 0, lower_inclusive = FALSE)

  ok <- is.finite(log_values)
  v <- log_values[ok]

  result <- list(
    log_base = if (is.null(raw_values)) NULL else .diag_log_base(log_values, raw_values),
    bimodality = .diag_bimodality(v, bimodality_threshold),
    robust_sd = .diag_robust_sd(v, sd_mad_ratio_threshold),
    fold_range = .diag_fold_range(v, max_fold_range),
    zero_lod = if (is.null(raw_values)) NULL else .diag_zero_lod(raw_values),
    within_between = if (is.null(group)) NULL else .diag_within_between(log_values, group)
  )
  structure(result, class = "prism_scale_diagnostics")
}

#' @export
print.prism_scale_diagnostics <- function(x, ...) {
  for (nm in names(x)) {
    entry <- x[[nm]]
    if (is.null(entry)) {
      cat(nm, ": skipped (required input not supplied)\n", sep = "")
      next
    }
    cat(nm, ": ", if (isTRUE(entry$flag)) "FLAGGED -- " else "ok -- ", entry$detail, "\n", sep = "")
  }
  invisible(x)
}

#' @keywords internal
.diag_log_base <- function(log_values, raw_values) {
  finite <- is.finite(log_values) & is.finite(raw_values) & raw_values > 0
  if (sum(finite) < 2L) {
    return(list(flag = NA, detail = "fewer than 2 finite, positive (log_values, raw_values) pairs -- cannot check.",
                ratio_median = NA_real_, ratio_sd = NA_real_, implied_base = NA_character_))
  }
  implied_natural_log <- log(raw_values[finite])
  ratio <- log_values[finite] / implied_natural_log
  ratio <- ratio[is.finite(ratio)]
  med <- stats::median(ratio)
  spread <- stats::mad(ratio)
  # If log_values is truly log_b(raw_values), then log_values = log(raw_values)/log(b),
  # so ratio = log_values/log(raw_values) = 1/log(b) -- e.g. 1/log(10) = 0.4343 for log10.
  candidates <- c(natural = 1, log2 = 1 / log(2), log10 = 1 / log(10))
  nearest <- names(candidates)[which.min(abs(candidates - med))]
  is_consistent_base <- spread < 0.05 * abs(med)
  mismatch <- is_consistent_base && abs(med - 1) > 0.05
  detail <- if (!is_consistent_base) {
    sprintf("log_values is not a fixed-base rescaling of log(raw_values) (ratio median %.3f, MAD %.3f across %d pairs) -- log_values likely encodes a different quantity than raw_values, not just a different log base.",
            med, spread, sum(finite))
  } else if (mismatch) {
    sprintf("log_values ~= %s * log(raw_values) (ratio %.4f, closest to '%s') -- log_values looks like it is on a %s scale, not natural log.",
            nearest, med, nearest, nearest)
  } else {
    sprintf("log_values ~= log(raw_values) (ratio %.4f) -- consistent with natural log.", med)
  }
  list(flag = isTRUE(mismatch) || !is_consistent_base, detail = detail,
       ratio_median = med, ratio_sd = spread, implied_base = nearest)
}

#' @keywords internal
.bimodality_coefficient <- function(v) {
  n <- length(v)
  m <- mean(v); s <- stats::sd(v)
  if (!is.finite(s) || s <= 0 || n < 4L) return(NA_real_)
  skew <- mean((v - m)^3) / s^3
  kurt <- mean((v - m)^4) / s^4
  (skew^2 + 1) / (kurt + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

#' @keywords internal
.diag_bimodality <- function(v, threshold) {
  bc <- .bimodality_coefficient(v)
  flag <- is.finite(bc) && bc > threshold
  detail <- if (!is.finite(bc)) {
    "not enough data to compute a bimodality coefficient."
  } else {
    sprintf("Sarle's bimodality coefficient = %.3f (threshold %.3f) -- %s.",
            bc, threshold, if (flag) "consistent with two or more mixed populations/units" else "no strong bimodality signal")
  }
  list(flag = flag, detail = detail, bimodality_coefficient = bc)
}

#' @keywords internal
.diag_robust_sd <- function(v, threshold) {
  sd_classical <- stats::sd(v)
  mad_robust <- 1.4826 * stats::mad(v)
  ratio <- if (mad_robust > 0) sd_classical / mad_robust else NA_real_
  flag <- is.finite(ratio) && ratio > threshold
  detail <- sprintf("sd = %.4g, 1.4826*mad = %.4g, ratio = %.3f (threshold %.2f) -- %s.",
                     sd_classical, mad_robust, ratio, threshold,
                     if (flag) "classical SD is inflated relative to the robust estimate; a small number of outliers likely dominate it" else "classical and robust spread agree")
  list(flag = flag, detail = detail, sd = sd_classical, mad_robust = mad_robust, ratio = ratio)
}

#' @keywords internal
.diag_fold_range <- function(v, max_fold_range) {
  fold_range <- exp(max(v) - min(v))
  flag <- fold_range > max_fold_range
  detail <- sprintf("exp(max - min) = %.4g-fold (threshold %.4g-fold) -- %s.",
                     fold_range, max_fold_range,
                     if (flag) "exceeds a physically plausible range for total biological load; check for control/blank samples or data-entry errors rather than treating this as real variation" else "within a physically plausible range")
  list(flag = flag, detail = detail, fold_range = fold_range)
}

#' @keywords internal
.diag_zero_lod <- function(raw_values) {
  n <- length(raw_values)
  n_nonpositive <- sum(is.finite(raw_values) & raw_values <= 0, na.rm = TRUE)
  n_nonfinite <- sum(!is.finite(raw_values))
  finite_positive <- raw_values[is.finite(raw_values) & raw_values > 0]
  floor_val <- if (length(finite_positive)) min(finite_positive) else NA_real_
  n_at_floor <- if (is.finite(floor_val)) sum(finite_positive == floor_val) else 0L
  floor_repeat_frac <- if (length(finite_positive)) n_at_floor / length(finite_positive) else NA_real_
  flag <- n_nonpositive > 0 || (is.finite(floor_repeat_frac) && floor_repeat_frac > 0.05 && n_at_floor > 1)
  detail <- sprintf(
    "%d/%d non-positive, %d/%d non-finite raw values; %d observation(s) (%.1f%% of finite-positive values) sit exactly at the observed minimum (%.4g) -- %s.",
    n_nonpositive, n, n_nonfinite, n, n_at_floor, 100 * (floor_repeat_frac %||% 0), floor_val,
    if (flag) "check whether zeros/repeats reflect a detection-limit floor substitution rather than genuine measurements" else "no evidence of LOD floor substitution"
  )
  list(flag = flag, detail = detail, n_nonpositive = n_nonpositive, n_nonfinite = n_nonfinite,
       floor_value = floor_val, n_at_floor = n_at_floor, floor_repeat_frac = floor_repeat_frac)
}

#' @keywords internal
.diag_within_between <- function(log_values, group) {
  ok <- is.finite(log_values) & !is.na(group)
  v <- log_values[ok]; g <- droplevels(as.factor(group[ok]))
  k <- nlevels(g)
  n <- length(v)
  if (k < 2L || n <= k) {
    return(list(flag = NA, detail = "need at least 2 groups and more observations than groups -- cannot decompose.",
                sigma_between = NA_real_, sigma_within = NA_real_, icc = NA_real_))
  }
  grand_mean <- mean(v)
  by_group <- split(v, g)
  n_i <- vapply(by_group, length, integer(1))
  means_i <- vapply(by_group, mean, numeric(1))
  sse <- sum(mapply(function(x, m) sum((x - m)^2), by_group, means_i))
  ssb <- sum(n_i * (means_i - grand_mean)^2)
  mse <- sse / (n - k)
  msb <- ssb / (k - 1)
  n0 <- (n - sum(n_i^2) / n) / (k - 1)
  sigma2_between <- max(0, (msb - mse) / n0)
  sigma2_within <- mse
  icc <- sigma2_between / (sigma2_between + sigma2_within)
  detail <- sprintf(
    "%d groups, n=%d: sigma_between = %.4g, sigma_within = %.4g, ICC = %.3f -- %s.",
    k, n, sqrt(sigma2_between), sqrt(sigma2_within), icc,
    if (icc > 0.8) "pooled variance is dominated by between-group differences (e.g. units/protocol), not within-group signal" else "within-group variance is not negligible relative to between-group"
  )
  list(flag = icc > 0.8, detail = detail, sigma_between = sqrt(sigma2_between), sigma_within = sqrt(sigma2_within), icc = icc)
}
