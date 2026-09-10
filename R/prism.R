#' Estimate absolute log-covariance bounds from compositional data
#'
#' @description
#' For features `i, j`, `Cov(log W_i, log W_j) = Sigma_rel[i,j] + sigma^2 +
#' sigma*(s_i*rho_i + s_j*rho_j)`, where `W` is absolute abundance,
#' `Sigma_rel` is the covariance of sample log-compositions, `sigma = SD(log
#' total scale)`, `s_d = sqrt(Sigma_rel[d,d])`, and `rho_d = Cor(log P_d,
#' log total scale)`. `sigma` and `rho` are generally not point-identified
#' from compositional data alone; `prism()` combines a composition
#' uncertainty estimator (default: Dirichlet-multinomial) with whatever
#' scale information is supplied and returns identification bounds,
#' bootstrap confidence intervals, and pairwise inference.
#'
#' @param counts Non-negative, integer-valued features-by-samples count matrix.
#' @param composition A `prism_composition_estimator` (see
#'   [prism_composition_dirichlet()], [prism_composition_fixed()],
#'   [prism_composition_estimator()]).
#' @param bootstrap `TRUE`/`FALSE` or `"both"`/`"none"`. `"both"` resamples
#'   samples (paired across counts and scale) for every draw; `"none"`
#'   holds samples fixed (composition uncertainty, if any, still varies
#'   across draws).
#' @param S Number of draws.
#' @param max_attempts_per_draw Maximum resample proposals per requested
#'   draw before PRISM gives up; see `diagnostics$sampling` for
#'   acceptance/rejection bookkeeping.
#' @param prevalence Minimum feature prevalence retained, in `[0, 1]`.
#' @param verbose Emit progress messages.
#' @param scale Optional positive sample-level scale: a length-N vector
#'   (one raw measurement per sample) or an N-row replicate matrix (one
#'   row per sample, one column per replicate measurement). Mutually
#'   exclusive with `scale_log` and with `sigma_L`/`sigma_U`/`rho_L`/`rho_U`.
#' @param scale_log Optional finite length-N vector of log total scale.
#'   Mutually exclusive with `scale` and with `sigma_L`/`sigma_U`/`rho_L`/`rho_U`.
#' @param scale_log_ci_level Coverage for the chi-square scale-SD and
#'   Fisher-z correlation intervals estimated from `scale`/`scale_log`
#'   within each draw. `0` uses point estimates.
#' @param scale_log_rho Estimate taxon-scale correlation from
#'   `scale`/`scale_log`. If `FALSE`, only scale-SD is estimated (the
#'   bound regime becomes bounded-scale-only).
#' @param scale_log_lower_zero If `TRUE`, force the lower scale-SD
#'   endpoint to `0` instead of the chi-square lower quantile (sigma
#'   cannot be negative, so the two-sided lower quantile asserts a
#'   positive floor that is an artifact of the two-sided interval, not a
#'   real constraint).
#' @param ci_alpha Tail probability for the final confidence intervals.
#' @param sigma_L,sigma_U Optional fixed scale-SD bounds,
#'   `0 <= sigma_L <= sigma_U`. Mutually exclusive with `scale`/`scale_log`.
#' @param rho_L,rho_U Optional fixed scalar or feature-level correlation
#'   bounds in `[-1, 1]`; require `sigma_L`/`sigma_U`. Mutually exclusive
#'   with `scale`/`scale_log`.
#' @param delta Half-width of the pairwise null region used for p-values.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param seed Reproducibility seed, or `NULL`.
#'
#' @return A `prism_result`: `ci_lower`/`ci_upper` (D x D matrices),
#'   `pairwise` (off-diagonal tested pairs with CI/p-value/q-value),
#'   `pairwise_rel` (relative log-covariance estimates and CIs, including
#'   the diagonal), and `diagnostics`/`parameters`.
#' @export
prism <- function(counts,
                   composition = prism_composition_dirichlet(),
                   bootstrap = TRUE,
                   S = 1000L,
                   max_attempts_per_draw = 100L,
                   prevalence = 0,
                   verbose = FALSE,
                   scale = NULL,
                   scale_log = NULL,
                   scale_log_ci_level = 0,
                   scale_log_rho = TRUE,
                   scale_log_lower_zero = FALSE,
                   ci_alpha = 0.05,
                   sigma_L = NULL,
                   sigma_U = NULL,
                   rho_L = NULL,
                   rho_U = NULL,
                   delta = 0.1,
                   p_adjust_method = "BH",
                   seed = 1L) {
  if (!inherits(composition, "prism_composition_estimator")) {
    stop(
      "composition must be a prism_composition_estimator; see ",
      "prism_composition_dirichlet(), prism_composition_fixed(), or ",
      "prism_composition_estimator().",
      call. = FALSE
    )
  }
  n_scale_inputs <- sum(!is.null(scale), !is.null(scale_log))
  if (n_scale_inputs > 1L) {
    stop("Supply at most one of scale or scale_log.", call. = FALSE)
  }
  has_fixed_bounds <- !is.null(sigma_L) || !is.null(sigma_U) || !is.null(rho_L) || !is.null(rho_U)
  if (n_scale_inputs == 1L && has_fixed_bounds) {
    stop("sigma_L/sigma_U/rho_L/rho_U cannot be combined with scale or scale_log.", call. = FALSE)
  }
  scale_log_rho <- isTRUE(scale_log_rho)
  scale_log_lower_zero <- isTRUE(scale_log_lower_zero)

  filtered <- .prism_filter_data(counts, scale, scale_log, prevalence, verbose)
  counts <- filtered$counts
  scale <- filtered$scale
  scale_log <- filtered$scale_log
  D <- nrow(counts)
  N <- ncol(counts)
  feature_names <- rownames(counts) %||% paste0("V", seq_len(D))

  scale_mode <- if (!is.null(scale)) {
    list(
      mode = "scale", scale_mat = .validate_scale_input(scale, N),
      ci_level = scale_log_ci_level, estimate_rho = scale_log_rho, lower_zero = scale_log_lower_zero
    )
  } else if (!is.null(scale_log)) {
    if (any(!is.finite(scale_log))) {
      stop("scale_log must be a finite length-N vector.", call. = FALSE)
    }
    list(
      mode = "scale_log", scale_log = as.numeric(scale_log),
      ci_level = scale_log_ci_level, estimate_rho = scale_log_rho, lower_zero = scale_log_lower_zero
    )
  } else {
    list(mode = "none")
  }

  draws <- .prism_bootstrap(
    counts, composition, scale_mode,
    sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U,
    bootstrap = bootstrap, S = S, max_attempts_per_draw = max_attempts_per_draw,
    seed = seed, verbose = verbose
  )
  agg <- .prism_aggregate(draws, feature_names, ci_alpha = ci_alpha, delta = delta, p_adjust_method = p_adjust_method)

  .new_prism_result(list(
    ci_lower = agg$ci_lower,
    ci_upper = agg$ci_upper,
    pairwise = agg$pairwise,
    pairwise_rel = agg$pairwise_rel,
    diagnostics = list(
      data = list(
        D = D, N = N,
        D_before = filtered$D_before, N_before = filtered$N_before,
        dropped_samples = filtered$N_before - N,
        dropped_features = filtered$D_before - D
      ),
      sampling = draws$diagnostics
    ),
    parameters = list(
      composition_estimator = composition$name,
      bootstrap = .resolve_bootstrap_mode(bootstrap),
      S = S,
      S_effective = dim(draws$lower)[3],
      ci_alpha = ci_alpha,
      delta = delta,
      p_adjust_method = p_adjust_method,
      seed = seed,
      scale_mode = scale_mode$mode
    )
  ))
}

#' Filter counts (and paired scale data) to a valid PRISM input
#' @keywords internal
.prism_filter_data <- function(counts, scale, scale_log, prevalence, verbose) {
  Y <- as.matrix(counts)
  if (!is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 2L) {
    stop("counts must be a numeric matrix with at least 2 features and 2 samples.", call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(Y < 0)) {
    stop("counts must be finite and non-negative.", call. = FALSE)
  }
  if (any(abs(Y - round(Y)) > 1e-8)) {
    stop("counts must be integer-valued.", call. = FALSE)
  }
  prevalence <- .assert_scalar_finite(prevalence, "prevalence", lower = 0, upper = 1)

  D_before <- nrow(Y)
  N_before <- ncol(Y)
  if (!is.null(scale) && (if (is.null(dim(scale))) length(scale) else nrow(scale)) != N_before) {
    stop("scale must have length N (one value per sample) or N rows.", call. = FALSE)
  }
  if (!is.null(scale_log) && length(scale_log) != N_before) {
    stop("scale_log must have length N (one value per sample).", call. = FALSE)
  }

  prevalence_frac <- rowMeans(Y > 0)
  keep_features <- prevalence_frac >= prevalence
  if (sum(keep_features) < 2L) {
    stop("Fewer than 2 features meet the prevalence threshold.", call. = FALSE)
  }
  Y <- Y[keep_features, , drop = FALSE]

  sample_totals <- colSums(Y)
  keep_samples <- is.finite(sample_totals) & sample_totals > 0
  if (verbose && any(!keep_samples)) {
    message("Dropping ", sum(!keep_samples), " zero-depth sample(s).")
  }
  if (sum(keep_samples) < 2L) {
    stop("Fewer than 2 samples remain after dropping zero-depth samples.", call. = FALSE)
  }
  Y <- Y[, keep_samples, drop = FALSE]
  if (!is.null(scale)) {
    scale <- if (is.null(dim(scale))) scale[keep_samples] else scale[keep_samples, , drop = FALSE]
  }
  if (!is.null(scale_log)) scale_log <- scale_log[keep_samples]

  feature_totals <- rowSums(Y)
  keep_features2 <- feature_totals > 0
  if (sum(keep_features2) < 2L) {
    stop("Fewer than 2 features remain after dropping zero-sum features.", call. = FALSE)
  }
  Y <- Y[keep_features2, , drop = FALSE]

  list(counts = Y, scale = scale, scale_log = scale_log, D_before = D_before, N_before = N_before)
}
