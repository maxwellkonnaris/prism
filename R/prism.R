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
#' @details
#' `scale` is the single argument for "what do I know about total scale?"
#' and accepts, in increasing order of how much it specifies:
#'
#' \describe{
#'   \item{`NULL` (default)}{no scale information (unbounded-scale regime),
#'     unless `sigma_L`/`sigma_U`/`rho_L`/`rho_U` are also given.}
#'   \item{a bare numeric vector (length N) or matrix (N rows, one column
#'     per replicate measurement)}{sample-level scale data; sigma/rho are
#'     estimated within each draw. Values are assumed raw (un-logged) and
#'     PRISM takes the log itself, *unless* any value is non-positive (a
#'     raw scale measurement cannot be), in which case the whole input is
#'     treated as already log-scale. All-positive values that are already
#'     log-scale are ambiguous by value alone and must be wrapped
#'     explicitly with [prism_scale_log()].}
#'   \item{[prism_scale_log()]}{sample-level scale data that is explicitly
#'     already on the log scale, with its own `ci_level`/`estimate_rho`/
#'     `lower_zero` settings.}
#'   \item{[prism_scale_bounds()]}{fixed `sigma_L`/`sigma_U`/`rho_L`/`rho_U`
#'     bounds, no per-sample data. `sigma_L`/`sigma_U`/`rho_L`/`rho_U` can
#'     be passed directly to `prism()` as a shorthand for this -- they are
#'     forwarded into a `prism_scale_bounds()` for you.}
#' }
#'
#' `scale_ci_level`/`scale_estimate_rho`/`scale_lower_zero` configure the
#' within-draw scale-SD/correlation estimate when `scale` is sample-level
#' data (a bare vector/matrix); they cannot be combined with
#' `scale = prism_scale_log(...)` (which carries its own settings) or with
#' `scale = prism_scale_bounds(...)`/`sigma_L` etc. (fixed, nothing to
#' estimate).
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
#' @param scale Scale information; see Details.
#' @param sigma_L,sigma_U Optional fixed scale-SD bounds,
#'   `0 <= sigma_L <= sigma_U`. Shorthand for
#'   `scale = prism_scale_bounds(sigma_L, sigma_U, rho_L, rho_U)`; cannot
#'   be combined with a non-`NULL` `scale`.
#' @param rho_L,rho_U Optional fixed scalar or feature-level correlation
#'   bounds in `[-1, 1]`; require `sigma_L`/`sigma_U`.
#' @param scale_ci_level Coverage for the chi-square scale-SD and
#'   Fisher-z correlation intervals estimated within each draw, when
#'   `scale` is a bare vector/matrix. `NULL` (default) uses `0` (point
#'   estimates).
#' @param scale_estimate_rho Estimate taxon-scale correlation, when
#'   `scale` is a bare vector/matrix. `NULL` (default) uses `TRUE`. If
#'   `FALSE`, only scale-SD is estimated (the bound regime becomes
#'   bounded-scale-only).
#' @param scale_lower_zero Force the lower scale-SD endpoint to `0`
#'   instead of the chi-square lower quantile, when `scale` is a bare
#'   vector/matrix (sigma cannot be negative, so the two-sided lower
#'   quantile asserts a positive floor that is an artifact of the
#'   two-sided interval, not a real constraint). `NULL` (default) uses
#'   `TRUE`.
#' @param ci_alpha Tail probability for the final confidence intervals.
#' @param delta Half-width of the pairwise null region used for p-values.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param seed Reproducibility seed, or `NULL`.
#' @param stream Stream draws to temporary disk files instead of holding
#'   them in memory. `FALSE` (the default) still streams automatically
#'   when the estimated cost of holding all draws in memory exceeds
#'   `stream_memory_limit`, so a large problem does not silently exhaust
#'   memory.
#' @param stream_memory_limit Byte threshold for the automatic streaming
#'   trigger; see `stream`. Default 2e9 (2 GB).
#'
#' @return A `prism_result`: `ci_lower`/`ci_upper` (D x D matrices),
#'   `pairwise` (off-diagonal tested pairs with CI/p-value/q-value),
#'   `pairwise_rel` (relative log-covariance estimates and CIs, including
#'   the diagonal), and `diagnostics`/`parameters`.
#' @export
prism <- function(counts,
                   composition = prism_composition_dirichlet(),
                   bootstrap = TRUE,
                   S = 2000L,
                   max_attempts_per_draw = 100L,
                   prevalence = 0,
                   verbose = FALSE,
                   scale = NULL,
                   sigma_L = NULL,
                   sigma_U = NULL,
                   rho_L = NULL,
                   rho_U = NULL,
                   scale_ci_level = NULL,
                   scale_estimate_rho = NULL,
                   scale_lower_zero = NULL,
                   ci_alpha = 0.05,
                   delta = 0,
                   p_adjust_method = "BH",
                   seed = 1L,
                   stream = FALSE,
                   stream_memory_limit = 2e9) {
  if (!inherits(composition, "prism_composition_estimator")) {
    stop(
      "composition must be a prism_composition_estimator; see ",
      "prism_composition_dirichlet(), prism_composition_fixed(), or ",
      "prism_composition_estimator().",
      call. = FALSE
    )
  }

  filtered <- .prism_filter_data(counts, scale, prevalence, verbose)
  counts <- filtered$counts
  scale <- filtered$scale
  D <- nrow(counts)
  N <- ncol(counts)
  feature_names <- rownames(counts) %||% paste0("V", seq_len(D))

  resolved <- .prism_resolve_scale(
    scale, N, sigma_L, sigma_U, rho_L, rho_U,
    scale_ci_level, scale_estimate_rho, scale_lower_zero
  )

  draws <- .prism_bootstrap(
    counts, composition, resolved$scale_mode,
    sigma_L = resolved$sigma_L, sigma_U = resolved$sigma_U, rho_L = resolved$rho_L, rho_U = resolved$rho_U,
    bootstrap = bootstrap, S = S, max_attempts_per_draw = max_attempts_per_draw,
    seed = seed, verbose = verbose, stream = stream, stream_memory_limit = stream_memory_limit
  )
  if (identical(draws$storage, "stream")) {
    on.exit(unlink(unlist(draws$files)), add = TRUE)
  }
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
      S_effective = draws$S_eff,
      ci_alpha = ci_alpha,
      delta = delta,
      p_adjust_method = p_adjust_method,
      seed = seed,
      scale_mode = resolved$scale_mode$mode,
      stream_active = draws$diagnostics$stream_active
    )
  ))
}

#' Resolve the `scale` argument (plus its shorthands) into a bootstrap-ready form
#'
#' @return `list(scale_mode, sigma_L, sigma_U, rho_L, rho_U)`, where
#'   `scale_mode` is ready to pass to `.prism_bootstrap()` and the
#'   `sigma_L`/etc. are the fixed-bounds values (or all `NULL` when
#'   `scale_mode$mode == "data"`).
#' @keywords internal
.prism_resolve_scale <- function(scale, N, sigma_L, sigma_U, rho_L, rho_U,
                                  scale_ci_level, scale_estimate_rho, scale_lower_zero) {
  flat_bounds <- !is.null(sigma_L) || !is.null(sigma_U) || !is.null(rho_L) || !is.null(rho_U)
  flat_settings <- !is.null(scale_ci_level) || !is.null(scale_estimate_rho) || !is.null(scale_lower_zero)
  if (!is.null(scale_ci_level)) {
    scale_ci_level <- .assert_scalar_finite(scale_ci_level, "scale_ci_level", lower = 0, upper = 1, upper_inclusive = FALSE)
  }
  if (!is.null(scale_estimate_rho) && (!is.logical(scale_estimate_rho) || length(scale_estimate_rho) != 1L || is.na(scale_estimate_rho))) {
    stop("scale_estimate_rho must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(scale_lower_zero) && (!is.logical(scale_lower_zero) || length(scale_lower_zero) != 1L || is.na(scale_lower_zero))) {
    stop("scale_lower_zero must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(scale)) {
    if (flat_settings) {
      stop(
        "scale_ci_level/scale_estimate_rho/scale_lower_zero require scale ",
        "to be sample-level data (a vector, matrix, or prism_scale_log()).",
        call. = FALSE
      )
    }
    if (!flat_bounds) {
      return(list(scale_mode = list(mode = "none"), sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL))
    }
    scale <- prism_scale_bounds(sigma_L, sigma_U, rho_L, rho_U)
  } else if (flat_bounds) {
    stop(
      "sigma_L/sigma_U/rho_L/rho_U cannot be combined with a non-NULL ",
      "scale argument; use scale = prism_scale_bounds(...) instead.",
      call. = FALSE
    )
  }

  if (inherits(scale, "prism_scale_bounds")) {
    if (flat_settings) {
      stop(
        "scale_ci_level/scale_estimate_rho/scale_lower_zero do not apply ",
        "to scale = prism_scale_bounds(...) (there is no data to estimate a CI from).",
        call. = FALSE
      )
    }
    return(list(
      scale_mode = list(mode = "none"),
      sigma_L = scale$sigma_L, sigma_U = scale$sigma_U, rho_L = scale$rho_L, rho_U = scale$rho_U
    ))
  }

  if (inherits(scale, "prism_scale_log")) {
    if (flat_settings) {
      stop(
        "scale_ci_level/scale_estimate_rho/scale_lower_zero are already ",
        "set on scale = prism_scale_log(...); do not also pass them to prism().",
        call. = FALSE
      )
    }
    mat <- .validate_scale_input(scale$values, N, require_positive = FALSE)
    return(list(
      scale_mode = list(
        mode = "data", scale_mat = mat, log_transform = FALSE,
        ci_level = scale$ci_level, estimate_rho = scale$estimate_rho, lower_zero = scale$lower_zero
      ),
      sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL
    ))
  }

  # A bare vector or matrix: raw sample-level scale, unless any value is
  # non-positive (impossible for a raw measurement), in which case it must
  # already be log-scale.
  raw_values <- if (is.null(dim(scale))) scale else as.matrix(scale)
  is_log <- is.numeric(raw_values) && any(raw_values <= 0, na.rm = TRUE)
  mat <- .validate_scale_input(scale, N, require_positive = !is_log)
  list(
    scale_mode = list(
      mode = "data", scale_mat = mat, log_transform = !is_log,
      ci_level = scale_ci_level %||% 0,
      estimate_rho = scale_estimate_rho %||% TRUE,
      lower_zero = scale_lower_zero %||% TRUE
    ),
    sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL
  )
}

#' Filter counts (and paired scale data) to a valid PRISM input
#' @keywords internal
.prism_filter_data <- function(counts, scale, prevalence, verbose) {
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

  scale_values <- if (is.null(scale) || inherits(scale, "prism_scale_bounds")) {
    NULL
  } else if (inherits(scale, "prism_scale_log")) {
    scale$values
  } else {
    scale
  }
  if (!is.null(scale_values)) {
    n_scale <- if (is.null(dim(scale_values))) length(scale_values) else nrow(scale_values)
    if (n_scale != N_before) {
      stop("scale must have length N (one value per sample) or N rows.", call. = FALSE)
    }
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
  if (!is.null(scale_values)) {
    scale_values <- if (is.null(dim(scale_values))) scale_values[keep_samples] else scale_values[keep_samples, , drop = FALSE]
  }

  feature_totals <- rowSums(Y)
  keep_features2 <- feature_totals > 0
  if (sum(keep_features2) < 2L) {
    stop("Fewer than 2 features remain after dropping zero-sum features.", call. = FALSE)
  }
  Y <- Y[keep_features2, , drop = FALSE]

  scale_filtered <- if (is.null(scale) || inherits(scale, "prism_scale_bounds")) {
    scale
  } else if (inherits(scale, "prism_scale_log")) {
    scale$values <- scale_values
    scale
  } else {
    scale_values
  }

  list(counts = Y, scale = scale_filtered, D_before = D_before, N_before = N_before)
}
