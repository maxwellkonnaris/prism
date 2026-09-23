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
#'   \item{numeric vector or matrix}{sample-level log-scale data. A matrix
#'     contains technical replicates in columns.}
#'   \item{[prism_scale_log()]}{sample-level log-scale data with non-default
#'     within-draw estimation settings.}
#'   \item{[prism_scale_bounds()]}{fixed `sigma_L`/`sigma_U`/`rho_L`/`rho_U`
#'     bounds, no per-sample data. `sigma_L`/`sigma_U`/`rho_L`/`rho_U` can
#'     be passed directly to `prism()` as a shorthand for this -- they are
#'     forwarded into a `prism_scale_bounds()` for you.}
#' }
#'
#' All sample-level scale data are assumed to be already log transformed;
#' PRISM never applies a logarithm. Bare numeric input uses point estimates
#' within each bootstrap draw. [prism_scale_log()] changes those settings;
#' analytic chi-square/Fisher-z intervals inside draws are available only as
#' an explicitly enabled experimental sensitivity analysis.
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
#' @param prevalence Minimum feature prevalence retained, in `[0, 1]`.
#' @param verbose Emit progress messages.
#' @param scale Scale information; see Details.
#' @param sigma_L,sigma_U Optional fixed scale-SD bounds,
#'   `0 <= sigma_L <= sigma_U`. Shorthand for
#'   `scale = prism_scale_bounds(sigma_L, sigma_U, rho_L, rho_U)`; cannot
#'   be combined with a non-`NULL` `scale`.
#' @param rho_L,rho_U Optional fixed scalar or feature-level correlation
#'   bounds in `[-1, 1]`; require `sigma_L`/`sigma_U`.
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
#' @param return_draws Return every draw's lower/upper/relative
#'   -covariance values, not just their aggregated summary. Adds `draws`
#'   to the result: `list(lower, upper, rel)`, each an S_effective x
#'   n_pairs matrix, plus `pair_index` mapping columns to `(i, j)` and
#'   feature names (the same pairs as `pairwise_rel`, i.e. every `i <= j`
#'   including the diagonal). For example, the endpoints in `ci_lower`/
#'   `ci_upper` are exactly `apply(draws$lower, 2, quantile, probs =
#'   ci_alpha / 2)` and the equivalent for `draws$upper`. This can be
#'   large: two S_effective x n_pairs matrices plus one more, in double
#'   precision.
#' @param workers Positive integer number of bootstrap draws evaluated
#'   concurrently. `1` is sequential. Parallel draws are collected in
#'   adaptive memory-bounded batches; only the parent process writes streamed
#'   output.
#' @param rho_backend Rho-support optimization backend: direct `"ecos"`
#'   (default) or the retained `"cvxr"` reference implementation.
#' @param solver Optional conic solver name. Direct `rho_backend = "ecos"`
#'   accepts `NULL` or `"ECOS"`; the CVXR reference backend accepts any
#'   installed CVXR conic solver.
#'
#' @return A `prism_result`: `ci_lower`/`ci_upper` (D x D matrices),
#'   `pairwise` (off-diagonal tested pairs with CI/p-value/q-value),
#'   `pairwise_rel` (relative log-covariance CIs, including
#'   the diagonal), `diagnostics`/`parameters`, and (if `return_draws`)
#'   `draws`.
#' @export
prism <- function(counts,
                   composition = prism_composition_dirichlet(),
                   bootstrap = TRUE,
                   S = 2000L,
                   prevalence = 0,
                   verbose = FALSE,
                   scale = NULL,
                   sigma_L = NULL,
                   sigma_U = NULL,
                   rho_L = NULL,
                   rho_U = NULL,
                   ci_alpha = 0.05,
                   delta = 0,
                   p_adjust_method = "BH",
                   seed = 1L,
                   stream = FALSE,
                   stream_memory_limit = 2e9,
                   return_draws = FALSE,
                   workers = 1L,
                   rho_backend = c("ecos", "cvxr"),
                   solver = NULL) {
  if (!inherits(composition, "prism_composition_estimator")) {
    stop(
      "composition must be a prism_composition_estimator; see ",
      "prism_composition_dirichlet(), prism_composition_fixed(), or ",
      "prism_composition_estimator().",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    old_rng <- .save_rng_state()
    on.exit(.restore_rng_state(old_rng), add = TRUE)
  }

  filtered <- .prism_filter_data(counts, scale, prevalence, verbose)
  counts <- filtered$counts
  scale <- filtered$scale
  D <- nrow(counts)
  N <- ncol(counts)
  feature_names <- rownames(counts) %||% paste0("V", seq_len(D))

  resolved <- .prism_resolve_scale(
    scale, N, sigma_L, sigma_U, rho_L, rho_U
  )
  rho_backend <- match.arg(rho_backend)

  draws <- .prism_bootstrap(
    counts, composition, resolved$scale_mode,
    sigma_L = resolved$sigma_L, sigma_U = resolved$sigma_U, rho_L = resolved$rho_L, rho_U = resolved$rho_U,
    bootstrap = bootstrap, S = S,
    seed = seed, verbose = verbose, stream = stream,
    stream_memory_limit = stream_memory_limit,
    workers = workers, rho_backend = rho_backend, solver = solver
  )
  if (identical(draws$storage, "stream")) {
    on.exit(unlink(unlist(draws$files)), add = TRUE)
  }
  bootstrap_mode <- .resolve_bootstrap_mode(bootstrap)
  agg <- .prism_aggregate(
    draws, feature_names, ci_alpha = ci_alpha, delta = delta,
    bootstrap_active = bootstrap_mode == "both",
    p_adjust_method = p_adjust_method
  )
  draws_out <- if (return_draws) .extract_prism_draws(draws, feature_names) else NULL

  .new_prism_result(list(
    ci_lower = agg$ci_lower,
    ci_upper = agg$ci_upper,
    pairwise = agg$pairwise,
    pairwise_rel = agg$pairwise_rel,
    draws = draws_out,
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
      bootstrap = bootstrap_mode,
      S = S,
      S_effective = draws$S_eff,
      ci_alpha = ci_alpha,
      delta = delta,
      p_adjust_method = p_adjust_method,
      seed = seed,
      scale_mode = resolved$scale_mode$mode,
      stream_active = draws$diagnostics$stream_active,
      workers_requested = draws$diagnostics$workers_requested,
      workers_used = draws$diagnostics$workers_used,
      parallel_active = draws$diagnostics$parallel_active,
      rho_backend = rho_backend,
      solver = solver
    )
  ))
}

#' Extract every draw's lower/upper/rel values, pair-major
#'
#' Reads from the in-memory arrays or, for a streamed run, straight from
#' the temporary files (still open at this point -- called before the
#' stream cleanup `on.exit()` fires).
#'
#' @param draws Output of `.prism_bootstrap()`.
#' @param feature_names Character vector of length D.
#' @return `list(lower, upper, rel)`, each an S_eff x n_pairs matrix, and
#'   `pair_index` (a data frame with `i`, `j`, `feature_i`, `feature_j`).
#' @keywords internal
.extract_prism_draws <- function(draws, feature_names) {
  pair_index <- .upper_pairs(draws$D)
  i <- pair_index[, 1]
  j <- pair_index[, 2]
  if (identical(draws$storage, "stream")) {
    n_pairs <- nrow(pair_index)
    lower <- t(.stream_read_block(draws$files$lower, n_pairs, draws$S_eff, 1L, n_pairs))
    upper <- t(.stream_read_block(draws$files$upper, n_pairs, draws$S_eff, 1L, n_pairs))
    rel <- t(.stream_read_block(draws$files$rel, n_pairs, draws$S_eff, 1L, n_pairs))
  } else {
    lower <- t(.pair_draws(draws$lower, i, j))
    upper <- t(.pair_draws(draws$upper, i, j))
    rel <- t(.pair_draws(draws$rel, i, j))
  }
  list(
    lower = lower, upper = upper, rel = rel,
    pair_index = data.frame(
      i = i, j = j, feature_i = feature_names[i], feature_j = feature_names[j],
      stringsAsFactors = FALSE
    )
  )
}

#' Resolve the `scale` argument (plus its shorthands) into a bootstrap-ready form
#'
#' @return `list(scale_mode, sigma_L, sigma_U, rho_L, rho_U)`, where
#'   `scale_mode` is ready to pass to `.prism_bootstrap()` and the
#'   `sigma_L`/etc. are the fixed-bounds values (or all `NULL` when
#'   `scale_mode$mode == "data"`).
#' @keywords internal
.prism_resolve_scale <- function(scale, N, sigma_L, sigma_U, rho_L, rho_U) {
  flat_bounds <- !is.null(sigma_L) || !is.null(sigma_U) || !is.null(rho_L) || !is.null(rho_U)

  if (is.null(scale)) {
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
    return(list(
      scale_mode = list(mode = "none"),
      sigma_L = scale$sigma_L, sigma_U = scale$sigma_U, rho_L = scale$rho_L, rho_U = scale$rho_U
    ))
  }

  if (inherits(scale, "prism_scale_log")) {
    mat <- .validate_scale_input(scale$values, N)
    return(list(
      scale_mode = list(
        mode = "data", scale_mat = mat,
        ci_level = scale$ci_level, estimate_rho = scale$estimate_rho, lower_zero = scale$lower_zero
      ),
      sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL
    ))
  }

  if (is.numeric(scale)) {
    return(list(
      scale_mode = list(
        mode = "data", scale_mat = .validate_scale_input(scale, N),
        ci_level = 0, estimate_rho = TRUE, lower_zero = TRUE
      ),
      sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL
    ))
  }

  stop("scale must be NULL, numeric log-scale data, prism_scale_log(), or prism_scale_bounds().", call. = FALSE)
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
  } else if (is.numeric(scale)) {
    scale
  } else {
    stop("scale must be NULL, numeric log-scale data, prism_scale_log(), or prism_scale_bounds().", call. = FALSE)
  }
  if (!is.null(scale_values)) {
    n_scale <- if (is.null(dim(scale_values))) length(scale_values) else nrow(scale_values)
    if (n_scale != N_before) {
      stop("scale must have length N (one value per sample) or N rows.", call. = FALSE)
    }
  }

  sample_totals <- colSums(Y)
  keep_samples <- is.finite(sample_totals) & sample_totals > 0
  if (any(!keep_samples)) {
    warning("Dropping ", sum(!keep_samples), " zero-depth sample(s).", call. = FALSE)
  }
  if (sum(keep_samples) < 2L) {
    stop("Fewer than 2 samples remain after dropping zero-depth samples.", call. = FALSE)
  }
  Y <- Y[, keep_samples, drop = FALSE]
  if (!is.null(scale_values)) {
    scale_values <- if (is.null(dim(scale_values))) scale_values[keep_samples] else scale_values[keep_samples, , drop = FALSE]
  }

  prevalence_frac <- rowMeans(Y > 0)
  feature_totals <- rowSums(Y)
  keep_features <- prevalence_frac >= prevalence & feature_totals > 0
  if (sum(keep_features) < 2L) {
    stop("Fewer than 2 features remain after prevalence and zero-sum filtering.", call. = FALSE)
  }
  if (any(!keep_features)) {
    warning(
      "Dropping ", sum(!keep_features),
      " feature(s) failing prevalence or positive-total requirements.",
      call. = FALSE
    )
  }
  Y <- Y[keep_features, , drop = FALSE]

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
