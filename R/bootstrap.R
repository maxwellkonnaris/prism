#' Normalize the bootstrap mode argument
#' @keywords internal
.resolve_bootstrap_mode <- function(bootstrap) {
  if (is.logical(bootstrap) && length(bootstrap) == 1L && !is.na(bootstrap)) {
    return(if (bootstrap) "both" else "none")
  }
  if (is.character(bootstrap) && length(bootstrap) == 1L && bootstrap %in% c("none", "both")) {
    return(bootstrap)
  }
  stop(
    "bootstrap must be TRUE, FALSE, \"none\", or \"both\". Independently ",
    "resampling only the composition or only the scale block is not ",
    "supported: those blocks need not jointly define a positive-",
    "semidefinite covariance unless resampled together.",
    call. = FALSE
  )
}

#' Is this cov_bounds()/relative-covariance error a retriable degenerate draw?
#' @keywords internal
.is_retriable_draw_error <- function(message) {
  patterns <- c("positive semidefinite", "Empty E intersection C", "zero-variance feature")
  any(vapply(patterns, function(p) grepl(p, message, fixed = TRUE), logical(1)))
}

#' Draw one composition realization
#'
#' Dispatches on estimator type. `fit_state` is `NULL` except for a custom
#' estimator, where it is the plug-in's fitted state (built once, outside
#' the draw loop).
#' @keywords internal
.draw_one_composition <- function(estimator, counts, fit_state, draw_index) {
  x <- switch(
    estimator$type,
    fixed = prism_closure(counts, estimator$options$pseudocount),
    dirichlet_multinomial = {
      shape <- estimator$options$concentration * (counts + estimator$options$pseudocount)
      matrix(stats::rgamma(length(shape), shape = shape, rate = 1), nrow = nrow(shape))
    },
    custom = estimator$options$draw(fit = fit_state, draw = draw_index)
  )
  .validate_composition_draw(x, nrow(counts), ncol(counts), estimator$name)
}

#' Resolve sigma/rho for one draw from whichever scale input mode is active
#'
#' `scale_mode` is either `list(mode = "none")` (fixed `sigma_L`/etc., or
#' no scale information at all) or `list(mode = "data", scale_mat = <N x
#' R matrix>, log_transform = <TRUE/FALSE>, ci_level = , estimate_rho = ,
#' lower_zero = )` -- one matrix shape and one estimation pathway
#' regardless of whether the original input was a plain vector, a
#' replicate matrix, raw, or already log-scale.
#' @keywords internal
.resolve_draw_scale <- function(scale_mode, sample_index, log_props,
                                 sigma_L, sigma_U, rho_L, rho_U) {
  if (scale_mode$mode == "none") {
    return(list(sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U))
  }
  R <- ncol(scale_mode$scale_mat)
  replicate_index <- if (R == 1L) NULL else sample.int(R, length(sample_index), replace = TRUE)
  scale_log_s <- .select_scale_log(scale_mode$scale_mat, sample_index, replicate_index, scale_mode$log_transform)
  est <- .estimate_scale_log_bounds(
    log_props, scale_log_s,
    ci_level = scale_mode$ci_level,
    estimate_rho = scale_mode$estimate_rho,
    lower_zero = scale_mode$lower_zero
  )
  list(sigma_L = est$sigma_L, sigma_U = est$sigma_U, rho_L = est$rho_L, rho_U = est$rho_U)
}

#' One bootstrap attempt: resample, draw a composition, estimate bounds
#' @keywords internal
.prism_draw_attempt <- function(counts, estimator, fit_state, scale_mode,
                                 sigma_L, sigma_U, rho_L, rho_U,
                                 resample_samples, draw_index) {
  N <- ncol(counts)
  sample_index <- if (resample_samples) sample.int(N, N, replace = TRUE) else seq_len(N)
  counts_s <- counts[, sample_index, drop = FALSE]

  x <- .draw_one_composition(estimator, counts_s, fit_state, draw_index)
  logP <- log(x)
  Sigma_rel_s <- prism_relative_covariance(logP)

  resolved <- .resolve_draw_scale(scale_mode, sample_index, logP, sigma_L, sigma_U, rho_L, rho_U)
  bounds <- cov_bounds(
    Sigma_rel_s,
    sigma_L = resolved$sigma_L, sigma_U = resolved$sigma_U,
    rho_L = resolved$rho_L, rho_U = resolved$rho_U
  )
  list(
    Sigma_rel = Sigma_rel_s, lower = bounds$lower, upper = bounds$upper,
    sigma_L = resolved$sigma_L, sigma_U = resolved$sigma_U,
    rho_L = resolved$rho_L, rho_U = resolved$rho_U
  )
}

#' Run the full PRISM bootstrap loop
#'
#' Resamples samples (paired across counts and scale, per `bootstrap`) and
#' draws compositions, retrying rejected draws up to
#' `max_attempts_per_draw` times, and returns every accepted draw's bound
#' matrices for `aggregate.R` to summarize.
#'
#' @param counts Features-by-samples count matrix (post-filtering).
#' @param estimator A `prism_composition_estimator`.
#' @param scale_mode A list describing the active scale-information mode:
#'   `list(mode = "none")` (fixed `sigma_L`/etc., or no scale information)
#'   or `list(mode = "data", scale_mat = <N x R>, log_transform = ,
#'   ci_level = , estimate_rho = , lower_zero = )`.
#' @param sigma_L,sigma_U,rho_L,rho_U Fixed bounds, used when `scale_mode$mode == "none"`.
#' @param bootstrap One of `"none"` or `"both"` (see `.resolve_bootstrap_mode()`).
#' @param S Requested number of draws.
#' @param max_attempts_per_draw Retry budget per requested draw.
#' @param seed Reproducibility seed, or `NULL`.
#' @param verbose Emit progress messages.
#' @param stream Stream draws to temporary disk files instead of holding
#'   them in memory. `FALSE` (the default) still streams automatically
#'   when the estimated in-memory cost exceeds `stream_memory_limit`.
#' @param stream_memory_limit Byte threshold for the automatic streaming
#'   trigger; see `stream`.
#' @return A list: `storage` (`"memory"` or `"stream"`), `D`, `S_eff`, and
#'   either `lower`/`upper`/`rel` (D x D x S_eff arrays, `storage ==
#'   "memory"`) or `files`/`pair_index` (`storage == "stream"`, see
#'   `.stream_read_block()`), plus `diagnostics` (sampling accept/reject
#'   bookkeeping).
#' @keywords internal
.prism_bootstrap <- function(counts, estimator, scale_mode,
                              sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL,
                              bootstrap = "both", S = 1000L, max_attempts_per_draw = 100L,
                              seed = 1L, verbose = FALSE,
                              stream = FALSE, stream_memory_limit = 2e9) {
  bootstrap <- .resolve_bootstrap_mode(bootstrap)
  bootstrap_active <- bootstrap == "both"
  S <- .assert_positive_int(S, "S")
  max_attempts_per_draw <- .assert_positive_int(max_attempts_per_draw, "max_attempts_per_draw")

  D <- nrow(counts)
  N <- ncol(counts)

  resample_samples <- bootstrap_active && estimator$type != "custom"
  has_randomness <- bootstrap_active || estimator$type %in% c("dirichlet_multinomial", "custom")
  S_eff <- if (has_randomness) S else 1L

  fit_state <- NULL
  if (estimator$type == "custom") {
    fit_state <- estimator$options$fit(counts = counts, n_draws = S_eff, seed = seed, verbose = verbose)
  }

  if (!is.null(seed)) set.seed(seed)
  rng_seeds <- matrix(
    sample.int(.Machine$integer.max, S_eff * max_attempts_per_draw),
    nrow = S_eff, ncol = max_attempts_per_draw
  )

  stream_cfg <- .resolve_stream(stream, D, S_eff, stream_memory_limit)
  pair_index <- .upper_pairs(D)
  success <- FALSE
  if (stream_cfg$active) {
    handle <- .stream_open()
    on.exit({
      .stream_close(handle)
      if (!success) unlink(unlist(handle$files))
    }, add = TRUE)
  } else {
    lower_arr <- array(NA_real_, c(D, D, S_eff))
    upper_arr <- array(NA_real_, c(D, D, S_eff))
    rel_arr <- array(NA_real_, c(D, D, S_eff))
  }
  attempts_used <- integer(S_eff)
  rejection_reasons <- character(0)
  accepted <- logical(S_eff)

  for (s in seq_len(S_eff)) {
    n_attempts <- if (bootstrap_active) max_attempts_per_draw else 1L
    for (attempt in seq_len(n_attempts)) {
      set.seed(rng_seeds[s, attempt])
      result <- tryCatch(
        .prism_draw_attempt(
          counts, estimator, fit_state, scale_mode, sigma_L, sigma_U, rho_L, rho_U,
          resample_samples, draw_index = s
        ),
        error = function(e) {
          if (.is_retriable_draw_error(conditionMessage(e))) {
            structure(list(reason = conditionMessage(e)), class = "prism_draw_rejected")
          } else {
            stop(e)
          }
        }
      )
      attempts_used[s] <- attempt
      if (!inherits(result, "prism_draw_rejected")) {
        if (stream_cfg$active) {
          .stream_write_draw(handle, pair_index, result$lower, result$upper, result$Sigma_rel)
        } else {
          lower_arr[, , s] <- result$lower
          upper_arr[, , s] <- result$upper
          rel_arr[, , s] <- result$Sigma_rel
        }
        accepted[s] <- TRUE
        break
      }
      rejection_reasons <- c(rejection_reasons, result$reason)
    }
    if (verbose && !accepted[s]) {
      message("Draw ", s, " exhausted ", n_attempts, " attempts without an accepted resample.")
    }
  }

  n_accepted <- sum(accepted)
  n_attempted <- sum(attempts_used)
  n_rejected <- n_attempted - n_accepted
  diagnostics <- list(
    requested_draws = S_eff,
    accepted_draws = n_accepted,
    attempted_draws = n_attempted,
    rejected_attempts = n_rejected,
    acceptance_fraction = if (n_attempted > 0) n_accepted / n_attempted else NA_real_,
    rejection_fraction = if (n_attempted > 0) n_rejected / n_attempted else NA_real_,
    draws_requiring_replacement = sum(attempts_used > 1L & accepted),
    max_attempts_per_draw = max_attempts_per_draw,
    rejection_reasons = table(rejection_reasons),
    stream_active = stream_cfg$active,
    stream_estimated_bytes = stream_cfg$estimated_bytes
  )

  if (bootstrap_active && n_accepted < S_eff) {
    stop(
      "PRISM could not obtain ", S_eff, " accepted bootstrap draws within ",
      max_attempts_per_draw, " attempts each (", n_accepted, " accepted, ",
      n_rejected, " rejected). Rejection reasons: ",
      paste(names(diagnostics$rejection_reasons), diagnostics$rejection_reasons, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  if (n_accepted < 1L) {
    stop("PRISM obtained no valid draws.", call. = FALSE)
  }

  success <- TRUE
  if (stream_cfg$active) {
    list(
      storage = "stream", D = D, S_eff = n_accepted,
      files = handle$files, pair_index = pair_index,
      diagnostics = diagnostics
    )
  } else {
    list(
      storage = "memory", D = D, S_eff = n_accepted,
      lower = lower_arr[, , accepted, drop = FALSE],
      upper = upper_arr[, , accepted, drop = FALSE],
      rel = rel_arr[, , accepted, drop = FALSE],
      diagnostics = diagnostics
    )
  }
}
