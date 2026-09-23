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
      draws <- matrix(stats::rgamma(length(shape), shape = shape, rate = 1), nrow = nrow(shape))
      # The true Dirichlet support is strictly positive; a rate-1 Gamma
      # draw can underflow to exactly 0 for a very small shape parameter
      # (extreme concentration/pseudocount, or a genuinely near-absent
      # feature), which is a floating-point artifact, not a real 0
      # probability. Floor it instead of failing the draw.
      pmax(draws, .Machine$double.xmin)
    },
    custom = estimator$options$draw(fit = fit_state, draw = draw_index)
  )
  .validate_composition_draw(x, nrow(counts), ncol(counts), estimator$name)
}

#' Resolve sigma/rho for one draw from whichever scale input mode is active
#'
#' `scale_mode` is either `list(mode = "none")` (fixed `sigma_L`/etc., or
#' no scale information at all) or `list(mode = "data", scale_mat = <N x
#' R matrix>, ci_level = , estimate_rho = , lower_zero = )`. `scale_mat`
#' is always already on the log scale.
#' @keywords internal
.resolve_draw_scale <- function(scale_mode, sample_index, log_props,
                                 sigma_L, sigma_U, rho_L, rho_U) {
  if (scale_mode$mode == "none") {
    return(list(
      sigma_L = sigma_L, sigma_U = sigma_U,
      rho_L = rho_L, rho_U = rho_U, rho_witness = NULL
    ))
  }
  scale_log_s <- .select_scale_log(scale_mode$scale_mat, sample_index)
  est <- .estimate_scale_log_bounds(
    log_props, scale_log_s,
    ci_level = scale_mode$ci_level,
    estimate_rho = scale_mode$estimate_rho,
    lower_zero = scale_mode$lower_zero
  )
  list(
    sigma_L = est$sigma_L, sigma_U = est$sigma_U,
    rho_L = est$rho_L, rho_U = est$rho_U,
    rho_witness = est$rho_witness
  )
}

#' Generate one bootstrap draw
#' @keywords internal
.prism_draw <- function(counts, estimator, fit_state, scale_mode,
                        sigma_L, sigma_U, rho_L, rho_U,
                        resample_samples, draw_index,
                        rho_backend = "ecos", solver = NULL) {
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
    rho_L = resolved$rho_L, rho_U = resolved$rho_U,
    rho_witness = resolved$rho_witness,
    rho_backend = rho_backend, solver = solver
  )
  list(
    Sigma_rel = Sigma_rel_s, lower = bounds$lower, upper = bounds$upper,
    sigma_L = resolved$sigma_L, sigma_U = resolved$sigma_U,
    rho_L = resolved$rho_L, rho_U = resolved$rho_U
  )
}

.prism_seeded_draw <- function(draw_index, rng_seeds, counts, estimator,
                               fit_state, scale_mode,
                               sigma_L, sigma_U, rho_L, rho_U,
                               resample_samples, rho_backend, solver) {
  set.seed(rng_seeds[draw_index])
  .prism_draw(
    counts, estimator, fit_state, scale_mode,
    sigma_L, sigma_U, rho_L, rho_U,
    resample_samples, draw_index = draw_index,
    rho_backend = rho_backend, solver = solver
  )
}

.prism_seeded_draw_safe <- function(...) {
  tryCatch(
    list(
      ok = TRUE,
      result = get(".prism_seeded_draw", envir = asNamespace("prism"), inherits = FALSE)(...)
    ),
    error = function(e) list(ok = FALSE, message = conditionMessage(e))
  )
}

#' Run the full PRISM bootstrap loop
#'
#' Resamples samples (paired across counts and scale, per `bootstrap`) and
#' draws compositions. Every requested draw is generated exactly once; an
#' invalid or infeasible draw stops with its original error rather than being
#' replaced by a new resample.
#'
#' @param counts Features-by-samples count matrix (post-filtering).
#' @param estimator A `prism_composition_estimator`.
#' @param scale_mode A list describing the active scale-information mode:
#'   `list(mode = "none")` (fixed `sigma_L`/etc., or no scale information)
#'   or `list(mode = "data", scale_mat = <N x R>, ci_level = ,
#'   estimate_rho = , lower_zero = )`, with `scale_mat` on the log scale.
#' @param sigma_L,sigma_U,rho_L,rho_U Fixed bounds, used when `scale_mode$mode == "none"`.
#' @param bootstrap One of `"none"` or `"both"` (see `.resolve_bootstrap_mode()`).
#' @param S Requested number of draws.
#' @param seed Reproducibility seed, or `NULL`.
#' @param verbose Emit progress messages.
#' @param stream Stream draws to temporary disk files instead of holding
#'   them in memory. `FALSE` (the default) still streams automatically
#'   when the estimated in-memory cost exceeds `stream_memory_limit`.
#' @param stream_memory_limit Byte threshold for the automatic streaming
#'   trigger; see `stream`.
#' @param workers Number of bootstrap draws to evaluate concurrently. Results
#'   are collected in adaptive bounded batches, so streaming memory is bounded
#'   independently of `S`.
#' @param rho_backend Rho-support backend, `"ecos"` or `"cvxr"`.
#' @param solver Optional solver name forwarded to [cov_bounds()].
#' @return A list: `storage` (`"memory"` or `"stream"`), `D`, `S_eff`, and
#'   either `lower`/`upper`/`rel` (D x D x S_eff arrays, `storage ==
#'   "memory"`) or `files`/`pair_index` (`storage == "stream"`, see
#'   `.stream_read_block()`), plus computational `diagnostics`.
#' @keywords internal
.prism_bootstrap <- function(counts, estimator, scale_mode,
                              sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL,
                              bootstrap = "both", S = 1000L,
                              seed = 1L, verbose = FALSE,
                              stream = FALSE, stream_memory_limit = 2e9,
                              workers = 1L, rho_backend = "ecos",
                              solver = NULL) {
  bootstrap <- .resolve_bootstrap_mode(bootstrap)
  bootstrap_active <- bootstrap == "both"
  S <- .assert_positive_int(S, "S")
  workers <- .assert_positive_int(workers, "workers")
  rho_backend <- match.arg(rho_backend, c("ecos", "cvxr"))

  D <- nrow(counts)
  N <- ncol(counts)

  resample_samples <- bootstrap_active && estimator$type != "custom"
  scale_randomness <- scale_mode$mode == "data" && ncol(scale_mode$scale_mat) > 1L
  has_randomness <- bootstrap_active || scale_randomness || estimator$type %in% c("dirichlet_multinomial", "custom")
  S_eff <- if (has_randomness) S else 1L

  fit_state <- NULL
  if (estimator$type == "custom") {
    fit_state <- estimator$options$fit(counts = counts, n_draws = S_eff, seed = seed, verbose = verbose)
  }

  if (!is.null(seed)) set.seed(seed)
  rng_seeds <- sample.int(.Machine$integer.max, S_eff)

  stream_cfg <- .resolve_stream(stream, D, S_eff, stream_memory_limit)
  workers_used <- .resolve_parallel_workers(
    D, S_eff, workers, stream_memory_limit
  )
  parallel_active <- workers_used > 1L
  parallel_batch_size <- .resolve_parallel_batch_size(
    D, S_eff, workers_used, stream_memory_limit
  )
  cluster_type <- if (parallel_active) {
    if (.Platform$OS.type == "windows") "PSOCK" else "multicore"
  } else {
    "sequential"
  }
  cluster <- NULL
  if (parallel_active && cluster_type == "PSOCK") {
    cluster <- parallel::makeCluster(workers_used, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_libpaths <- .libPaths()
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      loadNamespace("prism")
      invisible(NULL)
    }, worker_libpaths)
  }

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
  batch_starts <- seq.int(1L, S_eff, by = parallel_batch_size)
  for (batch_start in batch_starts) {
    draw_indices <- seq.int(batch_start, min(S_eff, batch_start + parallel_batch_size - 1L))
    results <- if (parallel_active && cluster_type == "PSOCK") {
      parallel::parLapply(
        cluster, draw_indices, .prism_seeded_draw_safe,
        rng_seeds = rng_seeds, counts = counts, estimator = estimator,
        fit_state = fit_state, scale_mode = scale_mode,
        sigma_L = sigma_L, sigma_U = sigma_U,
        rho_L = rho_L, rho_U = rho_U,
        resample_samples = resample_samples,
        rho_backend = rho_backend, solver = solver
      )
    } else if (parallel_active) {
      parallel::mclapply(
        draw_indices, .prism_seeded_draw_safe,
        rng_seeds = rng_seeds, counts = counts, estimator = estimator,
        fit_state = fit_state, scale_mode = scale_mode,
        sigma_L = sigma_L, sigma_U = sigma_U,
        rho_L = rho_L, rho_U = rho_U,
        resample_samples = resample_samples,
        rho_backend = rho_backend, solver = solver,
        mc.preschedule = TRUE, mc.set.seed = FALSE,
        mc.cores = min(workers_used, length(draw_indices))
      )
    } else {
      lapply(
        draw_indices, .prism_seeded_draw,
        rng_seeds = rng_seeds, counts = counts, estimator = estimator,
        fit_state = fit_state, scale_mode = scale_mode,
        sigma_L = sigma_L, sigma_U = sigma_U,
        rho_L = rho_L, rho_U = rho_U,
        resample_samples = resample_samples,
        rho_backend = rho_backend, solver = solver
      )
    }
    if (parallel_active) {
      failed <- which(!vapply(results, function(x) isTRUE(x$ok), logical(1)))
      if (length(failed)) {
        stop(results[[failed[1L]]]$message, call. = FALSE)
      }
      results <- lapply(results, `[[`, "result")
    }
    for (k in seq_along(draw_indices)) {
      s <- draw_indices[k]
      result <- results[[k]]
      if (stream_cfg$active) {
        .stream_write_draw(handle, pair_index, result$lower, result$upper, result$Sigma_rel)
      } else {
        lower_arr[, , s] <- result$lower
        upper_arr[, , s] <- result$upper
        rel_arr[, , s] <- result$Sigma_rel
      }
    }
  }

  diagnostics <- list(
    requested_draws = S_eff,
    completed_draws = S_eff,
    stream_active = stream_cfg$active,
    stream_estimated_bytes = stream_cfg$estimated_bytes,
    workers_requested = workers,
    workers_used = workers_used,
    parallel_active = parallel_active,
    parallel_batch_size = parallel_batch_size,
    parallel_cluster_type = cluster_type
  )

  success <- TRUE
  if (stream_cfg$active) {
    list(
      storage = "stream", D = D, S_eff = S_eff,
      files = handle$files, pair_index = pair_index,
      diagnostics = diagnostics
    )
  } else {
    list(
      storage = "memory", D = D, S_eff = S_eff,
      lower = lower_arr,
      upper = upper_arr,
      rel = rel_arr,
      diagnostics = diagnostics
    )
  }
}
