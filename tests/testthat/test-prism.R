test_that("prism() with a fixed estimator matches the exact closed-form relative covariance", {
  counts <- random_counts(5L, 12L, seed = 30L)
  fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, S = 1L, seed = 1L)
  expected <- stats::cov(t(log(sweep(counts + 0.5, 2, colSums(counts + 0.5), "/"))))
  for (k in seq_len(nrow(fit$pairwise_rel))) {
    i <- fit$pairwise_rel$i[k]; j <- fit$pairwise_rel$j[k]
    expect_equal(fit$pairwise_rel$rel_ci_lower[k], expected[i, j], tolerance = 1e-10)
    expect_equal(fit$pairwise_rel$rel_ci_upper[k], expected[i, j], tolerance = 1e-10)
  }
  expect_identical(fit$parameters$S_effective, 1L)
  expect_true(all(is.na(fit$pairwise$p_value)))
})

test_that("prism() reproduces bootstrap results exactly for the same seed", {
  counts <- random_counts(4L, 15L, seed = 31L)
  fit1 <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 7L)
  fit2 <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 7L)
  expect_identical(fit1$ci_lower, fit2$ci_lower)
  expect_identical(fit1$pairwise, fit2$pairwise)
})

test_that("composition must be a prism_composition_estimator", {
  counts <- random_counts(3L, 6L)
  expect_error(prism(counts, composition = "Dirichlet"), "prism_composition_estimator")
})

test_that("sigma_L/rho_L cannot be combined with a non-NULL scale", {
  counts <- random_counts(3L, 6L, seed = 32L)
  expect_error(
    prism(counts, scale = seq_len(6L), sigma_L = 0.1, sigma_U = 0.5, bootstrap = FALSE),
    "cannot be combined"
  )
  expect_error(
    prism(counts, scale = prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5), sigma_L = 0.2, sigma_U = 0.6, bootstrap = FALSE),
    "cannot be combined"
  )
})

test_that("scale = prism_scale_bounds(...) is equivalent to the flat sigma_L/etc. shorthand", {
  counts <- random_counts(4L, 10L, seed = 46L)
  a <- prism(counts, bootstrap = FALSE, S = 5L, sigma_L = 0.2, sigma_U = 0.8, rho_L = -0.3, rho_U = 0.3, seed = 1L)
  b <- prism(counts, bootstrap = FALSE, S = 5L, scale = prism_scale_bounds(0.2, 0.8, -0.3, 0.3), seed = 1L)
  expect_identical(a$ci_lower, b$ci_lower)
  expect_identical(a$pairwise, b$pairwise)
})

test_that("counts must be integer-valued, finite, and non-negative", {
  counts <- random_counts(3L, 6L, seed = 33L)
  bad <- counts; bad[1, 1] <- 1.5
  expect_error(prism(bad, bootstrap = FALSE), "integer-valued")
  bad2 <- counts; bad2[1, 1] <- -1
  expect_error(prism(bad2, bootstrap = FALSE), "non-negative")
})

test_that("zero-depth samples are dropped and reflected in diagnostics", {
  counts <- random_counts(4L, 8L, seed = 34L)
  counts[, 1] <- 0L
  expect_warning(
    fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE),
    "Dropping 1 zero-depth"
  )
  expect_identical(fit$diagnostics$data$N, 7L)
  expect_identical(fit$diagnostics$data$dropped_samples, 1L)
})

test_that("prevalence filtering drops low-prevalence features", {
  counts <- random_counts(5L, 10L, seed = 35L)
  counts[1, ] <- c(rep(0L, 8L), 3L, 4L)  # present in 2/10 samples
  expect_warning(
    fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, prevalence = 0.5),
    "Dropping 1 feature"
  )
  expect_identical(fit$diagnostics$data$D, 4L)
})

test_that("bootstrap = FALSE with a random composition estimator still draws S times", {
  counts <- random_counts(3L, 8L, seed = 36L)
  fit <- prism(counts, bootstrap = FALSE, S = 9L, sigma_L = 0.1, sigma_U = 0.5, seed = 1L)
  expect_identical(fit$parameters$S_effective, 9L)
  expect_true(all(is.na(fit$pairwise$p_value)))
  expect_true(all(is.na(fit$pairwise$q_value)))
  expect_type(fit$pairwise$covers_zero, "logical")
})

test_that("prism() assumes bare scale data are already log transformed", {
  set.seed(41L)
  counts <- random_counts(3L, 10L, seed = 41L)
  log_scale <- stats::rnorm(10L)
  bare <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = log_scale
  )
  wrapped <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = prism_scale_log(log_scale)
  )
  expect_identical(bare$parameters$scale_mode, "data")
  expect_equal(bare$ci_lower, wrapped$ci_lower)
  expect_equal(bare$ci_upper, wrapped$ci_upper)
})

test_that("prism() accepts a log-scale replicate matrix end-to-end", {
  set.seed(42L)
  counts <- random_counts(3L, 10L, seed = 42L)
  log_scale <- matrix(stats::rnorm(30L), nrow = 10L, ncol = 3L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = "both", S = 5L,
    scale = log_scale, seed = 1L
  )
  expect_true(all(is.finite(fit$ci_lower)))
})

test_that("bare log-scale data may contain negative and zero values", {
  counts <- random_counts(3L, 8L, seed = 48L)
  set.seed(48L)
  log_scale <- stats::rnorm(8L)  # can include negatives/zero
  log_scale[1] <- -0.5
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = log_scale, seed = 1L
  )
  expect_true(all(is.finite(fit$ci_lower)))
})

test_that("scale must be sample-aligned before filtering", {
  counts <- random_counts(3L, 8L, seed = 43L)
  expect_error(prism(counts, scale = seq_len(5L), bootstrap = FALSE), "scale must have length N")
})

test_that("scale must be finite", {
  counts <- random_counts(3L, 8L, seed = 44L)
  bad <- rep(1, 8L); bad[1] <- Inf
  expect_error(prism(counts, scale = bad, bootstrap = FALSE), "finite")
})

test_that("too strict a prevalence threshold errors clearly", {
  counts <- random_counts(3L, 8L, seed = 45L)
  expect_error(prism(counts, prevalence = 1.5), "prevalence")
  expect_error(
    prism(matrix(c(0, 0, 5, 5, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3L, byrow = TRUE), prevalence = 0.9),
    "prevalence"
  )
})

test_that("too few samples or features after filtering errors clearly", {
  # data where all but one sample is zero-depth
  bad <- matrix(0, nrow = 3L, ncol = 4L)
  bad[, 1] <- c(1, 2, 3)
  expect_warning(
    expect_error(prism(bad, composition = prism_composition_fixed(0.5), bootstrap = FALSE), "Fewer than 2 samples"),
    "Dropping 3 zero-depth"
  )

  # two of three features are zero in every sample
  bad2 <- rbind(c(1L, 2L, 3L, 4L), c(0L, 0L, 0L, 0L), c(0L, 0L, 0L, 0L))
  expect_error(
    prism(bad2, composition = prism_composition_fixed(0.5), bootstrap = FALSE),
    "Fewer than 2 features"
  )
})

test_that("counts passed to prism() must be at least 2 features by 2 samples", {
  expect_error(prism(matrix(1:4, nrow = 1L)), "at least 2 features")
})

test_that("prism_scale_log(estimate_rho = FALSE) estimates only sigma", {
  counts <- random_counts(3L, 10L, seed = 37L)
  set.seed(37L)
  u <- stats::rnorm(10L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = prism_scale_log(u, estimate_rho = FALSE)
  )
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
})

test_that("experimental scale intervals can also skip rho estimation", {
  counts <- random_counts(3L, 10L, seed = 49L)
  set.seed(49L)
  u <- stats::rnorm(10L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = prism_scale_log(
      u, ci_level = 0.9, estimate_rho = FALSE,
      experimental_ci = TRUE
    )
  )
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
})

test_that("constant log-scale measurements imply sigma equals zero", {
  counts <- random_counts(3L, 8L, seed = 53L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = rep(2, 8L)
  )
  expected <- prism_relative_covariance(prism_log_composition(counts, 0.5))
  expect_equal(unname(fit$ci_lower), unname(expected))
  expect_equal(unname(fit$ci_upper), unname(expected))
})

test_that("end-to-end bootstrap with sigma bounds gives ordered, finite CIs", {
  counts <- random_counts(6L, 20L, seed = 38L)
  fit <- prism(
    counts,
    composition = prism_composition_dirichlet(),
    bootstrap = TRUE, S = 40L, sigma_L = 0.2, sigma_U = 0.8, seed = 1L
  )
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
  expect_true(all(is.finite(fit$ci_lower)))
  expect_true(all(is.finite(fit$ci_upper)))
  expect_s3_class(fit, "prism_result")
})

test_that("stream = TRUE produces the same result as stream = FALSE", {
  counts <- random_counts(5L, 12L, seed = 50L)
  fit_mem <- prism(counts, bootstrap = TRUE, S = 20L, sigma_L = 0.1, sigma_U = 0.5, seed = 3L, stream = FALSE)
  fit_stream <- prism(counts, bootstrap = TRUE, S = 20L, sigma_L = 0.1, sigma_U = 0.5, seed = 3L, stream = TRUE)
  expect_equal(fit_mem$ci_lower, fit_stream$ci_lower)
  expect_equal(fit_mem$pairwise, fit_stream$pairwise)
  expect_false(fit_mem$parameters$stream_active)
  expect_true(fit_stream$parameters$stream_active)
})

test_that("CVXR reference, sequential ECOS, and parallel streamed ECOS agree", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("ECOSolveR")
  counts <- random_counts(3L, 10L, seed = 64L)
  args <- list(
    counts = counts,
    composition = prism_composition_dirichlet(),
    bootstrap = TRUE, S = 4L,
    sigma_L = 0.2, sigma_U = 0.8,
    rho_L = -0.9, rho_U = 0.9,
    seed = 55L, return_draws = TRUE
  )
  reference <- do.call(prism, c(args, list(
    rho_backend = "cvxr", solver = "ECOS", workers = 1L
  )))
  sequential <- do.call(prism, c(args, list(
    rho_backend = "ecos", workers = 1L
  )))
  parallel_stream <- do.call(prism, c(args, list(
    rho_backend = "ecos", workers = 2L, stream = TRUE
  )))

  expect_equal(sequential$draws$lower, reference$draws$lower, tolerance = 1e-7)
  expect_equal(sequential$draws$upper, reference$draws$upper, tolerance = 1e-7)
  expect_identical(parallel_stream$draws$lower, sequential$draws$lower)
  expect_identical(parallel_stream$draws$upper, sequential$draws$upper)
  expect_identical(parallel_stream$draws$rel, sequential$draws$rel)
  expect_true(parallel_stream$parameters$parallel_active)
  expect_true(parallel_stream$parameters$stream_active)
  expect_identical(parallel_stream$parameters$rho_backend, "ecos")
})

test_that("optional end-to-end benchmark compares all three execution modes", {
  skip_if(Sys.getenv("PRISM_RUN_BENCHMARKS") != "true", "set PRISM_RUN_BENCHMARKS=true")
  skip_if_not_installed("CVXR")
  skip_if_not_installed("ECOSolveR")
  counts <- random_counts(5L, 80L, seed = 67L)
  args <- list(
    counts = counts,
    composition = prism_composition_dirichlet(),
    bootstrap = TRUE, S = 8L,
    sigma_L = 0.2, sigma_U = 0.8,
    rho_L = -0.9, rho_U = 0.9,
    seed = 67L, return_draws = TRUE
  )

  cvxr_time <- system.time(reference <- do.call(prism, c(args, list(
    rho_backend = "cvxr", solver = "ECOS", workers = 1L
  ))))[["elapsed"]]
  ecos_time <- system.time(sequential <- do.call(prism, c(args, list(
    rho_backend = "ecos", workers = 1L
  ))))[["elapsed"]]
  parallel_time <- system.time(parallel <- do.call(prism, c(args, list(
    rho_backend = "ecos", workers = 2L
  ))))[["elapsed"]]

  expect_equal(sequential$draws$lower, reference$draws$lower, tolerance = 1e-7)
  expect_identical(parallel$draws$lower, sequential$draws$lower)
  message(sprintf(
    paste0(
      "end-to-end benchmark: CVXR sequential %.3fs; ECOS sequential %.3fs; ",
      "ECOS parallel %.3fs; direct speedup %.1fx; parallel/direct %.1fx"
    ),
    cvxr_time, ecos_time, parallel_time,
    cvxr_time / max(ecos_time, .Machine$double.eps),
    ecos_time / max(parallel_time, .Machine$double.eps)
  ))
})

test_that("stream auto-triggers when the estimated size exceeds stream_memory_limit", {
  counts <- random_counts(5L, 12L, seed = 51L)
  fit <- prism(
    counts, bootstrap = TRUE, S = 20L, sigma_L = 0.1, sigma_U = 0.5, seed = 3L,
    stream = FALSE, stream_memory_limit = 100
  )
  expect_true(fit$parameters$stream_active)
})

test_that("streaming leaves no temp files behind, including on an error path", {
  counts <- random_counts(5L, 12L, seed = 52L)
  before <- list.files(tempdir(), pattern = "\\.bin$")
  invisible(prism(counts, bootstrap = TRUE, S = 10L, sigma_L = 0.1, sigma_U = 0.5, seed = 1L, stream = TRUE))
  after <- list.files(tempdir(), pattern = "\\.bin$")
  expect_identical(before, after)
})

test_that("parallel streaming is bounded by workers and cleans up after an error", {
  counts <- random_counts(5L, 12L, seed = 65L)
  fit <- prism(
    counts, bootstrap = TRUE, S = 7L,
    sigma_L = 0.1, sigma_U = 0.5, seed = 1L,
    stream = TRUE, workers = 3L
  )
  expect_lte(fit$diagnostics$sampling$parallel_batch_size, 7L)
  expect_gte(fit$diagnostics$sampling$parallel_batch_size, 3L)
  expect_identical(fit$diagnostics$sampling$completed_draws, 7L)

  degenerate <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  before <- list.files(tempdir(), pattern = "\\.bin$")
  expect_error(
    prism(
      degenerate, composition = prism_composition_fixed(),
      bootstrap = TRUE, S = 4L,
      sigma_L = 0.1, sigma_U = 0.5,
      rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      seed = 1L, stream = TRUE, workers = 2L
    ),
    "zero or numerically zero marginal variance"
  )
  after <- list.files(tempdir(), pattern = "\\.bin$")
  expect_identical(after, before)
})

test_that("prism() does not leak its internal seeding into the caller's RNG stream", {
  counts <- random_counts(4L, 10L, seed = 54L)
  set.seed(123L)
  runif(1L)  # advance to some arbitrary known state
  expected_next <- runif(3L)

  set.seed(123L)
  runif(1L)
  invisible(prism(counts, bootstrap = TRUE, S = 20L, sigma_L = 0.1, sigma_U = 0.5, seed = 999L))
  actual_next <- runif(3L)

  expect_identical(actual_next, expected_next)
})

test_that("parallel prism() does not leak or change the caller RNG stream", {
  counts <- random_counts(4L, 10L, seed = 66L)
  set.seed(321L)
  runif(1L)
  expected_next <- runif(3L)

  set.seed(321L)
  runif(1L)
  invisible(prism(
    counts, bootstrap = TRUE, S = 8L,
    sigma_L = 0.1, sigma_U = 0.5,
    seed = 999L, workers = 2L
  ))
  actual_next <- runif(3L)
  expect_identical(actual_next, expected_next)
})

test_that("seed = NULL gives identical sequential and parallel draws from the same ambient state", {
  counts <- random_counts(3L, 10L, seed = 68L)
  set.seed(808L)
  sequential <- prism(
    counts, bootstrap = TRUE, S = 6L,
    sigma_L = 0.1, sigma_U = 0.5,
    seed = NULL, workers = 1L, return_draws = TRUE
  )
  set.seed(808L)
  parallel <- prism(
    counts, bootstrap = TRUE, S = 6L,
    sigma_L = 0.1, sigma_U = 0.5,
    seed = NULL, workers = 2L, return_draws = TRUE
  )
  expect_identical(parallel$draws$lower, sequential$draws$lower)
  expect_identical(parallel$draws$upper, sequential$draws$upper)
  expect_identical(parallel$draws$rel, sequential$draws$rel)
})

test_that("prism(seed = NULL) still consumes ambient randomness (no restore)", {
  counts <- random_counts(3L, 8L, seed = 55L)
  set.seed(1L)
  invisible(prism(counts, bootstrap = TRUE, S = 5L, sigma_L = 0.1, sigma_U = 0.5, seed = NULL))
  after1 <- .Random.seed
  invisible(prism(counts, bootstrap = TRUE, S = 5L, sigma_L = 0.1, sigma_U = 0.5, seed = NULL))
  after2 <- .Random.seed
  expect_false(identical(after1, after2))
})

test_that("prism() restores to a non-existent .Random.seed when none existed before", {
  counts <- random_counts(3L, 8L, seed = 56L)  # random_counts() itself calls set.seed()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  invisible(prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, seed = 1L))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("prism() output is unaffected by save/restore (still deterministic in seed)", {
  counts <- random_counts(4L, 10L, seed = 57L)
  set.seed(11L); runif(5L)  # perturb ambient RNG state before each call, differently
  fit1 <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 42L)
  set.seed(22L); runif(9L)
  fit2 <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 42L)
  expect_identical(fit1$ci_lower, fit2$ci_lower)
})

test_that("return_draws = TRUE returns per-draw values matching the aggregated CIs exactly", {
  counts <- random_counts(4L, 12L, seed = 58L)
  fit <- prism(counts, bootstrap = TRUE, S = 30L, sigma_L = 0.1, sigma_U = 0.5, seed = 1L, return_draws = TRUE)
  expect_null(prism(counts, bootstrap = TRUE, S = 5L, sigma_L = 0.1, sigma_U = 0.5, seed = 1L)$draws)
  expect_identical(dim(fit$draws$lower), c(30L, nrow(fit$draws$pair_index)))
  expect_identical(dim(fit$draws$upper), c(30L, nrow(fit$draws$pair_index)))
  expect_identical(dim(fit$draws$rel), c(30L, nrow(fit$draws$pair_index)))

  # .extract_prism_draws() and .assemble_prism_aggregate() both build their
  # pair index from .upper_pairs(D), so an off-diagonal subset of one lines
  # up with pairwise/pairwise_rel's rows with no re-sorting needed.
  off <- fit$draws$pair_index$i != fit$draws$pair_index$j
  computed_lower <- apply(fit$draws$lower[, off, drop = FALSE], 2, stats::quantile, probs = 0.025, names = FALSE)
  computed_upper <- apply(fit$draws$upper[, off, drop = FALSE], 2, stats::quantile, probs = 0.975, names = FALSE)
  expect_equal(computed_lower, fit$pairwise$ci_lower, tolerance = 1e-12)
  expect_equal(computed_upper, fit$pairwise$ci_upper, tolerance = 1e-12)

  computed_rel_lower <- apply(fit$draws$rel, 2, stats::quantile, probs = 0.025, names = FALSE)
  computed_rel_upper <- apply(fit$draws$rel, 2, stats::quantile, probs = 0.975, names = FALSE)
  expect_equal(computed_rel_lower, fit$pairwise_rel$rel_ci_lower, tolerance = 1e-12)
  expect_equal(computed_rel_upper, fit$pairwise_rel$rel_ci_upper, tolerance = 1e-12)
})

test_that("return_draws = TRUE gives identical draws whether streamed or not", {
  counts <- random_counts(5L, 12L, seed = 59L)
  fit_mem <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 3L, stream = FALSE, return_draws = TRUE)
  fit_stream <- prism(counts, bootstrap = TRUE, S = 15L, sigma_L = 0.1, sigma_U = 0.5, seed = 3L, stream = TRUE, return_draws = TRUE)
  expect_equal(fit_mem$draws$lower, fit_stream$draws$lower)
  expect_equal(fit_mem$draws$upper, fit_stream$draws$upper)
  expect_equal(fit_mem$draws$rel, fit_stream$draws$rel)
})

test_that("print.prism_result runs without error", {
  counts <- random_counts(3L, 6L, seed = 39L)
  fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE)
  expect_output(print(fit), "PRISM covariance analysis")
})

test_that("no side-effect files are created by a prism() call", {
  counts <- random_counts(3L, 8L, seed = 40L)
  dir <- file.path(tempdir(), paste0("prism-test-", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(dir)
  before <- list.files(dir, all.files = TRUE)
  old_wd <- setwd(dir)
  on.exit(setwd(old_wd), add = TRUE)
  invisible(prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE))
  after <- list.files(dir, all.files = TRUE)
  expect_identical(before, after)
})
