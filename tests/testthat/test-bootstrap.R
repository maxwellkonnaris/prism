no_scale <- function() list(mode = "none")

test_that("resolve_bootstrap_mode normalizes logical and string forms", {
  expect_identical(.resolve_bootstrap_mode(TRUE), "both")
  expect_identical(.resolve_bootstrap_mode(FALSE), "none")
  expect_identical(.resolve_bootstrap_mode("both"), "both")
  expect_identical(.resolve_bootstrap_mode("none"), "none")
  expect_error(.resolve_bootstrap_mode("composition"), "bootstrap must be")
  expect_error(.resolve_bootstrap_mode("scale"), "bootstrap must be")
})

test_that("extreme Dirichlet draws fail if flooring leaves zero marginal variance", {
  counts <- random_counts(4L, 10L, seed = 60L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_dirichlet(concentration = 1e-6), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, bootstrap = "none", S = 30L, seed = 1L
    ),
    "zero or numerically zero marginal variance"
  )
})

test_that("draw_one_composition floors an underflowed Dirichlet draw instead of erroring", {
  est <- prism_composition_dirichlet(concentration = 1e-6, pseudocount = 1e-6)
  counts <- matrix(0L, nrow = 3L, ncol = 4L)  # all-zero counts maximize underflow risk
  set.seed(1L)
  x <- .draw_one_composition(est, counts, fit_state = NULL, draw_index = 1L)
  expect_true(all(x > 0))
  expect_true(all(is.finite(x)))
})

test_that("well-behaved data completes every requested draw", {
  counts <- random_counts(4L, 20L, seed = 10L)
  fit <- .prism_bootstrap(
    counts, prism_composition_dirichlet(), no_scale(),
    sigma_L = 0.1, sigma_U = 0.5, bootstrap = "both", S = 25L, seed = 1L
  )
  d <- fit$diagnostics
  expect_identical(d$requested_draws, 25L)
  expect_identical(d$completed_draws, 25L)
  expect_identical(dim(fit$lower)[3], 25L)
})

test_that("a feature with exactly zero variance under fixed rho fails immediately", {
  # Every sample has identical composition (columns are scalar multiples of
  # one vector), so every feature's relative log-variance is exactly zero
  # under every possible resample.
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(0.5), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      bootstrap = "both", S = 2L, seed = 1L
    ),
    "zero or numerically zero marginal variance"
  )
})

test_that("without bootstrap, a deterministically degenerate draw preserves its original error", {
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(0.5), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      bootstrap = "none", S = 2L, seed = 1L
    ),
    "zero or numerically zero marginal variance"
  )
})

test_that("a streamed run cleans up its temp files when a draw errors", {
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  before <- list.files(tempdir(), pattern = "\\.bin$")
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(0.5), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      bootstrap = "both", S = 2L, seed = 1L, stream = TRUE
    ),
    "zero or numerically zero marginal variance"
  )
  after <- list.files(tempdir(), pattern = "\\.bin$")
  expect_identical(before, after)
})

test_that("bootstrap = 'none' with a fixed estimator collapses to a single deterministic draw", {
  counts <- random_counts(3L, 6L, seed = 5L)
  fit <- .prism_bootstrap(
    counts, prism_composition_fixed(0.5), no_scale(),
    sigma_L = 0.2, sigma_U = 0.2, bootstrap = "none", S = 500L, seed = 1L
  )
  expect_identical(dim(fit$lower)[3], 1L)
  expected <- prism_relative_covariance(prism_log_composition(counts, 0.5))
  expect_equal(fit$rel[, , 1], expected)
})

test_that("bootstrap = 'none' with a random estimator still draws S times", {
  counts <- random_counts(3L, 6L, seed = 5L)
  fit <- .prism_bootstrap(
    counts, prism_composition_dirichlet(), no_scale(),
    sigma_L = 0.2, sigma_U = 0.5, bootstrap = "none", S = 7L, seed = 1L
  )
  expect_identical(dim(fit$lower)[3], 7L)
  expect_false(isTRUE(all.equal(fit$rel[, , 1], fit$rel[, , 2])))
})

test_that("same seed reproduces identical draws exactly", {
  counts <- random_counts(4L, 15L, seed = 6L)
  fit1 <- .prism_bootstrap(
    counts, prism_composition_dirichlet(), no_scale(),
    sigma_L = 0.1, sigma_U = 0.5, bootstrap = "both", S = 10L, seed = 42L
  )
  fit2 <- .prism_bootstrap(
    counts, prism_composition_dirichlet(), no_scale(),
    sigma_L = 0.1, sigma_U = 0.5, bootstrap = "both", S = 10L, seed = 42L
  )
  expect_identical(fit1$lower, fit2$lower)
  expect_identical(fit1$upper, fit2$upper)
})

test_that("parallel bootstrap is exactly invariant to worker count", {
  counts <- random_counts(4L, 15L, seed = 61L)
  args <- list(
    counts = counts,
    estimator = prism_composition_dirichlet(),
    scale_mode = no_scale(),
    sigma_L = 0.2, sigma_U = 0.8,
    rho_L = -0.9, rho_U = 0.9,
    bootstrap = "both", S = 8L, seed = 42L
  )
  sequential <- do.call(.prism_bootstrap, c(args, list(workers = 1L)))
  parallel <- do.call(.prism_bootstrap, c(args, list(workers = 2L)))

  expect_identical(parallel$lower, sequential$lower)
  expect_identical(parallel$upper, sequential$upper)
  expect_identical(parallel$rel, sequential$rel)
  expect_false(sequential$diagnostics$parallel_active)
  expect_true(parallel$diagnostics$parallel_active)
  expect_identical(parallel$diagnostics$workers_used, 2L)
  expect_lte(parallel$diagnostics$parallel_batch_size, 8L)
  expect_gte(parallel$diagnostics$parallel_batch_size, 2L)
})

test_that("workers are validated and collapse to one for a single effective draw", {
  counts <- random_counts(3L, 8L, seed = 62L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(), no_scale(),
      bootstrap = "none", S = 1L, workers = 0L
    ),
    "workers must be one positive integer"
  )
  fit <- .prism_bootstrap(
    counts, prism_composition_fixed(), no_scale(),
    bootstrap = "none", S = 100L, workers = 4L
  )
  expect_identical(fit$S_eff, 1L)
  expect_identical(fit$diagnostics$workers_requested, 4L)
  expect_identical(fit$diagnostics$workers_used, 1L)
  expect_false(fit$diagnostics$parallel_active)
})

test_that("pre-logged scale data resamples paired with counts and yields finite bounds", {
  set.seed(21L)
  counts <- random_counts(3L, 15L, seed = 21L)
  u <- stats::rnorm(15L)
  scale_mat <- .validate_scale_input(u, N = 15L)
  scale_mode <- list(mode = "data", scale_mat = scale_mat, ci_level = 0, estimate_rho = TRUE, lower_zero = TRUE)
  fit <- .prism_bootstrap(counts, prism_composition_fixed(0.5), scale_mode, bootstrap = "both", S = 8L, seed = 3L)
  expect_true(all(is.finite(fit$lower)))
})

test_that("parallel paired scale resampling matches sequential draws", {
  set.seed(63L)
  counts <- random_counts(3L, 15L, seed = 63L)
  scale_mat <- .validate_scale_input(stats::rnorm(15L), N = 15L)
  scale_mode <- list(
    mode = "data", scale_mat = scale_mat,
    ci_level = 0, estimate_rho = TRUE, lower_zero = TRUE
  )
  sequential <- .prism_bootstrap(
    counts, prism_composition_fixed(), scale_mode,
    bootstrap = "both", S = 6L, seed = 3L, workers = 1L
  )
  parallel <- .prism_bootstrap(
    counts, prism_composition_fixed(), scale_mode,
    bootstrap = "both", S = 6L, seed = 3L, workers = 2L
  )
  expect_identical(parallel$lower, sequential$lower)
  expect_identical(parallel$upper, sequential$upper)
  expect_identical(parallel$rel, sequential$rel)
})

test_that("paired scale estimates carry their feasible rho witness into each draw", {
  counts <- random_counts(4L, 12L, seed = 23L)
  logP <- prism_log_composition(counts, pseudocount = 0.5)
  scale_log <- as.numeric(logP[1, ] + stats::rnorm(12L, sd = 0.1))
  scale_mode <- list(
    mode = "data",
    scale_mat = matrix(scale_log, ncol = 1L),
    ci_level = 0.95,
    estimate_rho = TRUE,
    lower_zero = TRUE
  )
  resolved <- .resolve_draw_scale(
    scale_mode, seq_len(12L), logP,
    sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL
  )
  expect_true(all(resolved$rho_L <= resolved$rho_witness))
  expect_true(all(resolved$rho_witness <= resolved$rho_U))
  expect_true(.rho_box_feasible(
    prism_relative_covariance(logP), resolved$rho_L, resolved$rho_U,
    witness = resolved$rho_witness
  ))
})

test_that("log-scale replicate data propagates within-subject replicate uncertainty", {
  set.seed(22L)
  counts <- random_counts(3L, 10L, seed = 22L)
  scale_mat <- .validate_scale_input(matrix(stats::rnorm(30L), nrow = 10L, ncol = 3L), N = 10L)
  scale_mode <- list(mode = "data", scale_mat = scale_mat, ci_level = 0, estimate_rho = TRUE, lower_zero = TRUE)
  fit <- .prism_bootstrap(counts, prism_composition_fixed(0.5), scale_mode, bootstrap = "both", S = 6L, seed = 4L)
  expect_true(all(is.finite(fit$lower)))
  expect_identical(dim(fit$lower)[3], 6L)
})

test_that("technical replicates create S draws even without subject bootstrap", {
  counts <- random_counts(3L, 10L, seed = 24L)
  scale_mat <- .validate_scale_input(
    matrix(stats::rnorm(30L), nrow = 10L, ncol = 3L), N = 10L
  )
  scale_mode <- list(
    mode = "data", scale_mat = scale_mat,
    ci_level = 0, estimate_rho = TRUE, lower_zero = TRUE
  )
  fit <- .prism_bootstrap(
    counts, prism_composition_fixed(0.5), scale_mode,
    bootstrap = "none", S = 6L, seed = 5L
  )
  expect_identical(fit$S_eff, 6L)
  expect_false(isTRUE(all.equal(fit$lower[, , 1], fit$lower[, , 2])))
})

test_that("a custom composition estimator's fit() runs once and draw() runs per index", {
  fit_calls <- 0L
  draw_calls <- integer(0)
  est <- prism_composition_estimator(
    name = "counter",
    fit = function(counts, n_draws, seed, verbose) {
      fit_calls <<- fit_calls + 1L
      list(counts = counts)
    },
    draw = function(fit, draw) {
      draw_calls <<- c(draw_calls, draw)
      prism_closure(fit$counts, 0.5)
    }
  )
  counts <- random_counts(3L, 8L, seed = 7L)
  result <- .prism_bootstrap(counts, est, no_scale(), sigma_L = 0.1, sigma_U = 0.5, bootstrap = "none", S = 4L, seed = 1L)
  expect_identical(fit_calls, 1L)
  expect_identical(sort(draw_calls), 1:4)
  expect_identical(dim(result$lower)[3], 4L)
})

test_that("a seeded custom estimator is invariant to parallel execution", {
  est <- prism_composition_estimator(
    name = "parallel-custom",
    fit = function(counts, n_draws, seed, verbose) list(counts = counts),
    draw = function(fit, draw) {
      shape <- fit$counts + 0.5
      matrix(stats::rgamma(length(shape), shape = shape), nrow = nrow(shape))
    }
  )
  counts <- random_counts(3L, 8L, seed = 70L)
  sequential <- .prism_bootstrap(
    counts, est, no_scale(), sigma_L = 0.1, sigma_U = 0.5,
    bootstrap = "none", S = 6L, seed = 70L, workers = 1L
  )
  parallel <- .prism_bootstrap(
    counts, est, no_scale(), sigma_L = 0.1, sigma_U = 0.5,
    bootstrap = "none", S = 6L, seed = 70L, workers = 2L
  )
  expect_identical(parallel$lower, sequential$lower)
  expect_identical(parallel$upper, sequential$upper)
  expect_identical(parallel$rel, sequential$rel)
})

test_that("an invalid custom draw fails immediately", {
  est <- prism_composition_estimator(
    name = "broken",
    fit = function(counts, n_draws, seed, verbose) list(counts = counts),
    draw = function(fit, draw) matrix(0, nrow = nrow(fit$counts), ncol = ncol(fit$counts))
  )
  counts <- random_counts(3L, 8L, seed = 8L)
  expect_error(
    .prism_bootstrap(counts, est, no_scale(), bootstrap = "none", S = 2L, seed = 1L),
    "returned an invalid draw"
  )
})
