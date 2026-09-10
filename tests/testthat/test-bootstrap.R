no_scale <- function() list(mode = "none")

test_that("resolve_bootstrap_mode normalizes logical and string forms", {
  expect_identical(.resolve_bootstrap_mode(TRUE), "both")
  expect_identical(.resolve_bootstrap_mode(FALSE), "none")
  expect_identical(.resolve_bootstrap_mode("both"), "both")
  expect_identical(.resolve_bootstrap_mode("none"), "none")
  expect_error(.resolve_bootstrap_mode("composition"), "bootstrap must be")
  expect_error(.resolve_bootstrap_mode("scale"), "bootstrap must be")
})

test_that("well-behaved data accepts every draw on the first attempt", {
  counts <- random_counts(4L, 20L, seed = 10L)
  fit <- .prism_bootstrap(
    counts, prism_composition_dirichlet(), no_scale(),
    sigma_L = 0.1, sigma_U = 0.5, bootstrap = "both", S = 25L, seed = 1L
  )
  d <- fit$diagnostics
  expect_identical(d$attempted_draws, d$accepted_draws + d$rejected_attempts)
  expect_equal(d$acceptance_fraction + d$rejection_fraction, 1)
  expect_identical(d$rejected_attempts, 0L)
  expect_identical(dim(fit$lower)[3], 25L)
})

test_that("a feature with exactly zero variance under fixed rho is rejected and eventually exhausts", {
  # Every sample has identical composition (columns are scalar multiples of
  # one vector), so every feature's relative log-variance is exactly zero
  # under every possible resample.
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(0.5), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      bootstrap = "both", S = 2L, max_attempts_per_draw = 3L, seed = 1L
    ),
    "zero-variance"
  )
})

test_that("without bootstrap, a deterministically degenerate draw still errors (no valid draws)", {
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  expect_error(
    .prism_bootstrap(
      counts, prism_composition_fixed(0.5), no_scale(),
      sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
      bootstrap = "none", S = 2L, seed = 1L
    ),
    "no valid draws"
  )
})

test_that("verbose messages on an exhausted draw", {
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), 12L), nrow = 4L)
  expect_message(
    tryCatch(
      .prism_bootstrap(
        counts, prism_composition_fixed(0.5), no_scale(),
        sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.3, 0, 0, 0), rho_U = c(0.3, 0, 0, 0),
        bootstrap = "both", S = 1L, max_attempts_per_draw = 2L, seed = 1L, verbose = TRUE
      ),
      error = function(e) NULL
    ),
    "exhausted"
  )
})

test_that("max_attempts_per_draw must be a positive integer", {
  counts <- random_counts(3L, 6L)
  expect_error(
    .prism_bootstrap(counts, prism_composition_dirichlet(), no_scale(), max_attempts_per_draw = 0),
    "max_attempts_per_draw"
  )
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

test_that("scale_log mode resamples paired with counts and yields finite bounds", {
  set.seed(21L)
  counts <- random_counts(3L, 15L, seed = 21L)
  u <- stats::rnorm(15L)
  scale_mode <- list(mode = "scale_log", scale_log = u, ci_level = 0, estimate_rho = TRUE, lower_zero = FALSE)
  fit <- .prism_bootstrap(counts, prism_composition_fixed(0.5), scale_mode, bootstrap = "both", S = 8L, seed = 3L)
  expect_true(all(is.finite(fit$lower)))
})

test_that("raw scale replicate mode resamples a replicate column per draw", {
  set.seed(22L)
  counts <- random_counts(3L, 10L, seed = 22L)
  scale_mat <- .validate_scale_input(matrix(stats::rlnorm(30L), nrow = 10L, ncol = 3L), N = 10L)
  scale_mode <- list(mode = "scale", scale_mat = scale_mat, ci_level = 0, estimate_rho = TRUE, lower_zero = FALSE)
  fit <- .prism_bootstrap(counts, prism_composition_fixed(0.5), scale_mode, bootstrap = "both", S = 6L, seed = 4L)
  expect_true(all(is.finite(fit$lower)))
  expect_identical(dim(fit$lower)[3], 6L)
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

test_that("an invalid custom draw fails immediately, without retrying", {
  est <- prism_composition_estimator(
    name = "broken",
    fit = function(counts, n_draws, seed, verbose) list(counts = counts),
    draw = function(fit, draw) matrix(0, nrow = nrow(fit$counts), ncol = ncol(fit$counts))
  )
  counts <- random_counts(3L, 8L, seed = 8L)
  expect_error(
    .prism_bootstrap(counts, est, no_scale(), bootstrap = "none", S = 2L, max_attempts_per_draw = 5L, seed = 1L),
    "returned an invalid draw"
  )
})
