test_that("disk-backed streaming matches full stochastic output", {
  counts <- matrix(
    c(
      10, 20, 30, 40,
      20, 10, 25, 15,
      5, 8, 12, 9
    ),
    nrow = 3,
    byrow = TRUE
  )
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 20,
    sigma_L = 0.2,
    sigma_U = 0.8,
    seed = 42,
    verbose = FALSE,
    return_rel_draws = TRUE,
    return_bound_draws = TRUE
  )

  before <- list.files(tempdir(), pattern = "^prism_.*[.]bin$", full.names = TRUE)
  full <- do.call(prism, c(common, list(stream = FALSE)))
  streamed <- do.call(prism, c(common, list(stream = 5)))
  after <- list.files(tempdir(), pattern = "^prism_.*[.]bin$", full.names = TRUE)

  expect_equal(streamed$ci_lower, full$ci_lower, tolerance = 1e-12)
  expect_equal(streamed$ci_upper, full$ci_upper, tolerance = 1e-12)
  expect_equal(streamed$pairwise, full$pairwise, tolerance = 1e-12)
  expect_equal(streamed$pairwise_rel, full$pairwise_rel, tolerance = 1e-12)
  expect_equal(streamed$rel_draws, full$rel_draws, tolerance = 1e-12)
  expect_equal(
    streamed$pairwise_bound_draws,
    full$pairwise_bound_draws,
    tolerance = 1e-12
  )
  expect_identical(
    streamed$pairwise_bound_draws$pair_index[, c("i", "j")],
    streamed$pairwise[, c("i", "j")]
  )
  expect_equal(
    apply(streamed$pairwise_bound_draws$lower, 2, stats::quantile, probs = 0.025),
    streamed$pairwise$ci_lower,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    apply(streamed$pairwise_bound_draws$upper, 2, stats::quantile, probs = 0.975),
    streamed$pairwise$ci_upper,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    streamed$id_region_width_median,
    full$id_region_width_median,
    tolerance = 1e-12
  )
  expect_equal(streamed$lower_iqr, full$lower_iqr, tolerance = 1e-12)
  expect_equal(streamed$upper_iqr, full$upper_iqr, tolerance = 1e-12)
  expect_true(streamed$parameters$stream_active)
  expect_equal(streamed$parameters$stream_pair_block_size, 1L)
  expect_identical(streamed$parameters$stream_storage, "temporary_binary_files")
  expect_setequal(after, before)
})

test_that("multicore streaming matches serial streaming", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not(future::supportsMulticore(), "forked workers are unavailable")
  counts <- matrix(
    c(
      10, 20, 30, 40,
      20, 10, 25, 15,
      5, 8, 12, 9
    ),
    nrow = 3,
    byrow = TRUE
  )
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 20,
    scale_log = c(-0.4, -0.1, 0.2, 0.5),
    scale_log_ci_level = 0,
    seed = 42,
    stream = 5
  )

  serial <- do.call(prism, c(common, list(workers = 1)))
  parallel <- do.call(prism, c(common, list(workers = 2)))

  expect_equal(parallel$ci_lower, serial$ci_lower, tolerance = 1e-12)
  expect_equal(parallel$ci_upper, serial$ci_upper, tolerance = 1e-12)
  expect_equal(parallel$pairwise, serial$pairwise, tolerance = 1e-12)
  expect_identical(parallel$diagnostics$sampling$rejected_attempts, 0L)
})

test_that("disk-backed streaming matches full high-resolution refinement", {
  skip_if_not_installed("CVXR")
  skip_if_not("ECOS" %in% CVXR::installed_solvers(), "ECOS is unavailable")
  counts <- rbind(
    c(10, 20, 15, 25, 12, 30, 18, 22),
    c(30, 10, 25, 15, 28, 10, 22, 18)
  )
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 4,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-1, -1),
    rho_U = c(0.2, 1),
    seed = 42,
    verbose = FALSE,
    high_resolution = TRUE,
    high_resolution_policy = "all",
    solver = "ECOS"
  )

  full <- do.call(prism, c(common, list(stream = FALSE)))
  streamed <- do.call(prism, c(common, list(stream = 2)))

  expect_equal(streamed$ci_lower, full$ci_lower, tolerance = 1e-7)
  expect_equal(streamed$ci_upper, full$ci_upper, tolerance = 1e-7)
  expect_equal(streamed$pairwise, full$pairwise, tolerance = 1e-7)
  expect_equal(
    streamed$diagnostics$high_resolution,
    full$diagnostics$high_resolution,
    tolerance = 1e-7
  )
})

test_that("scale-log rho witnesses preserve rank-deficient streamed equivalence", {
  counts <- matrix(
    c(
      10, 12, 14, 16, 18,
      18, 15, 13, 11, 9,
      5, 8, 6, 10, 7,
      20, 18, 24, 21, 26,
      7, 9, 12, 8, 11,
      13, 10, 16, 14, 19
    ),
    nrow = 6,
    byrow = TRUE
  )
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 8,
    scale_log = c(-0.6, -0.1, 0.2, 0.4, 0.8),
    scale_log_ci_level = 0.8,
    solver = "intentionally-unavailable",
    seed = 91,
    verbose = FALSE
  )

  full <- do.call(prism, c(common, list(stream = FALSE)))
  streamed <- do.call(prism, c(common, list(stream = 3)))

  expect_equal(streamed$ci_lower, full$ci_lower, tolerance = 1e-12)
  expect_equal(streamed$ci_upper, full$ci_upper, tolerance = 1e-12)
  expect_equal(streamed$pairwise, full$pairwise, tolerance = 1e-12)
  expect_true(streamed$parameters$stream_active)
})
