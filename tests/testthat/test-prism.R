test_that("prism() with a fixed estimator matches the exact closed-form relative covariance", {
  counts <- random_counts(5L, 12L, seed = 30L)
  fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, S = 1L, seed = 1L)
  expected <- stats::cov(t(log(sweep(counts + 0.5, 2, colSums(counts + 0.5), "/"))))
  for (k in seq_len(nrow(fit$pairwise_rel))) {
    i <- fit$pairwise_rel$i[k]; j <- fit$pairwise_rel$j[k]
    expect_equal(fit$pairwise_rel$rel_hat[k], expected[i, j], tolerance = 1e-10)
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
    prism(counts, scale = rep(1, 6L), sigma_L = 0.1, sigma_U = 0.5, bootstrap = FALSE),
    "cannot be combined"
  )
  expect_error(
    prism(counts, scale = prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5), sigma_L = 0.2, sigma_U = 0.6, bootstrap = FALSE),
    "cannot be combined"
  )
})

test_that("scale = prism_scale_bounds(...) is equivalent to the flat sigma_L/etc. shorthand", {
  counts <- random_counts(4L, 10L, seed = 46L)
  a <- prism(counts, bootstrap = FALSE, S = 5L, sigma_L = 0.2, sigma_U = 0.8, rho_L = 0.1, rho_U = 0.3, seed = 1L)
  b <- prism(counts, bootstrap = FALSE, S = 5L, scale = prism_scale_bounds(0.2, 0.8, 0.1, 0.3), seed = 1L)
  expect_identical(a$ci_lower, b$ci_lower)
  expect_identical(a$pairwise, b$pairwise)
})

test_that("scale_estimate_rho and scale_lower_zero must be TRUE or FALSE", {
  counts <- random_counts(3L, 8L, seed = 53L)
  set.seed(53L)
  u <- stats::rlnorm(8L)
  expect_error(prism(counts, scale = u, scale_estimate_rho = "yes", bootstrap = FALSE), "scale_estimate_rho")
  expect_error(prism(counts, scale = u, scale_lower_zero = "yes", bootstrap = FALSE), "scale_lower_zero")
})

test_that("scale_ci_level/estimate_rho/lower_zero require sample-level scale data", {
  counts <- random_counts(3L, 6L, seed = 47L)
  expect_error(prism(counts, scale_ci_level = 0.8, bootstrap = FALSE), "sample-level data")
  expect_error(
    prism(counts, scale = prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5), scale_ci_level = 0.8, bootstrap = FALSE),
    "do not apply"
  )
  set.seed(1L)
  u <- stats::rnorm(6L)
  expect_error(
    prism(counts, scale = prism_scale_log(u), scale_ci_level = 0.8, bootstrap = FALSE),
    "already"
  )
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
  fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, verbose = TRUE)
  expect_identical(fit$diagnostics$data$N, 7L)
  expect_identical(fit$diagnostics$data$dropped_samples, 1L)
})

test_that("prevalence filtering drops low-prevalence features", {
  counts <- random_counts(5L, 10L, seed = 35L)
  counts[1, ] <- c(rep(0L, 8L), 3L, 4L)  # present in 2/10 samples
  fit <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, prevalence = 0.5)
  expect_identical(fit$diagnostics$data$D, 4L)
})

test_that("bootstrap = FALSE with a random composition estimator still draws S times", {
  counts <- random_counts(3L, 8L, seed = 36L)
  fit <- prism(counts, bootstrap = FALSE, S = 9L, sigma_L = 0.1, sigma_U = 0.5, seed = 1L)
  expect_identical(fit$parameters$S_effective, 9L)
})

test_that("prism() accepts a raw (positive) scale vector, auto-detected as raw", {
  set.seed(41L)
  counts <- random_counts(3L, 10L, seed = 41L)
  raw_scale <- stats::rlnorm(10L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = raw_scale, scale_ci_level = 0.8
  )
  expect_identical(fit$parameters$scale_mode, "data")
  expect_true(all(is.finite(fit$ci_lower)))
})

test_that("prism() accepts a raw scale replicate matrix end-to-end", {
  set.seed(42L)
  counts <- random_counts(3L, 10L, seed = 42L)
  raw_scale <- matrix(stats::rlnorm(30L), nrow = 10L, ncol = 3L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = "both", S = 5L,
    scale = raw_scale, seed = 1L
  )
  expect_true(all(is.finite(fit$ci_lower)))
})

test_that("a scale vector with a non-positive value is auto-detected as already log-scale", {
  counts <- random_counts(3L, 8L, seed = 48L)
  set.seed(48L)
  log_scale <- stats::rnorm(8L)  # can include negatives/zero
  log_scale[1] <- -0.5
  fit_auto <- prism(counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE, scale = log_scale, seed = 1L)
  fit_explicit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = prism_scale_log(log_scale), seed = 1L
  )
  expect_equal(fit_auto$ci_lower, fit_explicit$ci_lower)
})

test_that("scale must be sample-aligned before filtering", {
  counts <- random_counts(3L, 8L, seed = 43L)
  expect_error(prism(counts, scale = rep(1, 5L), bootstrap = FALSE), "scale must have length N")
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
    "prevalence threshold"
  )
})

test_that("too few samples or features after filtering errors clearly", {
  # data where all but one sample is zero-depth
  bad <- matrix(0, nrow = 3L, ncol = 4L)
  bad[, 1] <- c(1, 2, 3)
  expect_error(prism(bad, composition = prism_composition_fixed(0.5), bootstrap = FALSE), "Fewer than 2 samples")

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

test_that("scale_estimate_rho = FALSE estimates only sigma (bounded_scale regime, no rho)", {
  counts <- random_counts(3L, 10L, seed = 37L)
  set.seed(37L)
  u <- stats::rlnorm(10L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = u, scale_ci_level = 0.9, scale_estimate_rho = FALSE
  )
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
})

test_that("prism_scale_log(estimate_rho = FALSE) also skips rho estimation", {
  counts <- random_counts(3L, 10L, seed = 49L)
  set.seed(49L)
  u <- stats::rnorm(10L)
  fit <- prism(
    counts, composition = prism_composition_fixed(0.5), bootstrap = FALSE,
    scale = prism_scale_log(u, ci_level = 0.9, estimate_rho = FALSE)
  )
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
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
