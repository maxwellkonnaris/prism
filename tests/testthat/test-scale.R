test_that("prism_scale_bounds validates and stores fixed bounds", {
  b <- prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5, rho_L = -0.2, rho_U = 0.4)
  expect_s3_class(b, "prism_scale_bounds")
  expect_identical(b$sigma_L, 0.1)
  expect_identical(b$rho_U, 0.4)

  expect_error(prism_scale_bounds(rho_L = 0, rho_U = 0.5), "require scale-SD bounds")
  expect_error(prism_scale_bounds(sigma_L = 0.5), "supplied together")
  expect_error(prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5, rho_L = 0.5), "supplied together")
  expect_error(prism_scale_bounds(sigma_L = 0.5, sigma_U = 0.1), "sigma_L must be")
  expect_error(prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5, rho_L = 0.6, rho_U = 0.4), "rho_L")
  expect_error(prism_scale_bounds(sigma_L = 0.1, sigma_U = 0.5, rho_L = 1.5, rho_U = 2), "rho_L")

  none <- prism_scale_bounds()
  expect_null(none$sigma_L)
})

test_that("prism_scale_log validates and stores settings", {
  s <- prism_scale_log(c(1, 2, 3), ci_level = 0.8, estimate_rho = FALSE, lower_zero = FALSE)
  expect_s3_class(s, "prism_scale_log")
  expect_identical(s$values, c(1, 2, 3))
  expect_identical(s$ci_level, 0.8)
  expect_false(s$estimate_rho)
  expect_false(s$lower_zero)

  expect_error(prism_scale_log(c(1, NA, 3)), "finite")
  expect_error(prism_scale_log(c(1, 2, 3), ci_level = 1), "ci_level")
  expect_error(prism_scale_log(c(1, 2, 3), estimate_rho = "yes"), "estimate_rho")
  expect_error(prism_scale_log(c(1, 2, 3), lower_zero = "yes"), "lower_zero")
})

test_that("prism_scale_log and .estimate_scale_log_bounds default lower_zero to TRUE", {
  formal_default <- eval(formals(.estimate_scale_log_bounds)$lower_zero)
  expect_true(formal_default)
  expect_true(prism_scale_log(c(1, 2, 3))$lower_zero)
})

test_that("validate_scale_input normalizes a vector to an N x 1 matrix", {
  m <- .validate_scale_input(c(1, 2, 3), N = 3L)
  expect_equal(dim(m), c(3L, 1L))
  expect_equal(as.numeric(m), c(1, 2, 3))
})

test_that("validate_scale_input accepts an N x R replicate matrix as-is", {
  raw <- matrix(1:6, nrow = 3L, ncol = 2L)
  m <- .validate_scale_input(raw, N = 3L)
  expect_equal(dim(m), c(3L, 2L))
})

test_that("validate_scale_input rejects wrong shape, non-positive, or non-finite values", {
  expect_error(.validate_scale_input(c(1, 2), N = 3L), "N rows")
  expect_error(.validate_scale_input(matrix(1:4, nrow = 2L), N = 3L), "N rows")
  expect_error(.validate_scale_input(c(1, -2, 3), N = 3L), "positive")
  expect_error(.validate_scale_input(c(1, Inf, 3), N = 3L), "finite")
})

test_that("select_scale_log picks the requested replicate and logs it", {
  scale_mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3L, ncol = 2L)
  out <- .select_scale_log(scale_mat, sample_index = c(1, 2, 3), replicate_index = c(2, 1, 2))
  expect_equal(out, log(c(10, 2, 30)))
})

test_that("select_scale_log ignores replicate_index for a single-column matrix", {
  scale_mat <- matrix(c(1, 2, 3), nrow = 3L, ncol = 1L)
  out <- .select_scale_log(scale_mat, sample_index = c(2, 1, 3), replicate_index = NULL)
  expect_equal(out, log(c(2, 1, 3)))
})

test_that("select_scale_log supports bootstrap resampling (repeated indices)", {
  scale_mat <- matrix(c(5, 6, 7), nrow = 3L, ncol = 1L)
  out <- .select_scale_log(scale_mat, sample_index = c(1, 1, 3))
  expect_equal(out, log(c(5, 5, 7)))
})

test_that("scale_log bounds recover exact +/-1 correlation for deterministic monotone data", {
  u <- c(1, 2, 3, 4, 5)
  X <- rbind(u, -u, 2 * u + 1)
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0)
  expect_equal(fit$rho_witness, c(1, -1, 1), tolerance = 1e-8)
})

test_that("scale_log CI always contains its own point-estimate witness", {
  set.seed(1)
  N <- 12L
  u <- stats::rnorm(N)
  X <- rbind(0.6 * u + stats::rnorm(N, sd = 0.3), stats::rnorm(N))
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0.9)
  expect_true(all(fit$rho_L <= fit$rho_witness + 1e-10))
  expect_true(all(fit$rho_U >= fit$rho_witness - 1e-10))
  expect_true(fit$sigma_L <= fit$sigma_U)
  expect_true(all(fit$rho_L >= -1) && all(fit$rho_U <= 1))
})

test_that("ci_level = 0 returns point estimates for both endpoints", {
  set.seed(2)
  N <- 10L
  u <- stats::rnorm(N)
  X <- matrix(stats::rnorm(3L * N), nrow = 3L)
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0)
  expect_equal(fit$sigma_L, fit$sigma_U)
  expect_equal(fit$rho_L, fit$rho_U)
  expect_equal(fit$rho_L, fit$rho_witness)
})

test_that("lower_zero forces the lower sigma endpoint to 0", {
  set.seed(3)
  N <- 10L
  u <- stats::rnorm(N)
  X <- matrix(stats::rnorm(2L * N), nrow = 2L)
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0.8, lower_zero = TRUE)
  expect_identical(fit$sigma_L, 0)
})

test_that("estimate_rho = FALSE skips rho entirely", {
  set.seed(4)
  N <- 10L
  u <- stats::rnorm(N)
  X <- matrix(stats::rnorm(2L * N), nrow = 2L)
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0.8, estimate_rho = FALSE)
  expect_null(fit$rho_L)
  expect_null(fit$rho_U)
})

test_that("estimate_scale_log_bounds rejects mismatched or non-finite input", {
  X <- matrix(stats::rnorm(6), nrow = 2L)
  expect_error(.estimate_scale_log_bounds(X, c(1, 2)), "sample-aligned")
  expect_error(.estimate_scale_log_bounds(X, c(1, NA, 3)), "finite")
})

test_that("with n < 4 samples, the rho interval falls back to the point estimate", {
  X <- matrix(stats::rnorm(2L * 3L), nrow = 2L)
  u <- stats::rnorm(3L)
  fit <- .estimate_scale_log_bounds(X, u, ci_level = 0.9)
  expect_identical(fit$rho_L, fit$rho_witness)
  expect_identical(fit$rho_U, fit$rho_witness)
})

test_that("feature_correlations is zero for a zero-variance shared vector", {
  X <- matrix(stats::rnorm(6), nrow = 2L)
  rho <- .feature_correlations(X, rep(1, 3L))
  expect_equal(rho, c(0, 0))
})
