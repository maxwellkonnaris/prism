test_that("prism_scale_diagnostics validates inputs", {
  expect_error(prism_scale_diagnostics(1), "at least 2 finite")
  expect_error(prism_scale_diagnostics(c(1, NA)), "at least 2 finite")
  expect_error(prism_scale_diagnostics(1:5, raw_values = 1:4), "same length")
  expect_error(prism_scale_diagnostics(1:5, group = c("a", "b")), "same length")
})

test_that("log_base check detects a base-10-vs-natural-log mismatch", {
  set.seed(1)
  raw <- exp(rnorm(200, 10, 1))
  log10_mislabeled_as_natural <- log10(raw)  # should be log(raw), i.e. off by a factor of log(10)
  d <- prism_scale_diagnostics(log10_mislabeled_as_natural, raw_values = raw)
  expect_true(d$log_base$flag)
  expect_equal(d$log_base$ratio_median, 1 / log(10), tolerance = 1e-6)
  expect_identical(d$log_base$implied_base, "log10")
})

test_that("log_base check passes when log_values really is log(raw_values)", {
  set.seed(1)
  raw <- exp(rnorm(200, 10, 1))
  d <- prism_scale_diagnostics(log(raw), raw_values = raw)
  expect_false(d$log_base$flag)
  expect_equal(d$log_base$ratio_median, 1, tolerance = 1e-6)
})

test_that("log_base check flags an unrelated quantity even though it is a valid log10", {
  set.seed(1)
  ct <- runif(100, 15, 35)                # a Ct cycle number, not a copy number
  raw_copy_number <- exp(rnorm(100, 20, 2))  # the real underlying total, unrelated to ct
  d <- prism_scale_diagnostics(log10(ct), raw_values = raw_copy_number)
  expect_true(d$log_base$flag)
  expect_match(d$log_base$detail, "not a fixed-base rescaling")
})

test_that("robust_sd flags a single dominating outlier", {
  set.seed(1)
  v <- rnorm(100, 10, 0.3)
  v[1] <- 40
  d <- prism_scale_diagnostics(v)
  expect_true(d$robust_sd$flag)
  expect_gt(d$robust_sd$ratio, 1.3)
})

test_that("robust_sd does not flag well-behaved Gaussian data", {
  set.seed(1)
  v <- rnorm(500, 10, 0.3)
  d <- prism_scale_diagnostics(v)
  expect_false(d$robust_sd$flag)
})

test_that("fold_range flags an implausible span and respects the threshold argument", {
  v <- c(rnorm(50, 10, 0.2), 40)  # one point ~ e^30-fold away from the rest
  d <- prism_scale_diagnostics(v, max_fold_range = 1e4)
  expect_true(d$fold_range$flag)
  d2 <- prism_scale_diagnostics(v, max_fold_range = 1e100)
  expect_false(d2$fold_range$flag)
})

test_that("zero_lod flags non-positive raw values and a repeated floor", {
  raw <- c(rep(0.001, 10), rexp(90, 1) + 1)
  d <- prism_scale_diagnostics(log(pmax(raw, 1e-6)), raw_values = raw)
  expect_true(d$zero_lod$flag)
  expect_equal(d$zero_lod$n_at_floor, 10L)

  raw2 <- c(-1, rexp(90, 1) + 1)
  d2 <- prism_scale_diagnostics(log(pmax(raw2, 1e-6)), raw_values = raw2)
  expect_true(d2$zero_lod$flag)
  expect_identical(d2$zero_lod$n_nonpositive, 1L)
})

test_that("zero_lod and log_base are NULL when raw_values is not supplied", {
  d <- prism_scale_diagnostics(rnorm(50))
  expect_null(d$log_base)
  expect_null(d$zero_lod)
})

test_that("within_between decomposes a known two-group design", {
  set.seed(1)
  g1 <- rnorm(200, 0, 0.1)   # small within-group noise
  g2 <- rnorm(200, 20, 0.1)  # large between-group offset, same small within-group noise
  v <- c(g1, g2)
  grp <- rep(c("a", "b"), each = 200)
  d <- prism_scale_diagnostics(v, group = grp)
  expect_true(d$within_between$flag)
  expect_gt(d$within_between$icc, 0.99)
  expect_lt(d$within_between$sigma_within, 0.2)
})

test_that("within_between is NULL when group is not supplied", {
  d <- prism_scale_diagnostics(rnorm(50))
  expect_null(d$within_between)
})

test_that("print.prism_scale_diagnostics runs without error", {
  d <- prism_scale_diagnostics(rnorm(50), raw_values = exp(rnorm(50)), group = rep(c("a", "b"), 25))
  expect_output(print(d), "log_base")
  expect_output(print(d), "within_between")
})
