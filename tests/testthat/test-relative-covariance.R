test_that("closure produces columns summing to one", {
  counts <- random_counts(5L, 8L)
  P <- prism_closure(counts, pseudocount = 0.5)
  expect_equal(unname(colSums(P)), rep(1, 8L))
  expect_true(all(P > 0))
  expect_true(all(is.finite(P)))
})

test_that("closure matches the hand-computed formula exactly", {
  counts <- matrix(c(0, 3, 5, 1, 2, 0, 4, 6), nrow = 2L)
  pseudocount <- 0.5
  expected <- sweep(counts + pseudocount, 2, colSums(counts + pseudocount), "/")
  expect_equal(prism_closure(counts, pseudocount), expected)
})

test_that("relative covariance is exactly cov(t(log(closure)))", {
  counts <- random_counts(6L, 10L, seed = 2L)
  pseudocount <- 0.3
  expected <- stats::cov(t(log(sweep(
    counts + pseudocount, 2, colSums(counts + pseudocount), "/"
  ))))
  logP <- prism_log_composition(counts, pseudocount)
  expect_equal(prism_relative_covariance(logP), expected)
})

test_that("zero pseudocount is valid only without zero counts", {
  counts_no_zero <- matrix(1:8, nrow = 2L)
  expect_silent(prism_closure(counts_no_zero, pseudocount = 0))

  counts_zero <- matrix(c(0, 3, 5, 1), nrow = 2L)
  expect_error(prism_closure(counts_zero, pseudocount = 0), "pseudocount = 0")
})

test_that("zero counts with a positive pseudocount stay finite", {
  counts <- matrix(c(0, 3, 5, 1, 0, 2), nrow = 2L)
  logP <- prism_log_composition(counts, pseudocount = 0.5)
  expect_true(all(is.finite(logP)))
})

test_that("counts must be a non-empty numeric matrix", {
  expect_error(prism_closure("not a matrix"), "non-empty numeric matrix")
  expect_error(prism_closure(matrix(numeric(0), nrow = 0)), "non-empty numeric matrix")
})

test_that("an overflowing column total after pseudocount errors clearly", {
  counts <- matrix(c(1e308, 1e308, 1, 1), nrow = 2L)
  expect_error(prism_closure(counts, pseudocount = 1e308), "positive, finite total")
})

test_that("counts must be finite, non-negative, and at least 2x2", {
  expect_error(prism_closure(matrix(c(-1, 2, 3, 4), nrow = 2L)), "non-negative")
  expect_error(prism_closure(matrix(c(Inf, 2, 3, 4), nrow = 2L)), "finite")
  expect_error(prism_closure(matrix(1:4, nrow = 1L)), "at least 2 features")
  expect_error(prism_closure(matrix(1:2, nrow = 2L)), "at least 2 features")
})

test_that("pseudocount must be a finite non-negative scalar", {
  counts <- random_counts(3L, 4L)
  expect_error(prism_closure(counts, pseudocount = -1), "pseudocount")
  expect_error(prism_closure(counts, pseudocount = c(0.5, 0.5)), "pseudocount")
  expect_error(prism_closure(counts, pseudocount = NA_real_), "pseudocount")
})

test_that("relative covariance requires at least 2 samples", {
  logP <- matrix(log(c(0.5, 0.5)), nrow = 2L, ncol = 1L)
  expect_error(prism_relative_covariance(logP), "At least 2 samples")
})

test_that("relative covariance rejects non-finite input", {
  logP <- matrix(c(0, -Inf, -1, -2), nrow = 2L)
  expect_error(prism_relative_covariance(logP), "finite")
})

test_that("relative covariance rejects zero marginal feature variance", {
  logP <- rbind(rep(-1, 4L), c(-2, -1, -3, -2))
  expect_error(
    prism_relative_covariance(logP),
    "zero or numerically zero marginal variance"
  )
})
