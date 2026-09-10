test_that("composition none adds alpha before closure and logging", {
  counts <- matrix(
    c(
      10, 20, 30, 40, 50,
      20, 18, 16, 14, 12,
      5, 7, 9, 11, 13
    ),
    nrow = 3,
    byrow = TRUE
  )

  fit <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    sigma_L = 0,
    sigma_U = 0,
    verbose = FALSE
  )
  adjusted_counts <- counts + 0.5
  expected_composition <- sweep(adjusted_counts, 2, colSums(adjusted_counts), "/")
  expected_covariance <- stats::cov(t(log(expected_composition)))
  pair_index <- cbind(fit$pairwise_rel$i, fit$pairwise_rel$j)

  expect_equal(fit$pairwise_rel$rel_hat, expected_covariance[pair_index])
  expect_equal(fit$parameters$alpha, 0.5)
})

test_that("composition none uses the pre-log pseudocount for zeros", {
  counts <- matrix(c(10, 0, 5, 8, 3, 7), nrow = 2, byrow = TRUE)

  fit <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    sigma_L = 0,
    sigma_U = 0,
    verbose = FALSE
  )

  expect_true(all(is.finite(fit$pairwise_rel$rel_hat)))
})
