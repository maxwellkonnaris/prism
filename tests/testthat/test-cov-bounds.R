test_that("correlation bounds require scale-SD bounds", {
  Sigma_rel <- matrix(c(1, 0.2, 0.2, 4), nrow = 2)

  expect_error(
    cov_bounds(Sigma_rel, rho_L = -0.5, rho_U = 0.5),
    "Correlation bounds require scale-SD bounds",
    fixed = TRUE
  )
  expect_error(
    prism:::cov_bounds_pairs(Sigma_rel, rho_L = -0.5, rho_U = 0.5),
    "Correlation bounds require scale-SD bounds",
    fixed = TRUE
  )
})

test_that("fixed correlations determine diagonal and off-diagonal covariance", {
  Sigma_rel <- matrix(c(1, 0.2, 0.2, 4), nrow = 2)
  sigma <- 0.5
  rho <- c(0.3, -0.4)

  bounds <- cov_bounds(
    Sigma_rel,
    sigma_L = sigma,
    sigma_U = sigma,
    rho_L = rho,
    rho_U = rho
  )

  expected <- Sigma_rel + sigma^2 +
    sigma * outer(sqrt(diag(Sigma_rel)) * rho, sqrt(diag(Sigma_rel)) * rho, "+")

  expect_s3_class(bounds, "prism_cov_bounds")
  expect_identical(bounds$regime, "bounded_scale_and_correlation")
  expect_true(bounds$sharp$off_diagonal)
  expect_equal(bounds$lower, expected)
  expect_equal(bounds$upper, expected)
})

test_that("correlation-aware diagonal bounds optimize the exact variance equation", {
  Sigma_rel <- diag(c(1, 4))
  sigma_L <- 0.2
  sigma_U <- 0.8
  rho_L <- c(-0.5, 0.1)
  rho_U <- c(0.25, 0.6)

  bounds <- cov_bounds(
    Sigma_rel,
    sigma_L = sigma_L,
    sigma_U = sigma_U,
    rho_L = rho_L,
    rho_U = rho_U
  )

  relative_sd <- sqrt(diag(Sigma_rel))
  sigma_star <- pmin(pmax(-relative_sd * rho_L, sigma_L), sigma_U)
  expected_lower <- diag(Sigma_rel) + sigma_star^2 +
    2 * sigma_star * relative_sd * rho_L
  expected_upper <- pmax(
    diag(Sigma_rel) + sigma_L^2 + 2 * sigma_L * relative_sd * rho_U,
    diag(Sigma_rel) + sigma_U^2 + 2 * sigma_U * relative_sd * rho_U
  )

  expect_equal(diag(bounds$lower), expected_lower)
  expect_equal(diag(bounds$upper), expected_upper)

  pairs <- prism:::cov_bounds_pairs(
    Sigma_rel,
    pair_index = cbind(1:2, 1:2),
    sigma_L = sigma_L,
    sigma_U = sigma_U,
    rho_L = rho_L,
    rho_U = rho_U
  )
  expect_equal(pairs$lower, expected_lower)
  expect_equal(pairs$upper, expected_upper)
})

test_that("fixed correlations must be compatible with relative covariance", {
  Sigma_rel <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)

  expect_error(
    cov_bounds(
      Sigma_rel,
      sigma_L = 0.5,
      sigma_U = 0.5,
      rho_L = c(1, -1),
      rho_U = c(1, -1)
    ),
    "augmented covariance matrix is not positive semidefinite",
    fixed = TRUE
  )
})

test_that("prism rejects explicit correlations without explicit or measured scale SD", {
  counts <- matrix(c(10, 20, 30, 20, 15, 10), nrow = 2, byrow = TRUE)

  expect_error(
    prism(
      counts,
      composition = "none",
      bootstrap = FALSE,
      rho_L = -0.5,
      rho_U = 0.5,
      verbose = FALSE
    ),
    "Correlation bounds require scale-SD bounds",
    fixed = TRUE
  )
})

test_that("genuine rho intervals require a globally feasible correlation vector", {
  Sigma_rel <- diag(3)

  expect_error(
    cov_bounds(
      Sigma_rel,
      sigma_L = 0.2,
      sigma_U = 0.8,
      rho_L = rep(0.6, 3),
      rho_U = rep(0.7, 3)
    ),
    "Empty E intersection C",
    fixed = TRUE
  )
  expect_error(
    prism:::cov_bounds_pairs(
      Sigma_rel,
      sigma_L = 0.2,
      sigma_U = 0.8,
      rho_L = rep(0.6, 3),
      rho_U = rep(0.7, 3)
    ),
    "Empty E intersection C",
    fixed = TRUE
  )
})

test_that("a valid witness certifies a rank-deficient rho interval without CVXR", {
  Sigma_rel <- matrix(1, nrow = 3, ncol = 3)
  rho_witness <- rep(0.3, 3)

  bounds <- cov_bounds(
    Sigma_rel,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = rep(0.2, 3),
    rho_U = rep(0.4, 3),
    rho_witness = rho_witness,
    solver = "intentionally-unavailable"
  )

  expect_s3_class(bounds, "prism_cov_bounds")
  expect_true(all(bounds$lower <= bounds$upper))
  expect_error(
    cov_bounds(
      Sigma_rel,
      sigma_L = 0.2,
      sigma_U = 0.8,
      rho_L = rep(0.2, 3),
      rho_U = rep(0.4, 3),
      rho_witness = c(0.3, 0.3)
    ),
    "rho_witness must be a finite length-D numeric vector",
    fixed = TRUE
  )
})

test_that("scale-log correlation intervals contain their empirical witness", {
  scale_log <- seq_len(6)
  log_proportions <- rbind(
    scale_log,
    rev(scale_log),
    c(1, 3, 2, 5, 4, 6)
  )

  bounds <- prism:::.estimate_scale_log_bounds(
    log_proportions,
    scale_log,
    ci_level = 0.8
  )

  expect_true(all(bounds$rho_L <= bounds$rho_witness))
  expect_true(all(bounds$rho_U >= bounds$rho_witness))
  expect_equal(bounds$rho_witness[1:2], c(1, -1), tolerance = 1e-12)
})
