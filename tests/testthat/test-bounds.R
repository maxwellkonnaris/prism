rel_cov_example <- function() {
  matrix(c(2, 0.3, 0.3, 1), nrow = 2L, dimnames = list(c("a", "b"), c("a", "b")))
}

test_that("rho bounds without sigma bounds is an error", {
  A <- rel_cov_example()
  expect_error(cov_bounds(A, rho_L = c(0, 0), rho_U = c(0.5, 0.5)), "Correlation bounds require scale-SD bounds")
})

test_that("sigma_L/sigma_U and rho_L/rho_U must each be supplied together", {
  A <- rel_cov_example()
  expect_error(cov_bounds(A, sigma_L = 0.1), "supplied together")
  expect_error(cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0, 0)), "supplied together")
})

test_that("unbounded scale: finite lower, infinite upper, zero diagonal lower", {
  A <- rel_cov_example()
  bounds <- cov_bounds(A)
  expect_identical(bounds$regime, "unbounded_scale")
  kappa <- sqrt(outer(diag(A), diag(A), "+") + 2 * A)
  expect_equal(bounds$lower, A - 0.25 * kappa^2, ignore_attr = TRUE)
  expect_true(all(is.infinite(bounds$upper)))
  expect_equal(unname(diag(bounds$lower)), c(0, 0))
})

test_that("bounded scale (no rho) matches the closed-form vertex-of-quadratic bound", {
  A <- rel_cov_example()
  sigma_L <- 0.1
  sigma_U <- 0.6
  bounds <- cov_bounds(A, sigma_L = sigma_L, sigma_U = sigma_U)
  expect_identical(bounds$regime, "bounded_scale")

  D <- nrow(A)
  expected_lower <- matrix(NA_real_, D, D)
  expected_upper <- matrix(NA_real_, D, D)
  for (i in 1:D) for (j in 1:D) {
    kappa <- sqrt(max(A[i, i] + 2 * A[i, j] + A[j, j], 0))
    sigma_star <- min(max(kappa / 2, sigma_L), sigma_U)
    expected_lower[i, j] <- A[i, j] - sigma_star * kappa + sigma_star^2
    expected_upper[i, j] <- max(A[i, j] + sigma_L^2 + kappa * sigma_L, A[i, j] + sigma_U^2 + kappa * sigma_U)
  }
  expect_equal(unname(bounds$lower), expected_lower, tolerance = 1e-10)
  expect_equal(unname(bounds$upper), expected_upper, tolerance = 1e-10)
})

test_that("fixed sigma and fixed rho give the exact algebraic covariance identity", {
  A <- rel_cov_example()
  sigma <- 0.4
  rho <- c(0.2, -0.1)
  bounds <- cov_bounds(A, sigma_L = sigma, sigma_U = sigma, rho_L = rho, rho_U = rho)
  s <- sqrt(diag(A))
  expected <- A + sigma^2 + sigma * outer(s * rho, s * rho, "+")
  expect_equal(unname(bounds$lower), unname(expected), tolerance = 1e-10)
  expect_equal(unname(bounds$upper), unname(expected), tolerance = 1e-10)
  expect_true(bounds$sharp$off_diagonal)
})

test_that("genuine rho interval matches the conservative box-support formula, and is not sharp", {
  A <- diag(2)
  sigma_L <- 0.2
  sigma_U <- 0.8
  rho_L <- c(-0.3, -0.1)
  rho_U <- c(0.4, 0.5)
  bounds <- cov_bounds(A, sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U)
  expect_false(bounds$sharp$off_diagonal)

  s <- sqrt(diag(A))
  D <- nrow(A)
  expected_lower <- matrix(NA_real_, D, D)
  expected_upper <- matrix(NA_real_, D, D)
  for (i in 1:D) for (j in 1:D) {
    kappa <- sqrt(max(A[i, i] + 2 * A[i, j] + A[j, j], 0))
    t_max <- min(kappa, s[i] * rho_U[i] + s[j] * rho_U[j])
    t_min <- -min(kappa, -(s[i] * rho_L[i] + s[j] * rho_L[j]))
    sigma_star <- min(max(-0.5 * t_min, sigma_L), sigma_U)
    expected_lower[i, j] <- A[i, j] + sigma_star^2 + t_min * sigma_star
    expected_upper[i, j] <- max(A[i, j] + sigma_L^2 + t_max * sigma_L, A[i, j] + sigma_U^2 + t_max * sigma_U)
  }
  expect_equal(unname(bounds$lower), expected_lower, tolerance = 1e-10)
  expect_equal(unname(bounds$upper), expected_upper, tolerance = 1e-10)
})

test_that("a fixed correlation incompatible with Sigma_rel errors on PSD compatibility", {
  A <- diag(2)
  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(1, 1), rho_U = c(1, 1)),
    "positive semidefinite"
  )
})

test_that("an empty rho box (box does not intersect the ellipsoid) errors", {
  A <- diag(3)
  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = rep(0.6, 3), rho_U = rep(0.7, 3)),
    "Empty E intersection C"
  )
})

test_that("a feasible rho_witness certifies a box not bracketing zero, without CVXR", {
  A <- diag(2)
  bounds <- cov_bounds(
    A, sigma_L = 0.1, sigma_U = 0.5,
    rho_L = c(0.1, 0.1), rho_U = c(0.5, 0.5),
    rho_witness = c(0.2, 0.2), solver = "intentionally-unavailable"
  )
  expect_true(all(bounds$lower <= bounds$upper))

  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(-0.5, -0.5), rho_U = c(0.5, 0.5), rho_witness = c(1, 1, 1)),
    "rho_witness"
  )
})

test_that("a non-bracketing box without a witness falls back to the full-rank optimizer", {
  A <- diag(2)
  bounds <- cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.1, 0.1), rho_U = c(0.5, 0.5))
  expect_true(all(bounds$lower <= bounds$upper))
})

test_that("rho_L/rho_U accept a scalar, broadcast to every feature", {
  A <- rel_cov_example()
  bounds_scalar <- cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = 0.1, rho_U = 0.5)
  bounds_vector <- cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.1, 0.1), rho_U = c(0.5, 0.5))
  expect_equal(bounds_scalar$lower, bounds_vector$lower)
  expect_equal(bounds_scalar$upper, bounds_vector$upper)
})

test_that("Sigma_rel must be square, finite, symmetric, and PSD", {
  expect_error(cov_bounds(matrix(1:6, nrow = 2L)), "square")
  expect_error(cov_bounds(matrix(c(1, NA, 0, 1), nrow = 2L)), "finite")
  expect_error(cov_bounds(matrix(c(1, 0.5, 0, 1), nrow = 2L)), "symmetric")
  expect_error(cov_bounds(matrix(c(1, 2, 2, 1), nrow = 2L)), "positive semidefinite")
})

test_that("Sigma_rel cannot contain a zero marginal variance", {
  A <- matrix(c(1, 0, 0, 0), nrow = 2L)
  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5),
    "zero or numerically zero marginal variance"
  )
})

test_that("rho_L must be <= rho_U and both within [-1, 1]", {
  A <- rel_cov_example()
  expect_error(cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(0.5, 0.5), rho_U = c(0.1, 0.1)), "rho_L must be")
  expect_error(cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(1.1, 0), rho_U = c(1.1, 0)), "rho_L")
})

test_that("sigma_L must be <= sigma_U and both non-negative", {
  A <- rel_cov_example()
  expect_error(cov_bounds(A, sigma_L = 0.5, sigma_U = 0.1), "sigma_L must be")
  expect_error(cov_bounds(A, sigma_L = -0.1, sigma_U = 0.5), "sigma_L")
})

test_that("rank-deficient Sigma_rel uses CVXR to certify a non-zero-bracketing rho box", {
  skip_if_not_installed("CVXR")
  solver <- CVXR::installed_solvers()
  skip_if(length(solver) == 0L, "no CVXR solver installed")
  A <- matrix(c(1, 0, 1, 0, 1, 1, 1, 1, 2), nrow = 3L)
  expect_equal(min(eigen(A, only.values = TRUE)$values), 0, tolerance = 1e-8)

  bounds <- cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = rep(0.01, 3), rho_U = rep(0.05, 3))
  expect_false(bounds$sharp$off_diagonal)
  expect_true(all(bounds$lower <= bounds$upper))
})

test_that("rank-deficient Sigma_rel with an invalid CVXR solver name errors clearly", {
  skip_if_not_installed("CVXR")
  A <- matrix(c(1, 0, 1, 0, 1, 1, 1, 1, 2), nrow = 3L)
  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = rep(0.01, 3), rho_U = rep(0.05, 3), solver = "not-a-solver"),
    "is not installed"
  )
})

test_that("rank-deficient Sigma_rel without CVXR or a witness errors clearly", {
  skip_if(requireNamespace("CVXR", quietly = TRUE), "CVXR is installed")
  A <- matrix(c(1, 0, 1, 0, 1, 1, 1, 1, 2), nrow = 3L)
  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = rep(0.01, 3), rho_U = rep(0.05, 3)),
    "requires the CVXR package"
  )
})

test_that("print.prism_cov_bounds runs without error", {
  bounds <- cov_bounds(rel_cov_example())
  expect_output(print(bounds), "PRISM covariance bounds")
})
