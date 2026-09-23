rel_cov_example <- function() {
  matrix(c(2, 0.3, 0.3, 1), nrow = 2L, dimnames = list(c("a", "b"), c("a", "b")))
}

cvxr_test_solver <- function() {
  available <- CVXR::installed_solvers()
  preferred <- intersect(c("ECOS", "CLARABEL", "SCS"), available)
  if (length(preferred)) preferred[1L] else NA_character_
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

test_that("genuine rho intervals use exact ellipsoid-box support values", {
  skip_if_not_installed("CVXR")
  solver <- cvxr_test_solver()
  skip_if(is.na(solver), "no supported conic solver is installed")

  A <- diag(c(1, 4))
  rho_L <- c(-1, -1)
  rho_U <- c(0.2, 1)
  support <- .rho_box_support(A, rho_L, rho_U, solver = solver)
  expected_M12 <- 0.2 + 2 * sqrt(1 - 0.2^2)
  expect_equal(support$M[1, 2], expected_M12, tolerance = 1e-5)
  expect_equal(support$m[1, 2], -sqrt(5), tolerance = 1e-5)

  sigma <- 0.4
  bounds <- cov_bounds(
    A, sigma_L = sigma, sigma_U = sigma,
    rho_L = rho_L, rho_U = rho_U, solver = solver
  )
  expect_true(bounds$sharp$off_diagonal)
  expect_equal(bounds$lower[1, 2], sigma^2 + sigma * support$m[1, 2], tolerance = 1e-5)
  expect_equal(bounds$upper[1, 2], sigma^2 + sigma * support$M[1, 2], tolerance = 1e-5)
})

test_that("zero upper scale bound returns Sigma_rel without imposing rho feasibility", {
  A <- diag(2)
  bounds <- cov_bounds(
    A, sigma_L = 0, sigma_U = 0,
    rho_L = c(1, 1), rho_U = c(1, 1)
  )
  expect_equal(bounds$lower, A)
  expect_equal(bounds$upper, A)
  expect_true(bounds$sharp$off_diagonal)
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

test_that("a feasible rho_witness certifies feasibility without CVXR", {
  A <- diag(2)
  expect_true(.rho_box_feasible(
    A, rho_L = c(0.1, 0.1), rho_U = c(0.5, 0.5),
    witness = c(0.2, 0.2), solver = "intentionally-unavailable"
  ))

  expect_error(
    cov_bounds(A, sigma_L = 0.1, sigma_U = 0.5, rho_L = c(-0.5, -0.5), rho_U = c(0.5, 0.5), rho_witness = c(1, 1, 1)),
    "rho_witness"
  )
})

test_that("full-rank feasibility does not require CVXR", {
  A <- diag(2)
  expect_true(.rho_box_feasible(A, rho_L = c(0.1, 0.1), rho_U = c(0.5, 0.5)))
})

test_that("rho_L/rho_U accept a scalar, broadcast to every feature", {
  skip_if_not_installed("CVXR")
  solver <- cvxr_test_solver()
  skip_if(is.na(solver), "no supported conic solver is installed")
  A <- rel_cov_example()
  scalar <- cov_bounds(
    A, sigma_L = 0.1, sigma_U = 0.5,
    rho_L = -0.9, rho_U = 0.9, solver = solver
  )
  vector <- cov_bounds(
    A, sigma_L = 0.1, sigma_U = 0.5,
    rho_L = rep(-0.9, 2), rho_U = rep(0.9, 2), solver = solver
  )
  expect_equal(scalar$lower, vector$lower)
  expect_equal(scalar$upper, vector$upper)
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
  solver <- cvxr_test_solver()
  skip_if(is.na(solver), "no supported conic solver is installed")
  A <- matrix(c(1, 0, 1, 0, 1, 1, 1, 1, 2), nrow = 3L)
  expect_equal(min(eigen(A, only.values = TRUE)$values), 0, tolerance = 1e-8)

  bounds <- cov_bounds(
    A, sigma_L = 0.1, sigma_U = 0.5,
    rho_L = rep(0.01, 3), rho_U = rep(0.05, 3), solver = solver
  )
  expect_true(bounds$sharp$off_diagonal)
  expect_true(all(bounds$lower <= bounds$upper))
})

test_that("singular feasibility rejects a box that violates the range condition", {
  skip_if_not_installed("CVXR")
  solver <- cvxr_test_solver()
  skip_if(is.na(solver), "no supported conic solver is installed")
  A <- matrix(1, 2, 2)
  expect_false(.rho_box_feasible(
    A, rho_L = c(0.2, -0.3), rho_U = c(0.3, -0.2), solver = solver
  ))
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
    "requires the optional CVXR package"
  )
})

test_that("print.prism_cov_bounds runs without error", {
  bounds <- cov_bounds(rel_cov_example())
  expect_output(print(bounds), "PRISM covariance bounds")
})
