test_that("dirichlet constructor validates and stores options", {
  est <- prism_composition_dirichlet(pseudocount = 0.5, concentration = 2)
  expect_s3_class(est, "prism_composition_estimator")
  expect_identical(est$type, "dirichlet_multinomial")
  expect_identical(est$options$pseudocount, 0.5)
  expect_identical(est$options$concentration, 2)

  expect_error(prism_composition_dirichlet(pseudocount = -1), "pseudocount")
  expect_error(prism_composition_dirichlet(concentration = 0), "concentration")
  expect_error(prism_composition_dirichlet(concentration = -1), "concentration")
})

test_that("fixed constructor validates and stores options", {
  est <- prism_composition_fixed(pseudocount = 0.3)
  expect_identical(est$type, "fixed")
  expect_identical(est$options$pseudocount, 0.3)
  expect_error(prism_composition_fixed(pseudocount = -0.1), "pseudocount")
})

test_that("custom estimator requires a name and fit/draw functions", {
  ok <- prism_composition_estimator(
    name = "my_model",
    fit = function(counts, n_draws, seed, verbose) list(counts = counts),
    draw = function(fit, draw) fit$counts
  )
  expect_identical(ok$type, "custom")
  expect_identical(ok$name, "my_model")

  expect_error(prism_composition_estimator("", fit = identity, draw = identity), "name")
  expect_error(prism_composition_estimator("x", fit = 1, draw = identity), "functions")
  expect_error(prism_composition_estimator("x", fit = identity, draw = 1), "functions")
})

test_that("validate_composition_draw closes valid draws and rejects invalid ones", {
  x <- matrix(c(1, 2, 3, 4), nrow = 2L)
  closed <- .validate_composition_draw(x, D = 2L, N = 2L, estimator_name = "test")
  expect_equal(unname(colSums(closed)), c(1, 1))

  expect_error(
    .validate_composition_draw(matrix(1:6, nrow = 2L), D = 2L, N = 2L, estimator_name = "test"),
    "returned an invalid draw"
  )
  expect_error(
    .validate_composition_draw(matrix(c(-1, 2, 3, 4), nrow = 2L), D = 2L, N = 2L, estimator_name = "test"),
    "returned an invalid draw"
  )
  expect_error(
    .validate_composition_draw(
      matrix(c(1e308, 1e308, 1, 1), nrow = 2L), D = 2L, N = 2L, estimator_name = "test"
    ),
    "returned invalid column sums"
  )
  expect_error(
    .validate_composition_draw(matrix(c(NA, 2, 3, 4), nrow = 2L), D = 2L, N = 2L, estimator_name = "test"),
    "returned an invalid draw"
  )
})
