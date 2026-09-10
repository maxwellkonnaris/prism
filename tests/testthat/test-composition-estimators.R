test_that("Dirichlet-multinomial is the default estimator", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11, 5, 7, 9), nrow = 3)
  fit <- prism(counts, bootstrap = FALSE, S = 4, sigma_L = 0, sigma_U = 0)

  expect_s3_class(fit, "prism_result")
  expect_identical(fit$parameters$composition_estimator, "dirichlet_multinomial")
  expect_identical(fit$parameters$composition, "Dirichlet")
  expect_named(
    fit$diagnostics,
    c("data", "scale", "composition", "sampling", "high_resolution", "simulation")
  )
})

test_that("paired bootstrap resamples empirical parameter blocks together", {
  set.seed(11)
  n <- 20L
  latent <- matrix(stats::rnorm(4L * n), nrow = 4L)
  counts <- matrix(
    stats::rpois(4L * n, lambda = exp(2 + latent / 3)),
    nrow = 4L
  )
  counts[counts < 1] <- 1
  scale_log <- 0.4 * latent[1, ] - 0.2 * latent[2, ] + stats::rnorm(n, sd = 0.5)
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    S = 20L,
    scale_log = scale_log,
    scale_log_ci_level = 0,
    delta = 0,
    stream = FALSE,
    seed = 31L
  )

  fits <- lapply(
    c("none", "both"),
    function(mode) do.call(prism, c(common, list(bootstrap = mode)))
  )
  names(fits) <- c("none", "both")

  expect_identical(fits$none$parameters$bootstrap_mode, "none")
  expect_identical(fits$both$parameters$bootstrap_mode, "both")
  expect_true(any(fits$both$pairwise_rel$rel_ci_width > 0))
  expect_identical(fits$both$diagnostics$sampling$accepted_draws, 20L)
  expect_equal(
    fits$both$diagnostics$sampling$attempted_draws,
    fits$both$diagnostics$sampling$accepted_draws +
      fits$both$diagnostics$sampling$rejected_attempts
  )
  expect_equal(
    fits$both$diagnostics$sampling$acceptance_fraction +
      fits$both$diagnostics$sampling$rejection_fraction,
    1
  )
  expect_identical(fits$both$diagnostics$sampling$rejected_attempts, 0L)
  expect_identical(fits$both$diagnostics$sampling$draws_requiring_replacement, 0L)
  expect_equal(
    sum(fits$both$diagnostics$sampling$rejection_reason_fractions),
    fits$both$diagnostics$sampling$rejection_fraction
  )
})

test_that("logical bootstrap values remain backward compatible", {
  counts <- matrix(c(10, 12, 8, 9, 15, 11, 7, 6, 5, 8, 9, 10), nrow = 3)
  scale_log <- c(-0.4, -0.1, 0.2, 0.5)
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    S = 8L,
    scale_log = scale_log,
    stream = FALSE,
    seed = 19L
  )

  logical_both <- do.call(prism, c(common, list(bootstrap = TRUE)))
  named_both <- do.call(prism, c(common, list(bootstrap = "both")))
  logical_none <- do.call(prism, c(common, list(bootstrap = FALSE)))
  named_none <- do.call(prism, c(common, list(bootstrap = "none")))

  expect_equal(logical_both$pairwise, named_both$pairwise, tolerance = 1e-12)
  expect_equal(logical_none$pairwise, named_none$pairwise, tolerance = 1e-12)
})

test_that("bootstrap replaces draws with zero-variance relative features", {
  proportions <- rbind(
    taxon_a = c(0.1, 0.1, 0.1, 0.2),
    taxon_b = c(0.3, 0.4, 0.5, 0.3),
    taxon_c = c(0.6, 0.5, 0.4, 0.5)
  )
  colnames(proportions) <- paste0("sample_", seq_len(ncol(proportions)))

  fit <- prism(
    proportions,
    counts_input = "proportions",
    composition = prism_composition_fixed(0),
    bootstrap = "both",
    S = 20L,
    scale_log = c(-0.6, -0.2, 0.3, 0.8),
    scale_log_ci_level = 0.5,
    alpha = 0,
    seed = 20260723L,
    stream = TRUE
  )

  sampling <- fit$diagnostics$sampling
  expect_identical(sampling$accepted_draws, 20L)
  expect_identical(sampling$max_attempts_per_draw, 100L)
  expect_gt(sampling$rejected_attempts, 0L)
  expect_gt(
    unname(sampling$rejection_reasons[["zero_variance_relative_feature"]]),
    0
  )
  expect_equal(
    sampling$attempted_draws,
    sampling$accepted_draws + sampling$rejected_attempts
  )
})

test_that("bootstrap proposal cap is validated", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  expect_error(
    prism(counts, S = 2L, max_attempts_per_draw = 0),
    "max_attempts_per_draw must be one positive integer"
  )
})

test_that("component-only bootstrap modes are unsupported", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  for (mode in c("composition", "scale")) {
    expect_error(
      prism(
        counts,
        composition = prism_composition_fixed(0.5),
        bootstrap = mode,
        S = 4,
        sigma_L = 0.2,
        sigma_U = 0.2
      ),
      "deprecated and unsupported"
    )
  }
})

test_that("custom composition estimator follows the fit/draw contract", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  estimator <- prism_composition_estimator(
    name = "test_plugin",
    fit = function(counts, n_draws, seed, verbose) {
      sweep(counts + 1, 2, colSums(counts + 1), "/")
    },
    draw = function(fit, draw) fit
  )
  fit <- prism(
    counts,
    composition = estimator,
    bootstrap = FALSE,
    S = 3,
    sigma_L = 0,
    sigma_U = 0,
    stream = FALSE
  )

  expect_identical(fit$parameters$composition_estimator, "test_plugin")
  expect_identical(fit$parameters$composition, "custom")
  expect_equal(fit$diagnostics$composition$name, "test_plugin")
})

test_that("invalid custom composition draws fail at the interface", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  estimator <- prism_composition_estimator(
    name = "invalid_plugin",
    fit = function(counts, n_draws, seed, verbose) NULL,
    draw = function(fit, draw) matrix(0, nrow = 2, ncol = 3)
  )

  expect_error(
    prism(
      counts,
      composition = estimator,
      bootstrap = FALSE,
      S = 2,
      sigma_L = 0,
      sigma_U = 0,
      stream = FALSE
    ),
    "returned an invalid draw"
  )
})

test_that("MLN is an optional estimator specification", {
  estimator <- prism_composition_mln(tune_draws = 25L, max_hessian_gb = 2)
  expect_s3_class(estimator, "prism_composition_estimator")
  expect_identical(estimator$type, "mln")
  expect_identical(estimator$options$tune_draws, 25L)
})

test_that("optional MLN estimator fits through its isolated module", {
  skip_if_not_installed("fido")
  counts <- matrix(
    c(
      20, 15, 12, 18, 21, 16, 19, 17,
      10, 14, 18, 11, 9, 13, 12, 15,
      7, 9, 11, 8, 10, 12, 9, 11
    ),
    nrow = 3,
    byrow = TRUE
  )
  estimator <- prism_composition_mln(
    tune_draws = 5L,
    center_scale_grid = 1,
    m_grid = 10L,
    omega_structure_grid = "ALR_I",
    gamma_grid = 1,
    max_hessian_gb = 1
  )
  fit <- prism(
    counts,
    composition = estimator,
    bootstrap = FALSE,
    S = 5,
    sigma_L = 0,
    sigma_U = 0,
    seed = 1
  )

  expect_s3_class(fit, "prism_result")
  expect_identical(
    fit$parameters$composition_estimator,
    "multinomial_logistic_normal"
  )
  expect_true(is.list(fit$diagnostics$composition$model))
})

test_that("zero pseudocount is allowed only when logarithms remain finite", {
  positive <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  expect_no_error(prism(
    positive,
    composition = prism_composition_fixed(0),
    bootstrap = FALSE,
    sigma_L = 0,
    sigma_U = 0
  ))

  with_zero <- positive
  with_zero[1, 1] <- 0
  expect_error(
    prism(
      with_zero,
      composition = prism_composition_fixed(0),
      bootstrap = FALSE,
      sigma_L = 0,
      sigma_U = 0
    ),
    "zero pseudocount is invalid"
  )
})

test_that("matrix scale remains sample-aligned after zero-depth filtering", {
  counts <- matrix(c(10, 12, 0, 0, 15, 11), nrow = 2)
  scale <- matrix(c(10, 11, 20, 21, 30, 31), nrow = 3, byrow = TRUE)

  expect_warning(
    fit <- prism(
      counts,
      composition = prism_composition_fixed(0.5),
      bootstrap = FALSE,
      scale = scale
    ),
    "Dropping 1 sample"
  )
  expect_s3_class(fit, "prism_result")
  expect_equal(fit$diagnostics$data$N, 2)
})

test_that("raw-scale uncertainty is reproducible in full and streamed modes", {
  counts <- matrix(c(10, 12, 8, 9, 15, 11, 7, 6, 5, 8, 9, 10), nrow = 3)
  scale <- matrix(c(10, 11, 12, 13, 20, 21, 18, 19), nrow = 4, byrow = TRUE)
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 5,
    scale = scale,
    scale_uncertainty_draws = 10,
    seed = 12
  )

  full <- do.call(prism, c(common, list(stream = FALSE)))
  streamed <- do.call(prism, c(common, list(stream = 4)))
  expect_equal(streamed$ci_lower, full$ci_lower, tolerance = 1e-12)
  expect_equal(streamed$ci_upper, full$ci_upper, tolerance = 1e-12)
  expect_equal(streamed$pairwise, full$pairwise, tolerance = 1e-12)
  expect_equal(streamed$diagnostics$scale, full$diagnostics$scale, tolerance = 1e-12)
  expect_identical(full$parameters$scale_uncertainty_draws, 10L)
})

test_that("deterministic scale_log intervals respect parameter domains", {
  counts <- matrix(c(10, 12, 8, 9, 15, 11, 7, 6, 5, 8, 9, 10), nrow = 3)
  fit <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    scale_log = c(-0.4, -0.1, 0.2, 0.5),
    scale_log_ci_level = 0.8
  )

  expect_gte(fit$parameters$scale_summary$sigma_L, 0)
  expect_gte(fit$parameters$scale_summary$sigma_U, fit$parameters$scale_summary$sigma_L)
  expect_true(all(fit$parameters$scale_summary$rho_L >= -1))
  expect_true(all(fit$parameters$scale_summary$rho_U <= 1))
})

test_that("scale_log can bootstrap scale SD without estimating rho", {
  counts <- matrix(
    c(3, 1, 2, 4, 5, 2, 1, 3),
    nrow = 2,
    dimnames = list(c("A", "B"), paste0("S", 1:4))
  )
  fit <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 12,
    scale_log = c(-1, -0.2, 0.4, 0.8),
    scale_log_ci_level = 0.8,
    scale_log_rho = FALSE,
    stream = TRUE,
    delta = 0,
    seed = 7
  )
  expect_equal(nrow(fit$pairwise), 1)
  expect_true(all(is.finite(fit$pairwise$ci_lower)))
  expect_true(all(is.finite(fit$pairwise$ci_upper)))
  expect_true(all(fit$pairwise$ci_lower <= fit$pairwise$ci_upper))
})

test_that("prism estimation does not create plot files", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  path <- tempfile("prism_side_effect_test_")
  dir.create(path)
  before <- list.files(path, all.files = TRUE)
  old <- setwd(path)
  on.exit(setwd(old), add = TRUE)

  prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    sigma_L = 0,
    sigma_U = 0
  )
  expect_setequal(list.files(path, all.files = TRUE), before)
})

test_that("deterministic fits can return their single identification region", {
  counts <- matrix(
    c(10, 20, 15, 12, 8, 11),
    nrow = 2,
    dimnames = list(c("A", "B"), paste0("S", 1:3))
  )
  fit <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    sigma_L = 0.2,
    sigma_U = 0.8,
    return_bound_draws = TRUE
  )

  expect_identical(fit$pairwise_bound_draws$draw_id, 1L)
  expect_identical(fit$pairwise_bound_draws$pair_index$taxon_i, "A")
  expect_identical(fit$pairwise_bound_draws$pair_index$taxon_j, "B")
  expect_equal(as.numeric(fit$pairwise_bound_draws$lower), fit$pairwise$ci_lower)
  expect_equal(as.numeric(fit$pairwise_bound_draws$upper), fit$pairwise$ci_upper)
})

test_that("legacy analysis arguments preserve current estimator behavior", {
  counts <- matrix(c(10, 20, 15, 12, 8, 11), nrow = 2)
  legacy <- prism(
    counts,
    composition = "none",
    bootstrap = FALSE,
    alpha = 0.5,
    return_mln_fit = FALSE,
    plot_diagnostics = FALSE,
    sigma_L = 0,
    sigma_U = 0
  )
  current <- prism(
    counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = FALSE,
    sigma_L = 0,
    sigma_U = 0
  )

  expect_equal(legacy$ci_lower, current$ci_lower)
  expect_equal(legacy$ci_upper, current$ci_upper)
  expect_equal(legacy$pairwise_rel, current$pairwise_rel)
  expect_identical(legacy$params, legacy$parameters)
})
