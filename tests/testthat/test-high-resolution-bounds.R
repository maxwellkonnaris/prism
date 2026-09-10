test_that("high resolution is off by default", {
  Sigma_rel <- diag(c(1, 4))
  standard <- prism:::cov_bounds_pairs(
    Sigma_rel,
    pair_index = matrix(c(1, 2), nrow = 1),
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-1, -1),
    rho_U = c(0.2, 1)
  )

  expect_identical(attr(standard, "resolution"), "standard")
  expect_false(attr(standard, "sharp_off_diagonal"))
})

test_that("high resolution fails clearly when CVXR is unavailable", {
  skip_if(requireNamespace("CVXR", quietly = TRUE), "CVXR is installed")
  Sigma_rel <- diag(c(1, 4))

  expect_error(
    prism:::cov_bounds_pairs(
      Sigma_rel,
      pair_index = matrix(c(1, 2), nrow = 1),
      sigma_L = 0.2,
      sigma_U = 0.8,
      rho_L = c(-1, -1),
      rho_U = c(0.2, 1),
      high_resolution = TRUE
    ),
    "requires the CVXR package",
    fixed = TRUE
  )
})

test_that("high-resolution support solves the ellipsoid-box intersection", {
  skip_if_not_installed("CVXR")
  available <- CVXR::installed_solvers()
  solver <- if ("ECOS" %in% available) "ECOS" else if ("CLARABEL" %in% available) "CLARABEL" else NA_character_
  skip_if(is.na(solver), "No supported conic solver is installed")

  Sigma_rel <- diag(c(1, 4))
  sharp <- prism:::cov_bounds_pairs(
    Sigma_rel,
    pair_index = matrix(c(1, 2), nrow = 1),
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-1, -1),
    rho_U = c(0.2, 1),
    high_resolution = TRUE,
    solver = solver
  )

  expected_t_max <- 0.2 + 2 * sqrt(1 - 0.2^2)
  expect_equal(attr(sharp, "support_t_max"), expected_t_max, tolerance = 1e-5)
  expect_identical(attr(sharp, "resolution"), "high")
  expect_true(attr(sharp, "sharp_off_diagonal"))
})

test_that("CI tightening distributions use conservative width as denominator", {
  diagnostics <- prism:::.high_resolution_ci_diagnostics(
    i = c(1, 1),
    j = c(2, 3),
    standard_lower = c(-1, -2),
    standard_upper = c(3, 2),
    sharp_lower = c(0, -1),
    sharp_upper = c(2, 1)
  )

  expect_equal(diagnostics$pairs$ci_lower_tightening_pct, c(25, 25))
  expect_equal(diagnostics$pairs$ci_upper_tightening_pct, c(25, 25))
  expect_equal(nrow(diagnostics$distribution), 4)
  expect_equal(nrow(diagnostics$endpoint_distribution), 4)
  expect_setequal(diagnostics$distribution$endpoint, c("Lower", "Upper"))
  expect_setequal(diagnostics$endpoint_distribution$endpoint, c("Lower", "Upper"))
  expect_equal(
    diagnostics$endpoint_distribution$ci_endpoint,
    c(0, -1, 2, 1)
  )
  expect_s3_class(plot_high_resolution_ci_distribution(diagnostics, delta = 0.1), "ggplot")
  expect_s3_class(
    plot_high_resolution_near_zero_changes(diagnostics, delta = 0.1),
    "ggplot"
  )
})

test_that("endpoint plot uses every tested pair in a PRISM result", {
  diagnostics <- prism:::.high_resolution_ci_diagnostics(
    i = 1,
    j = 2,
    standard_lower = -0.2,
    standard_upper = 0.3,
    sharp_lower = -0.1,
    sharp_upper = 0.2
  )
  result <- list(
    pairwise = data.frame(
      i = c(1, 1, 2),
      j = c(2, 3, 3),
      ci_lower = c(-0.1, -0.4, 0.2),
      ci_upper = c(0.2, 0.5, 0.8)
    ),
    high_resolution_diagnostics = diagnostics
  )

  plot <- plot_high_resolution_ci_distribution(result, delta = 0.1)
  expect_equal(nrow(plot$data), 6)
  expect_equal(sum(plot$data$endpoint == "Lower"), 3)
  expect_equal(sum(plot$data$endpoint == "Upper"), 3)
})

test_that("prism default does not request high-resolution refinement", {
  counts <- matrix(c(10, 20, 30, 20, 15, 10), nrow = 2, byrow = TRUE)
  fit <- prism(
    counts,
    composition = "none",
    bootstrap = FALSE,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-0.5, -0.5),
    rho_U = c(0.5, 0.5),
    verbose = FALSE
  )

  expect_false(fit$parameters$high_resolution)
  expect_null(fit$diagnostics$high_resolution)
})

test_that("prism deterministic high resolution returns CI tightening distributions", {
  skip_if_not_installed("CVXR")
  skip_if_not("ECOS" %in% CVXR::installed_solvers(), "ECOS is unavailable")
  counts <- rbind(
    c(10, 20, 15, 25, 12, 30, 18, 22),
    c(30, 10, 25, 15, 28, 10, 22, 18)
  )
  common <- list(
    counts = counts,
    composition = "none",
    bootstrap = FALSE,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-1, -1),
    rho_U = c(0.2, 1),
    verbose = FALSE
  )

  standard <- do.call(prism, common)
  sharp <- do.call(prism, c(common, list(
    high_resolution = TRUE,
    high_resolution_policy = "all",
    solver = "ECOS"
  )))

  expect_gte(sharp$pairwise$ci_lower, standard$pairwise$ci_lower - 1e-7)
  expect_lte(sharp$pairwise$ci_upper, standard$pairwise$ci_upper + 1e-7)
  expect_equal(nrow(sharp$diagnostics$high_resolution$distribution), 2)
  expect_equal(nrow(sharp$diagnostics$high_resolution$endpoint_distribution), 2)
  expect_setequal(
    names(sharp$diagnostics$high_resolution$pairs),
    c(
      "i", "j", "standard_lower", "standard_upper",
      "high_resolution_lower", "high_resolution_upper",
      "standard_width", "high_resolution_width",
      "ci_lower_tightening_pct", "ci_upper_tightening_pct"
    )
  )
})

test_that("bootstrap refinement occurs at draw level before full-family BH", {
  skip_if_not_installed("CVXR")
  skip_if_not("ECOS" %in% CVXR::installed_solvers(), "ECOS is unavailable")
  counts <- rbind(
    c(10, 20, 15, 25, 12, 30, 18, 22),
    c(30, 10, 25, 15, 28, 10, 22, 18)
  )
  common <- list(
    counts = counts,
    composition = "none",
    bootstrap = TRUE,
    S = 4,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = c(-1, -1),
    rho_U = c(0.2, 1),
    seed = 42,
    verbose = FALSE
  )

  standard <- do.call(prism, common)
  sharp <- do.call(prism, c(common, list(
    high_resolution = TRUE,
    high_resolution_policy = "all",
    solver = "ECOS"
  )))

  expect_gte(sharp$pairwise$ci_lower, standard$pairwise$ci_lower - 1e-7)
  expect_lte(sharp$pairwise$ci_upper, standard$pairwise$ci_upper + 1e-7)
  expect_equal(sharp$pairwise$q_value, stats::p.adjust(sharp$pairwise$p_value, method = "BH"))
  expect_equal(nrow(sharp$diagnostics$high_resolution$pairs), nrow(sharp$pairwise))
})

test_that("requested high-resolution pairs are validated", {
  counts <- rbind(
    c(10, 20, 15, 25),
    c(30, 10, 25, 15),
    c(12, 18, 14, 22)
  )
  common <- list(
    counts = counts,
    composition = "none",
    bootstrap = FALSE,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = rep(-1, 3),
    rho_U = rep(1, 3)
  )

  expect_error(
    do.call(prism, c(common, list(
      high_resolution_pair_index = matrix(c(1, 2), nrow = 1)
    ))),
    "requires high_resolution=TRUE"
  )
  expect_error(
    do.call(prism, c(common, list(
      high_resolution = TRUE,
      high_resolution_pair_index = matrix(c(1, 1), nrow = 1)
    ))),
    "off-diagonal pairs only"
  )
  expect_error(
    do.call(prism, c(common, list(
      high_resolution = TRUE,
      high_resolution_pair_index = rbind(c(1, 2), c(2, 1))
    ))),
    "duplicate unordered pairs"
  )
})

test_that("requested pairs alone are refined inside streamed stochastic draws", {
  skip_if_not_installed("CVXR")
  skip_if_not("ECOS" %in% CVXR::installed_solvers(), "ECOS is unavailable")
  counts <- rbind(
    c(10, 20, 15, 25, 12, 30, 18, 22),
    c(30, 10, 25, 15, 28, 10, 22, 18),
    c(12, 17, 21, 11, 19, 24, 16, 27)
  )
  common <- list(
    counts = counts,
    composition = prism_composition_fixed(0.5),
    bootstrap = TRUE,
    S = 4,
    sigma_L = 0.2,
    sigma_U = 0.8,
    rho_L = rep(-1, 3),
    rho_U = rep(1, 3),
    seed = 42,
    stream = 2,
    verbose = FALSE
  )
  standard <- do.call(prism, common)
  sharp <- do.call(prism, c(common, list(
    high_resolution = TRUE,
    high_resolution_policy = "all",
    high_resolution_pair_index = matrix(c(1, 2), nrow = 1),
    solver = "ECOS"
  )))

  selected <- standard$pairwise$i == 1L & standard$pairwise$j == 2L
  expect_equal(nrow(sharp$diagnostics$high_resolution$pairs), 1L)
  expect_equal(sharp$diagnostics$high_resolution$pairs[c("i", "j")], data.frame(i = 1L, j = 2L))
  expect_gte(sharp$pairwise$ci_lower[selected], standard$pairwise$ci_lower[selected] - 1e-7)
  expect_lte(sharp$pairwise$ci_upper[selected], standard$pairwise$ci_upper[selected] + 1e-7)
  expect_equal(sharp$pairwise[!selected, ], standard$pairwise[!selected, ], tolerance = 1e-12)
  expect_equal(sharp$pairwise$q_value, stats::p.adjust(sharp$pairwise$p_value, method = "BH"))
  expect_identical(sharp$parameters$high_resolution_selection, "requested_pairs")
  expect_identical(sharp$parameters$high_resolution_requested_pairs, 1L)
})
