test_that("the two high-resolution diagnostics produce the requested PDF plots", {
  skip_if_not_installed("ggplot2")
  diagnostics <- prism:::.high_resolution_ci_diagnostics(
    i = c(1, 1),
    j = c(2, 3),
    standard_lower = c(-0.12, -0.3),
    standard_upper = c(0.4, 0.08),
    sharp_lower = c(-0.08, -0.25),
    sharp_upper = c(0.35, 0.04)
  )
  fit <- list(
    pairwise = data.frame(
      i = c(1, 1, 2),
      j = c(2, 3, 3),
      ci_lower = c(-0.08, -0.25, 0.1),
      ci_upper = c(0.35, 0.04, 0.5)
    ),
    diagnostics = list(high_resolution = diagnostics)
  )

  endpoints <- plot_high_resolution_ci_distribution(fit, delta = 0.1)
  intervals <- plot_high_resolution_near_zero_changes(fit, delta = 0.1)
  expect_s3_class(endpoints, "ggplot")
  expect_s3_class(intervals, "ggplot")
  expect_equal(nrow(endpoints$data), 6L)
  expect_equal(nrow(intervals$data), 2L)
  expect_setequal(levels(endpoints$data$endpoint), c("Lower", "Upper"))

  endpoint_pdf <- tempfile(fileext = ".pdf")
  interval_pdf <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(endpoint_pdf, endpoints, width = 5, height = 4, device = "pdf")
  ggplot2::ggsave(interval_pdf, intervals, width = 6, height = 4, device = "pdf")
  expect_true(file.exists(endpoint_pdf) && file.info(endpoint_pdf)$size > 0)
  expect_true(file.exists(interval_pdf) && file.info(interval_pdf)$size > 0)
  expect_identical(readChar(endpoint_pdf, nchars = 4L, useBytes = TRUE), "%PDF")
  expect_identical(readChar(interval_pdf, nchars = 4L, useBytes = TRUE), "%PDF")
})
