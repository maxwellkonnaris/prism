test_that("loading prism has no analysis side effects", {
  script <- tempfile(fileext = ".R")
  package_path <- normalizePath(find.package("prism"), mustWork = TRUE)
  installed <- dir.exists(file.path(package_path, "Meta"))
  writeLines(c(
    "local({",
    "  before_objects <- ls(.GlobalEnv, all.names = TRUE)",
    "  before_files <- list.files('.', all.files = TRUE, recursive = TRUE, no.. = TRUE)",
    if (installed) {
      sprintf(
        "  suppressPackageStartupMessages(library(prism, lib.loc = %s))",
        dQuote(dirname(package_path))
      )
    } else {
      sprintf(
        "  suppressPackageStartupMessages(pkgload::load_all(%s, attach = FALSE, export_all = FALSE, helpers = FALSE, quiet = TRUE, warn_conflicts = FALSE))",
        dQuote(package_path)
      )
    },
    "  after_objects <- ls(.GlobalEnv, all.names = TRUE)",
    "  after_files <- list.files('.', all.files = TRUE, recursive = TRUE, no.. = TRUE)",
    "  created_objects <- setdiff(after_objects, before_objects)",
    "  created_files <- setdiff(after_files, before_files)",
    "  if (length(created_objects)) stop('Package load created global objects: ', paste(created_objects, collapse = ', '))",
    "  if (length(created_files)) stop('Package load created files: ', paste(created_files, collapse = ', '))",
    "  cat('PRISM_LOAD_OK\\n')",
    "})"
  ), script)

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = TRUE,
    stderr = TRUE
  )

  expect_identical(attr(output, "status"), NULL)
  expect_identical(output, "PRISM_LOAD_OK")
})
