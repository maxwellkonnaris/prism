set.seed(20260711)
D <- 40L
N <- 80L
S <- 40L
counts <- matrix(stats::rpois(D * N, lambda = 20), nrow = D, ncol = N)

profile_call <- function(label, stream) {
  path <- tempfile(fileext = ".mem")
  Rprofmem(path, threshold = 1000)
  elapsed <- system.time({
    fit <- prism::prism(
      counts,
      composition = prism::prism_composition_dirichlet(0.5),
      bootstrap = TRUE,
      S = S,
      sigma_L = 0.2,
      sigma_U = 0.8,
      seed = 1,
      stream = stream
    )
  })[["elapsed"]]
  Rprofmem(NULL)
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  unlink(path)
  data.frame(
    mode = label,
    elapsed_seconds = elapsed,
    largest_profiled_allocation_mb = max(bytes, na.rm = TRUE) / 1024^2,
    result_size_mb = as.numeric(utils::object.size(fit)) / 1024^2,
    stringsAsFactors = FALSE
  )
}

result <- rbind(
  profile_call("full", FALSE),
  profile_call("streamed", 5000)
)
print(result, row.names = FALSE)
