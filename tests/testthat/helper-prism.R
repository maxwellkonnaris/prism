random_counts <- function(D = 5L, N = 8L, seed = 1L, lambda = 20) {
  set.seed(seed)
  matrix(stats::rpois(D * N, lambda = lambda), nrow = D, ncol = N)
}
