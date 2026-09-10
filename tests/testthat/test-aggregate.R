make_draws <- function(D = 3L, S = 20L, seed = 1L) {
  set.seed(seed)
  lower <- array(NA_real_, c(D, D, S))
  upper <- array(NA_real_, c(D, D, S))
  rel <- array(NA_real_, c(D, D, S))
  for (s in seq_len(S)) {
    m <- matrix(stats::rnorm(D * D, sd = 0.3), D, D)
    m <- (m + t(m)) / 2
    lower[, , s] <- m - 0.2
    upper[, , s] <- m + 0.2
    rel[, , s] <- m
  }
  list(storage = "memory", D = D, S_eff = S, lower = lower, upper = upper, rel = rel)
}

make_stream_draws <- function(D = 3L, S = 20L, seed = 1L) {
  mem <- make_draws(D, S, seed)
  pair_index <- .upper_pairs(D)
  handle <- .stream_open()
  for (s in seq_len(S)) {
    .stream_write_draw(handle, pair_index, mem$lower[, , s], mem$upper[, , s], mem$rel[, , s])
  }
  .stream_close(handle)
  list(storage = "stream", D = D, S_eff = S, files = handle$files, pair_index = pair_index)
}

test_that("pair_draws returns a pairs x draws matrix for S == 1 and S > 1", {
  arr1 <- array(c(1, 2, 3, 4), c(2, 2, 1))
  out1 <- .pair_draws(arr1, i = c(1, 1, 2), j = c(1, 2, 2))
  expect_equal(dim(out1), c(3L, 1L))
  expect_equal(out1[, 1], c(arr1[1, 1, 1], arr1[1, 2, 1], arr1[2, 2, 1]))

  arr3 <- array(seq_len(2 * 2 * 3), c(2, 2, 3))
  out3 <- .pair_draws(arr3, i = c(1, 1, 2), j = c(1, 2, 2))
  expect_equal(dim(out3), c(3L, 3L))
  for (s in 1:3) expect_equal(out3[, s], c(arr3[1, 1, s], arr3[1, 2, s], arr3[2, 2, s]))
})

test_that("ci_lower/ci_upper are exactly the empirical quantiles of the draws", {
  draws <- make_draws(D = 3L, S = 25L, seed = 1L)
  out <- .prism_aggregate(draws, ci_alpha = 0.1)
  for (i in 1:3) for (j in 1:3) {
    expect_equal(
      out$ci_lower[i, j],
      unname(stats::quantile(draws$lower[i, j, ], probs = 0.05)),
      tolerance = 1e-12
    )
    expect_equal(
      out$ci_upper[i, j],
      unname(stats::quantile(draws$upper[i, j, ], probs = 0.95)),
      tolerance = 1e-12
    )
  }
})

test_that("pairwise holds only off-diagonal pairs; pairwise_rel holds i <= j including diagonal", {
  draws <- make_draws(D = 4L, S = 10L, seed = 2L)
  out <- .prism_aggregate(draws, feature_names = letters[1:4])
  expect_identical(nrow(out$pairwise), 6L)
  expect_true(all(out$pairwise$i < out$pairwise$j))
  expect_identical(nrow(out$pairwise_rel), 10L)
  expect_true(all(out$pairwise_rel$i <= out$pairwise_rel$j))
  expect_identical(out$pairwise$feature_i, letters[out$pairwise$i])
})

test_that("rel_hat is the mean of the rel draws and rel_ci_width is upper - lower", {
  draws <- make_draws(D = 3L, S = 12L, seed = 3L)
  out <- .prism_aggregate(draws)
  for (k in seq_len(nrow(out$pairwise_rel))) {
    i <- out$pairwise_rel$i[k]; j <- out$pairwise_rel$j[k]
    expect_equal(out$pairwise_rel$rel_hat[k], mean(draws$rel[i, j, ]), tolerance = 1e-12)
  }
  expect_equal(out$pairwise_rel$rel_ci_width, out$pairwise_rel$rel_ci_upper - out$pairwise_rel$rel_ci_lower)
})

test_that("a single deterministic draw gives NA p-values and the draw itself as the interval", {
  D <- 3L
  m <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), D, D)
  draws <- list(
    storage = "memory", D = D, S_eff = 1L,
    lower = array(m - 0.1, c(D, D, 1)),
    upper = array(m + 0.1, c(D, D, 1)),
    rel = array(m, c(D, D, 1))
  )
  out <- .prism_aggregate(draws)
  expect_true(all(is.na(out$pairwise$p_value)))
  expect_true(all(is.na(out$pairwise$q_value)))
  expect_equal(out$ci_lower, m - 0.1, ignore_attr = TRUE)
  expect_equal(out$ci_upper, m + 0.1, ignore_attr = TRUE)
  expect_equal(out$pairwise_rel$rel_hat, m[cbind(out$pairwise_rel$i, out$pairwise_rel$j)])
})

test_that("q_value is exactly BH-adjusted p_value", {
  draws <- make_draws(D = 5L, S = 30L, seed = 4L)
  out <- .prism_aggregate(draws, delta = 0.05, p_adjust_method = "BH")
  expect_equal(out$pairwise$q_value, stats::p.adjust(out$pairwise$p_value, method = "BH"))

  out_holm <- .prism_aggregate(draws, delta = 0.05, p_adjust_method = "holm")
  expect_equal(out_holm$pairwise$q_value, stats::p.adjust(out_holm$pairwise$p_value, method = "holm"))
})

test_that("p-values are small when draws consistently exclude the null region", {
  D <- 2L
  S <- 50L
  # every draw's lower bound is well above delta for pair (1,2)
  m <- array(0, c(D, D, S))
  lower <- array(0.5, c(D, D, S))
  upper <- array(0.9, c(D, D, S))
  draws <- list(storage = "memory", D = D, S_eff = S, lower = lower, upper = upper, rel = m)
  out <- .prism_aggregate(draws, delta = 0.1)
  expect_true(out$pairwise$p_value[1] < 0.05)
})

test_that("input validation rejects bad ci_alpha, delta, and p_adjust_method", {
  draws <- make_draws(D = 2L, S = 5L, seed = 5L)
  expect_error(.prism_aggregate(draws, ci_alpha = 1), "ci_alpha")
  expect_error(.prism_aggregate(draws, ci_alpha = -0.1), "ci_alpha")
  expect_error(.prism_aggregate(draws, delta = -1), "delta")
  expect_error(.prism_aggregate(draws, p_adjust_method = "not_a_method"), "p_adjust_method")
})

test_that("streamed and in-memory aggregation give identical results", {
  mem <- make_draws(D = 5L, S = 15L, seed = 7L)
  strm <- make_stream_draws(D = 5L, S = 15L, seed = 7L)
  out_mem <- .prism_aggregate(mem, ci_alpha = 0.1, delta = 0.05)
  out_strm <- .prism_aggregate(strm, ci_alpha = 0.1, delta = 0.05)
  expect_equal(out_mem$ci_lower, out_strm$ci_lower)
  expect_equal(out_mem$ci_upper, out_strm$ci_upper)
  expect_equal(out_mem$pairwise, out_strm$pairwise)
  expect_equal(out_mem$pairwise_rel, out_strm$pairwise_rel)
  unlink(unlist(strm$files))
})

test_that("streamed aggregation handles a block boundary correctly", {
  strm <- make_stream_draws(D = 8L, S = 6L, seed = 8L)  # 36 pairs
  mem <- make_draws(D = 8L, S = 6L, seed = 8L)
  out_strm <- .stream_pair_stats(strm, .upper_pairs(8L), ci_alpha = 0.05, delta = 0, block_size = 5L)
  out_mem <- .memory_pair_stats(mem, .upper_pairs(8L), ci_alpha = 0.05, delta = 0)
  expect_equal(out_strm, out_mem)
  unlink(unlist(strm$files))
})

test_that("default delta is 0", {
  draws <- make_draws(D = 2L, S = 30L, seed = 9L)
  out <- .prism_aggregate(draws)
  out_delta0 <- .prism_aggregate(draws, delta = 0)
  expect_identical(out$pairwise$p_value, out_delta0$pairwise$p_value)
})

test_that("feature_names default to V1..VD when not supplied", {
  draws <- make_draws(D = 3L, S = 4L, seed = 6L)
  out <- .prism_aggregate(draws)
  expect_identical(colnames(out$ci_lower), c("V1", "V2", "V3"))
})
