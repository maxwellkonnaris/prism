test_that("estimate_draw_bytes matches the direct formula", {
  expect_equal(.estimate_draw_bytes(D = 10L, S = 100L), 3 * (10 * 11 / 2) * 100 * 8)
})

test_that("resolve_stream: TRUE always activates, FALSE stays off below the limit", {
  small <- .resolve_stream(FALSE, D = 5L, S = 10L, memory_limit_bytes = 2e9)
  expect_false(small$active)
  forced <- .resolve_stream(TRUE, D = 5L, S = 10L, memory_limit_bytes = 2e9)
  expect_true(forced$active)
})

test_that("resolve_stream auto-triggers when estimated size exceeds the limit", {
  big <- .resolve_stream(FALSE, D = 300L, S = 3000L, memory_limit_bytes = 1e6)
  expect_true(big$active)
  expect_gt(big$estimated_bytes, 1e6)
})

test_that("resolve_stream validates its inputs", {
  expect_error(.resolve_stream("yes", D = 3L, S = 3L), "stream must be")
  expect_error(.resolve_stream(NA, D = 3L, S = 3L), "stream must be")
  expect_error(.resolve_stream(FALSE, D = 3L, S = 3L, memory_limit_bytes = -1), "stream_memory_limit")
})

test_that("write then read a stream round-trips exactly, including a partial final block", {
  D <- 4L
  S <- 7L
  pair_index <- .upper_pairs(D)
  n_pairs <- nrow(pair_index)
  set.seed(1L)
  lower_true <- matrix(stats::rnorm(n_pairs * S), n_pairs, S)
  upper_true <- lower_true + 1
  rel_true <- matrix(stats::rnorm(n_pairs * S), n_pairs, S)

  handle <- .stream_open()
  for (s in seq_len(S)) {
    lm <- matrix(0, D, D)
    lm[cbind(pair_index[, 1], pair_index[, 2])] <- lower_true[, s]
    um <- matrix(0, D, D)
    um[cbind(pair_index[, 1], pair_index[, 2])] <- upper_true[, s]
    rm_ <- matrix(0, D, D)
    rm_[cbind(pair_index[, 1], pair_index[, 2])] <- rel_true[, s]
    .stream_write_draw(handle, pair_index, lm, um, rm_)
  }
  .stream_close(handle)
  on.exit(unlink(unlist(handle$files)))

  block1 <- .stream_read_block(handle$files$lower, n_pairs, S, 1L, 3L)
  block2 <- .stream_read_block(handle$files$lower, n_pairs, S, 4L, n_pairs)
  expect_equal(rbind(block1, block2), lower_true)

  upper_block <- .stream_read_block(handle$files$upper, n_pairs, S, 1L, n_pairs)
  expect_equal(upper_block, upper_true)
})
