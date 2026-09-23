#' Estimate the memory cost of holding all draws in memory
#'
#' Three D x D x S arrays (lower, upper, rel bounds), 8 bytes per double.
#' @keywords internal
.estimate_draw_bytes <- function(D, S) {
  3 * D * D * S * 8
}

.estimate_parallel_result_bytes <- function(D) {
  3 * D * D * 8
}

.parallel_result_memory_budget <- function(memory_limit_bytes) {
  min(memory_limit_bytes, 256 * 1024^2)
}

.resolve_parallel_workers <- function(D, S, workers, memory_limit_bytes) {
  bytes_per_result <- .estimate_parallel_result_bytes(D)
  max_by_memory <- max(1, floor(.parallel_result_memory_budget(memory_limit_bytes) / bytes_per_result))
  as.integer(min(S, workers, max_by_memory))
}

.resolve_parallel_batch_size <- function(D, S, workers, memory_limit_bytes) {
  if (workers <= 1L || S <= 1L) return(1L)
  bytes_per_result <- .estimate_parallel_result_bytes(D)
  memory_budget <- .parallel_result_memory_budget(memory_limit_bytes)
  max_by_memory <- max(1, floor(memory_budget / bytes_per_result))
  as.integer(min(S, workers * 8L, max_by_memory))
}

#' Decide whether draws should be streamed to disk instead of held in memory
#'
#' `stream = TRUE` always streams. `stream = FALSE` (the default) streams
#' anyway when the estimated in-memory cost exceeds `memory_limit_bytes`,
#' so a large problem does not silently exhaust memory; raise
#' `memory_limit_bytes` to opt out of that safety margin.
#'
#' @keywords internal
.resolve_stream <- function(stream, D, S, memory_limit_bytes = 2e9) {
  if (!is.logical(stream) || length(stream) != 1L || is.na(stream)) {
    stop("stream must be TRUE or FALSE.", call. = FALSE)
  }
  memory_limit_bytes <- .assert_scalar_finite(memory_limit_bytes, "stream_memory_limit", lower = 0, lower_inclusive = FALSE)
  estimated_bytes <- .estimate_draw_bytes(D, S)
  list(active = stream || estimated_bytes > memory_limit_bytes, estimated_bytes = estimated_bytes)
}

#' Open temporary binary files for streamed lower/upper/rel draw storage
#' @keywords internal
.stream_open <- function() {
  files <- list(lower = tempfile(fileext = ".bin"), upper = tempfile(fileext = ".bin"), rel = tempfile(fileext = ".bin"))
  cons <- lapply(files, file, open = "wb")
  list(files = files, cons = cons)
}

#' Append one draw's flattened pair values to the stream
#' @keywords internal
.stream_write_draw <- function(handle, pair_index, lower, upper, rel) {
  ij <- cbind(pair_index[, 1], pair_index[, 2])
  writeBin(as.double(lower[ij]), handle$cons$lower, size = 8L, endian = "little")
  writeBin(as.double(upper[ij]), handle$cons$upper, size = 8L, endian = "little")
  writeBin(as.double(rel[ij]), handle$cons$rel, size = 8L, endian = "little")
  invisible(NULL)
}

#' Close the streamed write connections
#' @keywords internal
.stream_close <- function(handle) {
  for (con in handle$cons) close(con)
  invisible(NULL)
}

#' Read a contiguous block of pairs across every draw from a stream file
#'
#' File layout is draw-major: draw `s`'s `n_pairs` values are contiguous,
#' draws are written in generation order. Reading one pair-block requires
#' one seek + read per draw, bounding peak memory to
#' `(pair_end - pair_start + 1) x S` doubles regardless of `n_pairs`.
#' @keywords internal
.stream_read_block <- function(path, n_pairs, S, pair_start, pair_end) {
  width <- pair_end - pair_start + 1L
  con <- file(path, open = "rb")
  on.exit(close(con))
  out <- matrix(NA_real_, nrow = width, ncol = S)
  for (s in seq_len(S)) {
    seek(con, where = ((s - 1L) * n_pairs + (pair_start - 1L)) * 8, origin = "start")
    out[, s] <- readBin(con, what = "double", n = width, size = 8L, endian = "little")
  }
  out
}
