`%||%` <- function(x, y) if (is.null(x)) y else x

.assert_scalar_finite <- function(x, name, lower = -Inf, upper = Inf,
                                   lower_inclusive = TRUE, upper_inclusive = TRUE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (lower_inclusive) x >= lower else x > lower
    ok <- ok && if (upper_inclusive) x <= upper else x < upper
  }
  if (!ok) {
    stop(name, " must be one finite number in ",
         if (lower_inclusive) "[" else "(", lower, ", ", upper,
         if (upper_inclusive) "]" else ")", ".", call. = FALSE)
  }
  as.numeric(x)
}

.assert_positive_int <- function(x, name) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 1 && x == round(x)
  if (!ok) stop(name, " must be one positive integer.", call. = FALSE)
  as.integer(x)
}

#' Save/restore the global RNG state
#'
#' Lets a function call `set.seed()` internally for reproducibility
#' without leaking that reseed into the caller's subsequent random draws.
#' `.restore_rng_state()` removes `.Random.seed` entirely when it did not
#' exist before saving, rather than leaving a stray seed behind.
#' @keywords internal
.save_rng_state <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
}

#' Restore the global RNG state saved by `.save_rng_state()`
#' @keywords internal
.restore_rng_state <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old, envir = .GlobalEnv)
  }
}
