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
