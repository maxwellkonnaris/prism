#' Specify a compositional uncertainty estimator
#'
#' PRISM separates uncertainty in sample compositions from uncertainty in
#' scale. These constructors define how sample compositions are drawn
#' inside the PRISM bootstrap loop.
#'
#' @param pseudocount Non-negative value added to counts before a logarithm
#'   or before forming Dirichlet parameters.
#' @param concentration Positive multiplier for Dirichlet parameters.
#' @param name Short, unique estimator name.
#' @param fit Function with arguments `counts`, `n_draws`, `seed`, `verbose`.
#'   Must return an arbitrary fitted state object.
#' @param draw Function with arguments `fit` and `draw`. Must return a
#'   finite, strictly positive features-by-samples composition matrix.
#'
#' @return A `prism_composition_estimator` specification.
#' @export
prism_composition_dirichlet <- function(pseudocount = 0.5, concentration = 1) {
  .new_composition_estimator(
    type = "dirichlet_multinomial",
    name = "dirichlet_multinomial",
    options = list(
      pseudocount = .assert_scalar_finite(pseudocount, "pseudocount", lower = 0),
      concentration = .assert_scalar_finite(
        concentration, "concentration", lower = 0, lower_inclusive = FALSE
      )
    )
  )
}

#' @rdname prism_composition_dirichlet
#' @export
prism_composition_fixed <- function(pseudocount = 0.5) {
  .new_composition_estimator(
    type = "fixed",
    name = "fixed",
    options = list(
      pseudocount = .assert_scalar_finite(pseudocount, "pseudocount", lower = 0)
    )
  )
}

#' @rdname prism_composition_dirichlet
#' @export
prism_composition_estimator <- function(name, fit, draw) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("name must be one non-empty string.", call. = FALSE)
  }
  if (!is.function(fit) || !is.function(draw)) {
    stop("fit and draw must be functions.", call. = FALSE)
  }
  .new_composition_estimator(
    type = "custom",
    name = name,
    options = list(fit = fit, draw = draw)
  )
}

.new_composition_estimator <- function(type, name, options) {
  structure(
    list(type = type, name = name, options = options),
    class = "prism_composition_estimator"
  )
}

#' Validate a composition draw
#'
#' Every draw (built-in or from a custom plug-in) is checked and reclosed
#' before covariance estimation.
#'
#' @param x Candidate features-by-samples composition matrix.
#' @param D,N Expected dimensions.
#' @param estimator_name Estimator name, for error messages.
#' @return A finite, strictly positive features-by-samples proportion
#'   matrix whose columns sum to one.
#' @keywords internal
.validate_composition_draw <- function(x, D, N, estimator_name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), c(D, N)) || any(!is.finite(x)) || any(x <= 0)) {
    stop(
      "Composition estimator '", estimator_name,
      "' returned an invalid draw; expected a finite, strictly positive ",
      D, " x ", N, " matrix.",
      call. = FALSE
    )
  }
  totals <- colSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Composition estimator '", estimator_name, "' returned invalid column sums.", call. = FALSE)
  }
  sweep(x, 2, totals, "/")
}
