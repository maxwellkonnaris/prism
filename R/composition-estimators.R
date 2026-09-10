#' Specify a compositional uncertainty estimator
#'
#' PRISM separates uncertainty in sample compositions from uncertainty in scale.
#' These constructors define how sample compositions are supplied to the PRISM
#' covariance-bound engine.
#'
#' @param pseudocount Non-negative value added to counts before a logarithm or
#'   before forming Dirichlet parameters.
#' @param concentration Positive multiplier for Dirichlet parameters.
#' @param ... Model-specific options retained by an MLN estimator.
#' @param name Short, unique estimator name.
#' @param fit Function with arguments `counts`, `n_draws`, `seed`, and `verbose`.
#'   It must return an arbitrary fitted state object.
#' @param draw Function with arguments `fit` and `draw`. It must return a finite,
#'   strictly positive features-by-samples composition matrix whose columns sum
#'   to one (within numerical tolerance).
#'
#' @return A `prism_composition_estimator` specification.
#' @export
prism_composition_dirichlet <- function(pseudocount = 1, concentration = 1) {
  .new_composition_estimator(
    type = "dirichlet_multinomial",
    name = "dirichlet_multinomial",
    options = list(
      pseudocount = .validate_pseudocount(pseudocount),
      concentration = .validate_positive_scalar(concentration, "concentration")
    )
  )
}

#' @rdname prism_composition_dirichlet
#' @export
prism_composition_fixed <- function(pseudocount = 0.5) {
  .new_composition_estimator(
    type = "fixed",
    name = "fixed",
    options = list(pseudocount = .validate_pseudocount(pseudocount))
  )
}

#' @rdname prism_composition_dirichlet
#' @export
prism_composition_mln <- function(pseudocount = 1, ...) {
  .new_composition_estimator(
    type = "mln",
    name = "multinomial_logistic_normal",
    options = c(
      list(pseudocount = .validate_pseudocount(pseudocount)),
      list(...)
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

.validate_pseudocount <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop("pseudocount must be one finite number >= 0.", call. = FALSE)
  }
  as.numeric(x)
}

.validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(name, " must be one positive finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.resolve_composition_estimator <- function(composition, alpha, concentration) {
  if (inherits(composition, "prism_composition_estimator")) {
    return(composition)
  }
  if (!is.character(composition) || length(composition) < 1L) {
    stop(
      "composition must be a PRISM composition estimator or a supported legacy string.",
      call. = FALSE
    )
  }
  choice <- match.arg(composition, c("Dirichlet", "none", "MLN"))
  if (choice == "Dirichlet") {
    return(prism_composition_dirichlet(
      pseudocount = if (is.null(alpha)) 1 else alpha,
      concentration = if (is.null(concentration)) 1 else concentration
    ))
  }
  if (choice == "none") {
    return(prism_composition_fixed(
      pseudocount = if (is.null(alpha)) 0.5 else alpha
    ))
  }
  prism_composition_mln(pseudocount = if (is.null(alpha)) 1 else alpha)
}

.validate_composition_draw <- function(x, D, N, estimator_name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), c(D, N)) || any(!is.finite(x)) || any(x <= 0)) {
    stop(
      "Composition estimator '", estimator_name,
      "' returned an invalid draw; expected a finite, strictly positive D x N matrix.",
      call. = FALSE
    )
  }
  totals <- colSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Composition estimator returned invalid column sums.", call. = FALSE)
  }
  sweep(x, 2, totals, "/")
}
