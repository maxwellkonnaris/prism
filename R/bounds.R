.validate_cov_bounds_inputs <- function(
  Sigma_rel, sigma_L, sigma_U, rho_L, rho_U,
  solver = "ECOS",
  eig_tol = 1e-10,
  rho_witness = NULL,
  check_rho_feasibility = TRUE
) {
  A <- as.matrix(Sigma_rel)
  if (!is.numeric(A) || length(dim(A)) != 2L || nrow(A) != ncol(A) || nrow(A) < 1L) {
    stop("Sigma_rel must be a non-empty numeric square matrix.", call. = FALSE)
  }
  if (any(!is.finite(A))) {
    stop("Sigma_rel must contain only finite values.", call. = FALSE)
  }

  scale <- max(1, max(abs(A)))
  tolerance <- 100 * .Machine$double.eps * scale
  if (max(abs(A - t(A))) > tolerance) {
    stop("Sigma_rel must be symmetric.", call. = FALSE)
  }
  A <- (A + t(A)) / 2
  if (any(diag(A) < -tolerance)) {
    stop("Sigma_rel cannot have negative diagonal variances.", call. = FALSE)
  }
  diag(A) <- pmax(diag(A), 0)
  eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  psd_tolerance <- 1e-8 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -psd_tolerance) {
    stop("Sigma_rel must be positive semidefinite.", call. = FALSE)
  }

  has_sigma <- !is.null(sigma_L) || !is.null(sigma_U)
  has_rho <- !is.null(rho_L) || !is.null(rho_U)
  if (has_sigma && (is.null(sigma_L) || is.null(sigma_U))) {
    stop("Provide both sigma_L and sigma_U, or neither.", call. = FALSE)
  }
  if (has_rho && (is.null(rho_L) || is.null(rho_U))) {
    stop("Provide both rho_L and rho_U, or neither.", call. = FALSE)
  }
  if (has_rho && !has_sigma) {
    stop(
      "Correlation bounds require scale-SD bounds: provide sigma_L and sigma_U with rho_L and rho_U.",
      call. = FALSE
    )
  }

  D <- nrow(A)
  if (!is.null(rho_witness) &&
      (!is.numeric(rho_witness) || length(rho_witness) != D ||
       any(!is.finite(rho_witness)))) {
    stop("rho_witness must be a finite length-D numeric vector.", call. = FALSE)
  }
  if (!has_sigma) {
    return(list(
      A = A,
      D = D,
      sd_rel = sqrt(diag(A)),
      sigma_L = NULL,
      sigma_U = NULL,
      rho_L = NULL,
      rho_U = NULL,
      regime = "unbounded_scale"
    ))
  }

  sigma_L <- as.numeric(sigma_L)
  sigma_U <- as.numeric(sigma_U)
  if (length(sigma_L) != 1L || length(sigma_U) != 1L ||
      !is.finite(sigma_L) || !is.finite(sigma_U) ||
      sigma_L < 0 || sigma_U < sigma_L) {
    stop("Need scalar finite 0 <= sigma_L <= sigma_U.", call. = FALSE)
  }

  if (!has_rho) {
    return(list(
      A = A,
      D = D,
      sd_rel = sqrt(diag(A)),
      sigma_L = sigma_L,
      sigma_U = sigma_U,
      rho_L = NULL,
      rho_U = NULL,
      regime = "bounded_scale"
    ))
  }

  if (length(rho_L) == 1L) rho_L <- rep(as.numeric(rho_L), D)
  if (length(rho_U) == 1L) rho_U <- rep(as.numeric(rho_U), D)
  rho_L <- as.numeric(rho_L)
  rho_U <- as.numeric(rho_U)
  if (length(rho_L) != D || length(rho_U) != D ||
      any(!is.finite(rho_L)) || any(!is.finite(rho_U))) {
    stop("rho_L/rho_U must be finite scalars or length-D vectors.", call. = FALSE)
  }
  if (any(rho_L < -1) || any(rho_U > 1) || any(rho_L > rho_U)) {
    stop("Need -1 <= rho_L <= rho_U <= 1 elementwise.", call. = FALSE)
  }

  sd_rel <- sqrt(diag(A))
  zero_sd <- sd_rel == 0
  if (any(zero_sd & (rho_L != 0 | rho_U != 0))) {
    stop("Correlation with scale is undefined for a zero-variance relative feature; its rho bounds must both be 0.", call. = FALSE)
  }
  if (all(rho_L == rho_U)) {
    covariance_with_unit_scale <- sd_rel * rho_L
    augmented <- rbind(
      cbind(A, covariance_with_unit_scale),
      c(covariance_with_unit_scale, 1)
    )
    augmented_eigenvalues <- eigen(augmented, symmetric = TRUE, only.values = TRUE)$values
    augmented_tolerance <- 1e-8 * max(1, max(abs(augmented_eigenvalues)))
    if (min(augmented_eigenvalues) < -augmented_tolerance) {
      stop(
        "The supplied fixed correlations are incompatible with Sigma_rel: the augmented covariance matrix is not positive semidefinite.",
        call. = FALSE
      )
    }
  } else if (isTRUE(check_rho_feasibility)) {
    # Genuine interval box: pairwise (2-coordinate) support checks done later
    # in cov_bounds()/cov_bounds_pairs() cannot detect infeasibility caused by
    # coordinates outside a given pair (i, j). Certify the FULL D-dimensional
    # E \u2229 C once, up front, exactly as Lemma 3 / Corollary 5 require.
    # check_rho_feasibility = FALSE lets a caller that already certified this
    # once (cov_bounds_pairs(), before its high-resolution per-pair loop) skip
    # repeating the global solve for every pair.
    .check_rho_box_feasible(
      A, rho_L, rho_U, sd_rel,
      eig_tol = eig_tol, solver = solver, witness = rho_witness
    )
  }

  list(
    A = A,
    D = D,
    sd_rel = sd_rel,
    sigma_L = sigma_L,
    sigma_U = sigma_U,
    rho_L = rho_L,
    rho_U = rho_U,
    regime = "bounded_scale_and_correlation"
  )
}

#' Certify that a rho interval box intersects the SPSD feasibility ellipsoid
#'
#' Checks whether there exists a single correlation vector rho satisfying
#' both `rho_L <= rho <= rho_U` (the box C) and the Lemma 3 SPSD condition
#' (the ellipsoid E), i.e. whether `E \eqn{\cap} C` is nonempty. This is a
#' single D-dimensional (rank(Sigma_rel)-dimensional, via the u = M y
#' reparameterization) convex feasibility problem, solved once regardless of
#' how many (i, j) pairs are later queried. It is intentionally distinct from
#' the per-pair `.sharp_EcapC_pair_bound()` used under `high_resolution =
#' TRUE`: pairwise support optimization over a 2-coordinate projection of C
#' cannot certify infeasibility caused by coordinates outside that pair, so
#' it must not be relied on as a feasibility check.
#'
#' @param solver CVXR solver name, used only for rank-deficient `Sigma_rel`
#'   when no witness resolves feasibility.
#' @param witness Optional finite length-D correlation vector known (or
#'   suspected) to lie in `E \eqn{\cap} C`. If it checks out, feasibility is
#'   certified immediately with no optimization or CVXR call at all. An
#'   infeasible or invalid witness does not itself signal infeasibility -- it
#'   just skips the shortcut and falls through to the full check.
#' @keywords internal
# CVXR API break, same family as the earlier solve()->psolve() rename: in
# CVXR <1.x, psolve() returned a rich `cvxr_result` list with $status/
# $value directly. In CVXR 1.9.2 (confirmed via ?psolve and empirically,
# by reproducing the exact "$ operator is invalid for atomic vectors"
# failure against a local CVXR 1.9.2 install matching the HPC container),
# psolve() returns ONLY the bare numeric objective value -- status must
# come from the separate CVXR::status(problem) accessor instead, which
# isn't exported at all in CVXR 1.0.15 (calling it there would itself
# error). Branches on the ACTUAL shape of what psolve() returned, not a
# hardcoded version check, so it works correctly against either CVXR
# generation without needing to know which one is installed:
# CVXR::status() is only referenced (and only needs to exist) inside the
# branch that's taken when psolve() did NOT return the old rich object.
.psolve_status_value <- function(problem, solver) {
  fit <- CVXR::psolve(problem, solver = solver)
  if (is.list(fit) && !is.null(fit[["status"]])) {
    list(status = fit$status, value = fit$value)
  } else {
    list(status = CVXR::status(problem), value = fit)
  }
}

.check_rho_box_feasible <- function(
  Sigma_rel, rho_L, rho_U, sd_rel,
  eig_tol = 1e-10,
  solver = "ECOS",
  witness = NULL
) {
  D <- nrow(Sigma_rel)
  u_L <- sd_rel * rho_L
  u_U <- sd_rel * rho_U

  infeasible_stop <- function() {
    stop(
      "Empty E intersection C: no correlation vector rho within [rho_L, rho_U] is ",
      "compatible with Sigma_rel's positive semidefiniteness constraint (Lemma 3). ",
      "Widen rho_L/rho_U, or check Sigma_rel.",
      call. = FALSE
    )
  }

  if (all(u_L <= 1e-12) && all(u_U >= -1e-12)) {
    return(invisible(TRUE))
  }

  eig <- eigen(Sigma_rel, symmetric = TRUE)
  keep <- eig$values > eig_tol * max(1, max(abs(eig$values)))
  rank_A <- sum(keep)
  V_r <- eig$vectors[, keep, drop = FALSE]
  lambda_r <- eig$values[keep]

  u_is_feasible <- function(u) {
    if (any(u < u_L - 1e-8) || any(u > u_U + 1e-8)) return(FALSE)
    if (rank_A == 0L) return(all(abs(u) <= 1e-8))
    coord <- as.numeric(crossprod(V_r, u))
    if (rank_A < D) {
      residual <- u - as.numeric(V_r %*% coord)
      if (max(abs(residual)) > 1e-6 * max(1, max(abs(u)))) return(FALSE)
    }
    sum(coord^2 / lambda_r) <= 1 + 1e-8
  }

  # Shortcut 1: u = 0 is always in Range(Sigma_rel) and always satisfies the
  # ellipsoid trivially (0 <= 1). If the box brackets zero -- the common case
  # for a symmetric uncertainty interval around an unknown correlation -- rho
  # = 0 is an immediate, solver-free witness.
  # Shortcut 2: caller-supplied witness (e.g. the empirical rho estimated
  # from scale_log in a given bootstrap draw, which is PSD-compatible with
  # that draw's Sigma_rel by construction). Checking it needs the covariance
  # eigendecomposition above, but no numerical optimizer or CVXR solve.
  if (!is.null(witness)) {
    if (u_is_feasible(sd_rel * witness)) return(invisible(TRUE))
  }

  if (rank_A == 0L) {
    # Sigma_rel is (numerically) zero: its range is {0}, so u = 0 is the only
    # feasible point, regardless of the rho box. Shortcut 1 already returned
    # TRUE if that held, so reaching here means it does not.
    infeasible_stop()
  }

  if (rank_A == D) {
    # Full-rank Sigma_rel (the common case): E = {u : u' Sigma_rel^-1 u <= 1}
    # is a box-constrained convex quadratic feasibility problem directly in
    # u-space, solvable exactly with base R (no CVXR dependency needed) via
    # min_u u' Sigma_rel^-1 u  s.t.  u_L <= u <= u_U.
    Sigma_inv <- eig$vectors %*% (t(eig$vectors) / eig$values)
    Sigma_inv <- (Sigma_inv + t(Sigma_inv)) / 2
    objective <- function(u) as.numeric(t(u) %*% Sigma_inv %*% u)
    gradient  <- function(u) as.numeric(2 * (Sigma_inv %*% u))
    u0 <- pmin(pmax(0, u_L), u_U)
    fit <- stats::optim(
      par = u0, fn = objective, gr = gradient, method = "L-BFGS-B",
      lower = u_L, upper = u_U,
      control = list(factr = 1e7, pgtol = 1e-14, maxit = 2000L)
    )
    if (!is.finite(fit$value) || fit$convergence != 0L) {
      stop(
        "E intersection C feasibility solve did not converge (optim status ",
        fit$convergence, ": ", fit$message, "). Bounds near this boundary ",
        "cannot be trusted without a converged solve.",
        call. = FALSE
      )
    }
    if (fit$value > 1 + 1e-6) infeasible_stop()
    return(invisible(TRUE))
  }

  # Rank-deficient Sigma_rel: u is additionally constrained to Range(Sigma_rel),
  # which the box-constrained optim() above cannot enforce directly, and
  # neither shortcut above resolved it. Certify exactly via the same
  # u = M y parameterization used by the high-resolution solver. Rank
  # deficiency is common (more taxa than samples), so failing to certify is
  # an error, not a warning -- silently proceeding would defeat the purpose
  # of this check.
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop(
      "Sigma_rel is rank-deficient, so certifying that rho_L/rho_U is compatible with ",
      "its correlation structure (E intersection C) requires the CVXR package (no ",
      "witness resolved it and the box does not bracket zero). Install CVXR, or pass ",
      "rho_witness with a known-feasible correlation vector, to proceed.",
      call. = FALSE
    )
  }
  available_solvers <- CVXR::installed_solvers()
  if (!(solver %in% available_solvers)) {
    stop(
      "Requested CVXR solver '", solver, "' is unavailable for the E intersection C ",
      "feasibility check. Installed solvers: ", paste(available_solvers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  M <- sweep(V_r, 2, sqrt(lambda_r), `*`)

  y <- CVXR::Variable(ncol(M))
  u <- M %*% y
  problem <- CVXR::Problem(
    CVXR::Minimize(CVXR::sum_squares(y)),
    constraints = list(u >= u_L, u <= u_U)
  )
  solved <- .psolve_status_value(problem, solver)
  fit_value <- solved$value
  fit_status <- solved$status

  infeasible <- isTRUE(fit_status %in% c("infeasible", "infeasible_inaccurate")) ||
    (fit_status %in% c("optimal", "optimal_inaccurate") &&
      length(fit_value) == 1L && is.finite(fit_value) && fit_value > 1 + 1e-8)

  if (infeasible) infeasible_stop()
  if (!(fit_status %in% c("optimal", "optimal_inaccurate")) ||
      length(fit_value) != 1L || !is.finite(fit_value)) {
    stop("E intersection C feasibility solve failed: ", fit_status, call. = FALSE)
  }

  invisible(TRUE)
}

.diagonal_cov_bounds <- function(relative_variance, relative_sd, sigma_L, sigma_U, rho_L = NULL, rho_U = NULL) {
  if (is.null(rho_L)) {
    sigma_star <- pmin(pmax(relative_sd, sigma_L), sigma_U)
    return(list(
      lower = (sigma_star - relative_sd)^2,
      upper = (sigma_U + relative_sd)^2
    ))
  }

  sigma_star <- pmin(pmax(-relative_sd * rho_L, sigma_L), sigma_U)
  lower <- relative_variance + sigma_star^2 + 2 * sigma_star * relative_sd * rho_L
  upper <- pmax(
    relative_variance + sigma_L^2 + 2 * sigma_L * relative_sd * rho_U,
    relative_variance + sigma_U^2 + 2 * sigma_U * relative_sd * rho_U
  )
  list(lower = pmax(lower, 0), upper = pmax(upper, 0))
}

.new_prism_cov_bounds <- function(lower, upper, validated) {
  structure(
    list(
      lower = lower,
      upper = upper,
      regime = validated$regime,
      inputs = list(
        sigma_L = validated$sigma_L,
        sigma_U = validated$sigma_U,
        rho_L = validated$rho_L,
        rho_U = validated$rho_U
      ),
      sharp = list(
        diagonal = TRUE,
        off_diagonal = validated$regime != "bounded_scale_and_correlation" ||
          all(validated$rho_L == validated$rho_U)
      )
    ),
    class = "prism_cov_bounds"
  )
}

#' Print covariance bounds
#'
#' @param x A `prism_cov_bounds` object.
#' @param ... Unused.
#' @export
print.prism_cov_bounds <- function(x, ...) {
  cat("PRISM covariance bounds\n")
  cat("  regime: ", x$regime, "\n", sep = "")
  cat("  dimensions: ", nrow(x$lower), " x ", ncol(x$lower), "\n", sep = "")
  cat("  sharp diagonal: yes\n")
  cat("  sharp off-diagonal: ", if (isTRUE(x$sharp$off_diagonal)) "yes" else "no (conservative E-intersection-C bound)", "\n", sep = "")
  invisible(x)
}

.new_prism_cov_bounds_pairs <- function(
  i,
  j,
  lower,
  upper,
  validated,
  resolution = "standard",
  support_t_min = NULL,
  support_t_max = NULL
) {
  out <- data.frame(i = i, j = j, lower = lower, upper = upper)
  attr(out, "regime") <- validated$regime
  attr(out, "sharp_off_diagonal") <- validated$regime != "bounded_scale_and_correlation" ||
    all(validated$rho_L == validated$rho_U)
  attr(out, "resolution") <- resolution
  if (!is.null(support_t_min)) attr(out, "support_t_min") <- support_t_min
  if (!is.null(support_t_max)) attr(out, "support_t_max") <- support_t_max
  class(out) <- c("prism_cov_bounds_pairs", "data.frame")
  out
}

.solve_EcapC_support <- function(
  Sigma_rel,
  rho_L,
  rho_U,
  c_vec,
  solver = "ECOS",
  eig_tol = 1e-10,
  check_feasibility = TRUE
) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("high_resolution=TRUE requires the CVXR package.", call. = FALSE)
  }
  if (!is.character(solver) || length(solver) != 1L || !nzchar(solver)) {
    stop("solver must be a single non-empty string.", call. = FALSE)
  }
  available_solvers <- CVXR::installed_solvers()
  if (!(solver %in% available_solvers)) {
    stop(
      "Requested CVXR solver '", solver, "' is unavailable. Installed solvers: ",
      paste(available_solvers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.numeric(eig_tol) || length(eig_tol) != 1L ||
      !is.finite(eig_tol) || eig_tol <= 0) {
    stop("eig_tol must be a single positive finite number.", call. = FALSE)
  }

  validated <- .validate_cov_bounds_inputs(
    Sigma_rel,
    sigma_L = 0,
    sigma_U = 0,
    rho_L = rho_L,
    rho_U = rho_U,
    solver = solver,
    eig_tol = eig_tol,
    check_rho_feasibility = check_feasibility
  )
  A <- validated$A
  rho_L <- validated$rho_L
  rho_U <- validated$rho_U
  D <- validated$D
  if (!is.numeric(c_vec) || length(c_vec) != D || any(!is.finite(c_vec))) {
    stop("c_vec must be a finite numeric vector of length D.", call. = FALSE)
  }

  sd_rel <- validated$sd_rel
  u_L <- sd_rel * rho_L
  u_U <- sd_rel * rho_U

  eig <- eigen(A, symmetric = TRUE)
  keep <- eig$values > eig_tol * max(1, max(abs(eig$values)))
  if (!any(keep)) {
    if (any(u_L > 0 | u_U < 0)) {
      stop("Empty E intersection C: only u=0 is feasible, but the rho box excludes 0.", call. = FALSE)
    }
    return(c(t_min = 0, t_max = 0))
  }

  V_r <- eig$vectors[, keep, drop = FALSE]
  lambda_r <- eig$values[keep]
  M <- sweep(V_r, 2, sqrt(lambda_r), `*`)
  rank_A <- ncol(M)

  solve_one <- function(sign) {
    y <- CVXR::Variable(rank_A)
    u <- M %*% y
    objective <- CVXR::sum_entries(
      CVXR::multiply(as.numeric(sign * c_vec), u)
    )
    problem <- CVXR::Problem(
      CVXR::Maximize(objective),
      constraints = list(
        CVXR::p_norm(y, 2) <= 1,
        u >= u_L,
        u <= u_U
      )
    )
    solved <- .psolve_status_value(problem, solver)
    fit_value <- solved$value
    fit_status <- solved$status
    if (!(fit_status %in% c("optimal", "optimal_inaccurate")) ||
        length(fit_value) != 1L || !is.finite(fit_value)) {
      stop("High-resolution E intersection C solve failed: ", fit_status, call. = FALSE)
    }
    as.numeric(fit_value)
  }

  t_max <- solve_one(1)
  t_min <- -solve_one(-1)
  c(t_min = t_min, t_max = t_max)
}

.sharp_EcapC_pair_bound <- function(
  Sigma_rel,
  i,
  j,
  sigma_L,
  sigma_U,
  rho_L,
  rho_U,
  solver = "ECOS",
  eig_tol = 1e-10,
  check_feasibility = TRUE
) {
  A <- as.matrix(Sigma_rel)
  D <- nrow(A)
  if (length(i) != 1L || length(j) != 1L ||
      !is.finite(i) || !is.finite(j) ||
      i < 1L || i > D || j < 1L || j > D) {
    stop("i and j must be scalar indices in 1:D.", call. = FALSE)
  }
  i <- as.integer(i)
  j <- as.integer(j)
  c_vec <- numeric(D)
  if (i == j) {
    c_vec[i] <- 2
  } else {
    c_vec[i] <- 1
    c_vec[j] <- 1
  }

  support <- .solve_EcapC_support(
    Sigma_rel = A,
    rho_L = rho_L,
    rho_U = rho_U,
    c_vec = c_vec,
    solver = solver,
    eig_tol = eig_tol,
    check_feasibility = check_feasibility
  )
  t_min <- support[["t_min"]]
  t_max <- support[["t_max"]]
  sigma_star <- min(max(-0.5 * t_min, sigma_L), sigma_U)
  lower <- A[i, j] + sigma_star^2 + sigma_star * t_min
  upper <- A[i, j] + max(
    sigma_L^2 + sigma_L * t_max,
    sigma_U^2 + sigma_U * t_max
  )
  c(lower = lower, upper = upper, t_min = t_min, t_max = t_max)
}

.resolve_high_resolution <- function(
  high_resolution,
  high_resolution_policy,
  high_resolution_eps
) {
  if (!is.logical(high_resolution) || length(high_resolution) != 1L || is.na(high_resolution)) {
    stop("high_resolution must be a single TRUE or FALSE.", call. = FALSE)
  }
  high_resolution_policy <- match.arg(
    high_resolution_policy,
    c("borderline", "all", "never")
  )
  if (!is.numeric(high_resolution_eps) || length(high_resolution_eps) != 1L ||
      !is.finite(high_resolution_eps) || high_resolution_eps < 0) {
    stop("high_resolution_eps must be a single finite number >= 0.", call. = FALSE)
  }
  list(
    enabled = isTRUE(high_resolution) && high_resolution_policy != "never",
    policy = if (isTRUE(high_resolution)) high_resolution_policy else "never",
    eps = as.numeric(high_resolution_eps)
  )
}

.select_high_resolution_pairs <- function(lower, upper, delta, policy, eps_fraction) {
  if (policy == "all") return(rep(TRUE, length(lower)))
  if (policy == "never") return(rep(FALSE, length(lower)))
  width <- pmax(upper - lower, 0)
  eps <- eps_fraction * width
  near_positive <- lower <= delta & upper > delta & (delta - lower) <= eps
  near_negative <- upper >= -delta & lower < -delta & (upper + delta) <= eps
  near_positive | near_negative
}

.validate_high_resolution_pair_index <- function(pair_index, D) {
  if (is.null(pair_index)) return(NULL)
  pair_index <- as.matrix(pair_index)
  if (!is.numeric(pair_index) || ncol(pair_index) != 2L || nrow(pair_index) < 1L) {
    stop(
      "high_resolution_pair_index must be a non-empty two-column numeric matrix or data frame.",
      call. = FALSE
    )
  }
  if (any(!is.finite(pair_index)) || any(pair_index != floor(pair_index))) {
    stop("high_resolution_pair_index must contain finite integer indices.", call. = FALSE)
  }
  pair_index <- matrix(as.integer(pair_index), ncol = 2L)
  if (any(pair_index < 1L | pair_index > D)) {
    stop("high_resolution_pair_index contains indices outside 1:D after filtering.", call. = FALSE)
  }
  if (any(pair_index[, 1L] == pair_index[, 2L])) {
    stop("high_resolution_pair_index may contain off-diagonal pairs only.", call. = FALSE)
  }
  pair_index <- cbind(
    pmin(pair_index[, 1L], pair_index[, 2L]),
    pmax(pair_index[, 1L], pair_index[, 2L])
  )
  key <- paste(pair_index[, 1L], pair_index[, 2L], sep = ":")
  if (anyDuplicated(key)) {
    stop("high_resolution_pair_index contains duplicate unordered pairs.", call. = FALSE)
  }
  colnames(pair_index) <- c("i", "j")
  pair_index
}

.select_requested_high_resolution_pairs <- function(
  all_pairs,
  requested_pairs,
  lower,
  upper,
  delta,
  policy,
  eps_fraction
) {
  if (is.null(requested_pairs)) {
    return(.select_high_resolution_pairs(lower, upper, delta, policy, eps_fraction))
  }
  all_key <- paste(all_pairs[, 1L], all_pairs[, 2L], sep = ":")
  requested_key <- paste(requested_pairs[, 1L], requested_pairs[, 2L], sep = ":")
  matched <- match(requested_key, all_key)
  if (anyNA(matched)) {
    stop("A requested high-resolution pair is absent after filtering.", call. = FALSE)
  }
  selected <- rep(FALSE, nrow(all_pairs))
  selected[matched] <- TRUE
  selected
}

.validate_high_resolution_tightening <- function(
  standard_lower,
  standard_upper,
  sharp_lower,
  sharp_upper,
  tolerance = 1e-7
) {
  scale <- pmax(1, abs(standard_lower), abs(standard_upper))
  if (any(sharp_lower < standard_lower - tolerance * scale) ||
      any(sharp_upper > standard_upper + tolerance * scale)) {
    stop(
      "High-resolution bounds were wider than standard conservative bounds beyond solver tolerance.",
      call. = FALSE
    )
  }
  list(
    lower = pmax(sharp_lower, standard_lower),
    upper = pmin(sharp_upper, standard_upper)
  )
}

.high_resolution_ci_diagnostics <- function(i, j, standard_lower, standard_upper, sharp_lower, sharp_upper) {
  width <- standard_upper - standard_lower
  lower_gain <- sharp_lower - standard_lower
  upper_gain <- standard_upper - sharp_upper
  lower_pct <- ifelse(width > 0, 100 * lower_gain / width, NA_real_)
  upper_pct <- ifelse(width > 0, 100 * upper_gain / width, NA_real_)
  pairs <- data.frame(
    i = i,
    j = j,
    standard_lower = standard_lower,
    standard_upper = standard_upper,
    high_resolution_lower = sharp_lower,
    high_resolution_upper = sharp_upper,
    standard_width = width,
    high_resolution_width = sharp_upper - sharp_lower,
    ci_lower_tightening_pct = lower_pct,
    ci_upper_tightening_pct = upper_pct
  )
  distribution <- rbind(
    data.frame(i = i, j = j, endpoint = "Lower", tightening_pct = lower_pct),
    data.frame(i = i, j = j, endpoint = "Upper", tightening_pct = upper_pct)
  )
  endpoint_distribution <- rbind(
    data.frame(i = i, j = j, endpoint = "Lower", ci_endpoint = sharp_lower),
    data.frame(i = i, j = j, endpoint = "Upper", ci_endpoint = sharp_upper)
  )
  list(pairs = pairs, distribution = distribution, endpoint_distribution = endpoint_distribution)
}

#' Plot high-resolution CI endpoint distributions
#'
#' Displays lower and upper CI endpoint distributions across all tested pairs
#' in a PRISM result. When a diagnostics list is supplied directly, only the
#' refined pairs available in that list can be shown.
#'
#' @param x A PRISM result containing `high_resolution_diagnostics`, or the
#'   diagnostics list itself.
#' @param bins Ignored; retained for backward compatibility.
#' @param delta Optional non-negative null-region half-width. When supplied,
#'   dashed horizontal lines mark `-delta` and `delta`.
#'
#' @return A `ggplot` object.
#' @export
plot_high_resolution_ci_distribution <- function(x, bins = 30L, delta = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_high_resolution_ci_distribution() requires ggplot2.", call. = FALSE)
  }
  if (!is.null(delta) &&
      (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0)) {
    stop("delta must be a single finite number >= 0.", call. = FALSE)
  }
  diagnostics <- if (!is.null(x$diagnostics$high_resolution)) {
    x$diagnostics$high_resolution
  } else if (!is.null(x$high_resolution_diagnostics)) {
    x$high_resolution_diagnostics
  } else {
    x
  }
  distribution <- NULL
  if (is.data.frame(x$pairwise) &&
      all(c("i", "j", "ci_lower", "ci_upper") %in% names(x$pairwise))) {
    distribution <- rbind(
      data.frame(
        i = x$pairwise$i,
        j = x$pairwise$j,
        endpoint = "Lower",
        ci_endpoint = x$pairwise$ci_lower
      ),
      data.frame(
        i = x$pairwise$i,
        j = x$pairwise$j,
        endpoint = "Upper",
        ci_endpoint = x$pairwise$ci_upper
      )
    )
  }
  if (is.null(distribution)) distribution <- diagnostics$endpoint_distribution
  if (is.null(distribution) && is.data.frame(diagnostics$pairs)) {
    distribution <- rbind(
      data.frame(
        i = diagnostics$pairs$i,
        j = diagnostics$pairs$j,
        endpoint = "Lower",
        ci_endpoint = diagnostics$pairs$high_resolution_lower
      ),
      data.frame(
        i = diagnostics$pairs$i,
        j = diagnostics$pairs$j,
        endpoint = "Upper",
        ci_endpoint = diagnostics$pairs$high_resolution_upper
      )
    )
  }
  if (!is.data.frame(distribution) ||
      !all(c("endpoint", "ci_endpoint") %in% names(distribution)) ||
      nrow(distribution) == 0L) {
    stop("No CI endpoint distribution is available.", call. = FALSE)
  }
  distribution <- distribution[is.finite(distribution$ci_endpoint), , drop = FALSE]
  if (nrow(distribution) == 0L) {
    stop("CI endpoint values are all non-finite.", call. = FALSE)
  }
  distribution$endpoint <- factor(distribution$endpoint, levels = c("Lower", "Upper"))

  p <- ggplot2::ggplot(
    distribution,
    ggplot2::aes(x = endpoint, y = ci_endpoint, color = endpoint)
  ) +
    ggplot2::geom_boxplot(
      width = 0.42,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
      size = 1.9,
      alpha = 0.78
    ) +
    ggplot2::scale_color_manual(values = c("Lower" = "#C43C39", "Upper" = "#2E5E9E"), guide = "none") +
    ggplot2::labs(
      x = NULL,
      y = "Covariance CI endpoint"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  if (is.null(delta)) {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4)
  } else {
    p <- p +
      ggplot2::geom_hline(
        yintercept = unique(c(-delta, delta)),
        color = "grey35",
        linetype = "dashed",
        linewidth = 0.45
      )
  }
  p
}

#' Plot standard and high-resolution intervals for refined pairs
#'
#' Shows conservative and high-resolution intervals for taxon pairs refined by
#' the high-resolution solver.
#'
#' @param x A PRISM result containing `high_resolution_diagnostics`, or the
#'   diagnostics list itself.
#' @param near_zero_threshold Optional non-negative covariance-distance
#'   threshold used to restrict the plot to pairs whose conservative interval
#'   has an endpoint near zero. Default `NULL` displays all refined pairs.
#' @param delta Optional non-negative null-region half-width. When supplied,
#'   dashed vertical lines mark `-delta` and `delta`.
#'
#' @return A `ggplot` object.
#' @export
plot_high_resolution_near_zero_changes <- function(x, near_zero_threshold = NULL, delta = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_high_resolution_near_zero_changes() requires ggplot2.", call. = FALSE)
  }
  if (!is.null(near_zero_threshold) &&
      (!is.numeric(near_zero_threshold) || length(near_zero_threshold) != 1L ||
       !is.finite(near_zero_threshold) || near_zero_threshold < 0)) {
    stop("near_zero_threshold must be a single finite number >= 0.", call. = FALSE)
  }
  if (!is.null(delta) &&
      (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0)) {
    stop("delta must be a single finite number >= 0.", call. = FALSE)
  }
  diagnostics <- if (!is.null(x$diagnostics$high_resolution)) {
    x$diagnostics$high_resolution
  } else if (!is.null(x$high_resolution_diagnostics)) {
    x$high_resolution_diagnostics
  } else {
    x
  }
  pairs <- diagnostics$pairs
  required <- c(
    "i", "j", "standard_lower", "standard_upper",
    "high_resolution_lower", "high_resolution_upper"
  )
  if (!is.data.frame(pairs) || !all(required %in% names(pairs)) || nrow(pairs) == 0L) {
    stop("No high-resolution pair diagnostics are available.", call. = FALSE)
  }
  pairs$near_zero_distance <- pmin(abs(pairs$standard_lower), abs(pairs$standard_upper))
  if (!is.null(near_zero_threshold)) {
    pairs <- pairs[pairs$near_zero_distance <= near_zero_threshold, , drop = FALSE]
    if (nrow(pairs) == 0L) {
      stop("No refined pair has a conservative endpoint within near_zero_threshold of zero.", call. = FALSE)
    }
  }
  pairs$pair <- paste(pairs$i, pairs$j, sep = "-")
  pairs <- pairs[order(pairs$near_zero_distance, pairs$pair), , drop = FALSE]
  pairs$pair <- factor(pairs$pair, levels = rev(pairs$pair))

  p <- ggplot2::ggplot(pairs, ggplot2::aes(y = pair)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = standard_lower, xmax = standard_upper, color = "Conservative"),
      orientation = "y",
      width = 0.26,
      linewidth = 1.15
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = high_resolution_lower,
        xmax = high_resolution_upper,
        color = "High resolution"
      ),
      orientation = "y",
      width = 0.14,
      linewidth = 0.75
    ) +
    ggplot2::scale_color_manual(
      values = c("Conservative" = "#A6A6A6", "High resolution" = "#111111"),
      breaks = c("Conservative", "High resolution")
    ) +
    ggplot2::labs(
      x = "Covariance CI",
      y = "High-resolution refined pair",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "top",
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  if (is.null(delta)) {
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "grey25", linewidth = 0.45)
  } else {
    p <- p +
      ggplot2::geom_vline(
        xintercept = unique(c(-delta, delta)),
        color = "grey35",
        linetype = "dashed",
        linewidth = 0.45
      )
  }
  p
}

#' Bound absolute log-covariance from relative covariance and scale information
#'
#' Selects the identification equations from the supplied information:
#' unbounded scale when no scale information is supplied, bounded scale when
#' only `sigma_L` and `sigma_U` are supplied, and bounded scale with restricted
#' correlations when both sigma and rho bounds are supplied. Correlation bounds
#' without scale-SD bounds are invalid.
#'
#' @param Sigma_rel Finite symmetric positive-semidefinite relative
#'   log-covariance matrix.
#' @param sigma_L,sigma_U Optional finite scalar lower and upper bounds for the
#'   SD of log scale, satisfying `0 <= sigma_L <= sigma_U`.
#' @param rho_L,rho_U Optional finite scalar or length-D lower and upper bounds
#'   for correlations between each relative log component and log scale. Values
#'   must lie in `[-1, 1]` and require sigma bounds. A fixed correlation
#'   (`rho_L == rho_U`) coming from the same composition and scale data is
#'   PSD-compatible by construction; a genuine interval is only meaningful if
#'   some correlation vector in that interval is simultaneously compatible
#'   with `Sigma_rel` (Lemma 3) -- this is verified once, up front, before any
#'   bound is computed, and an incompatible interval is a hard error, not a
#'   silently-returned (and meaningless) bound.
#' @param rho_witness Optional finite length-D vector: a correlation vector
#'   believed to lie in the feasible set (e.g. an empirical rho estimated
#'   alongside `Sigma_rel`). If it checks out, feasibility is certified with
#'   no optimizer or CVXR call. A malformed witness is an error; a well-formed
#'   but infeasible witness falls through to the full feasibility check.
#' @param solver,eig_tol Only used to certify feasibility when `Sigma_rel` is
#'   rank-deficient and no witness (nor the box bracketing zero) already
#'   resolves it; in that case certification requires the CVXR package. This
#'   is unrelated to `cov_bounds_pairs()`'s `high_resolution` argument.
#'
#' @return A `prism_cov_bounds` object. Its `lower` and `upper` elements are
#'   covariance-bound matrices; `regime` records the selected equations;
#'   `inputs` records validated bounds; and `sharp` states whether diagonal and
#'   off-diagonal bounds are sharp. Correlation-interval off-diagonal bounds use
#'   the package's conservative ellipsoid-intersection-box support bound.
#' @export
cov_bounds <- function(
  Sigma_rel, sigma_L = NULL, sigma_U = NULL, rho_L = NULL, rho_U = NULL,
  rho_witness = NULL, solver = "ECOS", eig_tol = 1e-10
) {
  validated <- .validate_cov_bounds_inputs(
    Sigma_rel, sigma_L, sigma_U, rho_L, rho_U,
    solver = solver, eig_tol = eig_tol, rho_witness = rho_witness
  )
  A <- validated$A
  D <- validated$D
  d <- diag(A)
  sigma_rel_sd <- validated$sd_rel

  # kappa_ij^2 = Var(log W_i^|| + log W_j^||) = A_ii + 2*A_ij + A_jj
  kappa2 <- outer(d, d, "+") + 2 * A
  kappa2 <- pmax(kappa2, 0)
  kappa_ellipsoid <- sqrt(kappa2)

  if (validated$regime == "unbounded_scale") {
    # Theorem 1 (unbounded-scale). rho bounds do not change the +Inf upper bound.
    lower <- A - 0.25 * kappa2
    upper <- matrix(Inf, D, D)
    diag(lower) <- 0
    return(.new_prism_cov_bounds(lower, upper, validated))
  }

  sigma_L <- validated$sigma_L
  sigma_U <- validated$sigma_U

  if (validated$regime == "bounded_scale") {
    # Corollary 3 (bounded scale, unbounded rho)
    kappa <- kappa_ellipsoid
    sigma_star <- pmin(pmax(kappa / 2, sigma_L), sigma_U)

    lower <- A - sigma_star * kappa + sigma_star^2
    upper <- A + sigma_U * kappa + sigma_U^2

    diagonal <- .diagonal_cov_bounds(d, sigma_rel_sd, sigma_L, sigma_U)
    diag(lower) <- diagonal$lower
    diag(upper) <- diagonal$upper

    return(.new_prism_cov_bounds(lower, upper, validated))
  }

  # Corollary 4 (bounded scale, bounded rho in E intersection C)
  rho_L <- validated$rho_L
  rho_U <- validated$rho_U

  # Conservative (valid) kappa_sup and kappa_inf using the paper's simple box upper bound.
  # For w = sigma_i e_i + sigma_j e_j (sigma >= 0), the box max is sigma_i*rho_U[i] + sigma_j*rho_U[j].
  tU <- sigma_rel_sd * rho_U
  box_max <- outer(tU, tU, "+")
  kappa_sup <- pmin(kappa_ellipsoid, box_max)

  # For -w, the box max is -(sigma_i*rho_L[i] + sigma_j*rho_L[j]).
  tL <- sigma_rel_sd * rho_L
  box_max_negw <- -outer(tL, tL, "+")
  ub_sup_negw <- pmin(kappa_ellipsoid, box_max_negw)
  kappa_inf <- -ub_sup_negw
  if (any(kappa_inf > kappa_sup + 1e-12)) {
    stop("The supplied rho bounds are incompatible with the relative-covariance geometry for at least one pair.", call. = FALSE)
  }

  # Lower bound: choose sigma closest to vertex at -0.5*kappa_inf, then plug in
  sigma_inf_star <- pmin(pmax(-0.5 * kappa_inf, sigma_L), sigma_U)
  lower <- A + sigma_inf_star^2 + kappa_inf * sigma_inf_star

  # Upper bound: max at the boundaries
  upper <- pmax(
    A + sigma_L^2 + kappa_sup * sigma_L,
    A + sigma_U^2 + kappa_sup * sigma_U
  )

  diagonal <- .diagonal_cov_bounds(d, sigma_rel_sd, sigma_L, sigma_U, rho_L, rho_U)
  diag(lower) <- diagonal$lower
  diag(upper) <- diagonal$upper

  .new_prism_cov_bounds(lower, upper, validated)
}


#' Pairwise/block version of cov_bounds()
#'
#' Computes covariance bounds only for selected (i, j) pairs.
#' This is intended to match cov_bounds(Sigma_rel, ...)$lower/$upper
#' on the requested indices, without constructing full D x D lower/upper outputs.
#'
#' @param Sigma_rel D x D relative covariance matrix.
#' @param pair_index Two-column matrix/data.frame with integer columns i, j.
#'   If NULL, all upper-triangle pairs including diagonal are used.
#' @param sigma_L,sigma_U Optional scalar lower/upper bounds for scale SD.
#' @param rho_L,rho_U Optional scalar or length-D lower/upper bounds for rho.
#'   A genuine interval (`rho_L != rho_U` for some feature) is checked once,
#'   up front, for a nonempty intersection with the Lemma 3 feasibility
#'   ellipsoid; an incompatible interval is a hard error rather than a
#'   silently-returned bound. `high_resolution` below is unrelated to this
#'   check and does not skip it.
#' @param rho_witness Optional finite length-D vector believed to lie in the
#'   feasible set; if it checks out, the feasibility check needs no optimizer
#'   or CVXR call. Passed through to the same check
#'   `cov_bounds()` uses.
#' @param high_resolution Logical. If `TRUE`, use an exact numerical
#'   ellipsoid-box support solve for nonzero-width rho intervals, on top of
#'   (not instead of) the upfront feasibility check above. Default `FALSE`
#'   uses the fast conservative equations. The upfront feasibility check is
#'   performed exactly once regardless of this setting or the number of pairs
#'   requested -- it is not repeated inside the per-pair high-resolution loop.
#' @param solver CVXR solver name used when `high_resolution = TRUE`, and also
#'   for the upfront feasibility check if `Sigma_rel` turns out to be
#'   rank-deficient and no witness resolves it.
#' @param eig_tol Relative eigenvalue tolerance for the high-resolution
#'   range-space parameterization (and for the same parameterization inside
#'   the upfront feasibility check).
#'
#' @return data.frame with columns i, j, lower, upper.
cov_bounds_pairs <- function(
  Sigma_rel,
  pair_index = NULL,
  sigma_L = NULL,
  sigma_U = NULL,
  rho_L = NULL,
  rho_U = NULL,
  rho_witness = NULL,
  high_resolution = FALSE,
  solver = "ECOS",
  eig_tol = 1e-10
) {
  if (!is.logical(high_resolution) || length(high_resolution) != 1L || is.na(high_resolution)) {
    stop("high_resolution must be a single TRUE or FALSE.", call. = FALSE)
  }
  validated <- .validate_cov_bounds_inputs(
    Sigma_rel, sigma_L, sigma_U, rho_L, rho_U,
    solver = solver, eig_tol = eig_tol, rho_witness = rho_witness
  )
  if (isTRUE(high_resolution) && validated$regime != "bounded_scale_and_correlation") {
    warning("high_resolution=TRUE has no effect outside bounded rho interval mode.", call. = FALSE)
    high_resolution <- FALSE
  }
  A <- validated$A
  D <- validated$D
  d <- diag(A)
  sigma_rel_sd <- validated$sd_rel

  if (is.null(pair_index)) {
    pair_index <- which(upper.tri(matrix(0, D, D), diag = TRUE), arr.ind = TRUE)
  } else {
    pair_index <- as.matrix(pair_index)
  }

  if (ncol(pair_index) != 2L) {
    stop("pair_index must have exactly two columns: i and j.", call. = FALSE)
  }

  i <- as.integer(pair_index[, 1])
  j <- as.integer(pair_index[, 2])

  if (length(i) == 0L) {
    return(.new_prism_cov_bounds_pairs(integer(), integer(), numeric(), numeric(), validated))
  }

  if (any(!is.finite(i)) || any(!is.finite(j))) {
    stop("pair_index contains non-finite indices.", call. = FALSE)
  }

  if (any(i < 1L | i > D | j < 1L | j > D)) {
    stop("pair_index contains indices outside 1:D.", call. = FALSE)
  }

  Aij <- A[cbind(i, j)]
  di <- d[i]
  dj <- d[j]
  sdi <- sigma_rel_sd[i]
  sdj <- sigma_rel_sd[j]

  # kappa_ij^2 = A_ii + 2 A_ij + A_jj
  kappa2 <- pmax(di + 2 * Aij + dj, 0)
  kappa <- sqrt(kappa2)

  is_diag <- i == j

  # Case 1: unbounded scale.
  if (validated$regime == "unbounded_scale") {
    lower <- Aij - 0.25 * kappa2
    upper <- rep(Inf, length(Aij))

    lower[is_diag] <- 0

    return(.new_prism_cov_bounds_pairs(i, j, lower, upper, validated))
  }

  sigma_L <- validated$sigma_L
  sigma_U <- validated$sigma_U

  # Case 2: bounded scale, unbounded rho.
  if (validated$regime == "bounded_scale") {
    sigma_star <- pmin(pmax(kappa / 2, sigma_L), sigma_U)

    lower <- Aij - sigma_star * kappa + sigma_star^2
    upper <- Aij + sigma_U * kappa + sigma_U^2

    diagonal <- .diagonal_cov_bounds(di, sdi, sigma_L, sigma_U)
    lower[is_diag] <- diagonal$lower[is_diag]
    upper[is_diag] <- diagonal$upper[is_diag]

    return(.new_prism_cov_bounds_pairs(i, j, lower, upper, validated))
  }

  # Case 3: bounded scale, bounded rho.
  rho_L <- validated$rho_L
  rho_U <- validated$rho_U

  if (isTRUE(high_resolution)) {
    if (all(rho_L == rho_U)) {
      warning("high_resolution=TRUE is unnecessary when rho_L == rho_U; using exact closed-form bounds.", call. = FALSE)
    } else {
      out_lower <- numeric(length(i))
      out_upper <- numeric(length(i))
      out_t_min <- numeric(length(i))
      out_t_max <- numeric(length(i))
      for (k in seq_along(i)) {
        bound <- .sharp_EcapC_pair_bound(
          Sigma_rel = A,
          i = i[k],
          j = j[k],
          sigma_L = sigma_L,
          sigma_U = sigma_U,
          rho_L = rho_L,
          rho_U = rho_U,
          solver = solver,
          eig_tol = eig_tol,
          # The global E \u2229 C feasibility check already ran once, above,
          # in this function's own .validate_cov_bounds_inputs() call. Do not
          # repeat that (expensive) global solve inside every one of the
          # per-pair .solve_EcapC_support() calls below.
          check_feasibility = FALSE
        )
        out_lower[k] <- bound[["lower"]]
        out_upper[k] <- bound[["upper"]]
        out_t_min[k] <- bound[["t_min"]]
        out_t_max[k] <- bound[["t_max"]]
      }
      out <- .new_prism_cov_bounds_pairs(
        i,
        j,
        out_lower,
        out_upper,
        validated,
        resolution = "high",
        support_t_min = out_t_min,
        support_t_max = out_t_max
      )
      attr(out, "sharp_off_diagonal") <- TRUE
      return(out)
    }
  }

  tU <- sigma_rel_sd * rho_U
  box_max <- tU[i] + tU[j]
  kappa_sup <- pmin(kappa, box_max)

  tL <- sigma_rel_sd * rho_L
  box_max_negw <- -(tL[i] + tL[j])
  ub_sup_negw <- pmin(kappa, box_max_negw)
  kappa_inf <- -ub_sup_negw
  if (any(kappa_inf > kappa_sup + 1e-12)) {
    stop("The supplied rho bounds are incompatible with the relative-covariance geometry for at least one requested pair.", call. = FALSE)
  }

  sigma_inf_star <- pmin(pmax(-0.5 * kappa_inf, sigma_L), sigma_U)

  lower <- Aij + sigma_inf_star^2 + kappa_inf * sigma_inf_star

  upper <- pmax(
    Aij + sigma_L^2 + kappa_sup * sigma_L,
    Aij + sigma_U^2 + kappa_sup * sigma_U
  )

  diagonal <- .diagonal_cov_bounds(di, sdi, sigma_L, sigma_U, rho_L[i], rho_U[i])
  lower[is_diag] <- diagonal$lower[is_diag]
  upper[is_diag] <- diagonal$upper[is_diag]

  .new_prism_cov_bounds_pairs(i, j, lower, upper, validated)
}
