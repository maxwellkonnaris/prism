#' Identification bounds for absolute log-covariance
#'
#' For each pair of features `(i, j)`,
#' `Cov(log W_i, log W_j) = Sigma_rel[i,j] + sigma^2 + sigma*(s_i*rho_i + s_j*rho_j)`,
#' where `sigma = SD(log total scale)`, `s_d = sqrt(Sigma_rel[d,d])`, and
#' `rho_d = Cor(log P_d, log total scale)`. `sigma` and `rho` are generally
#' not point-identified from compositional data alone, so this returns the
#' tightest bound on the covariance implied by whatever scale information is
#' supplied:
#'
#' \describe{
#'   \item{no `sigma_L`/`sigma_U`}{unbounded scale: finite lower bound,
#'     infinite upper bound}
#'   \item{`sigma_L`, `sigma_U`}{bounded scale: sharp closed-form bounds}
#'   \item{`sigma_L`, `sigma_U`, and a fixed `rho_L == rho_U`}{sharp
#'     closed-form bounds, after checking the fixed correlation is
#'     compatible with `Sigma_rel`}
#'   \item{`sigma_L`, `sigma_U`, and a genuine `rho_L < rho_U` interval}{
#'     sharp numerical bounds obtained by optimizing over the intersection of
#'     the correlation box and the positive-semidefinite feasibility ellipsoid}
#' }
#'
#' Supplying `rho_L`/`rho_U` without `sigma_L`/`sigma_U` is an error, since
#' correlation information alone does not bound `sigma`. A genuine rho
#' interval requires the optional CVXR package. If `sigma_U = 0`, scale terms
#' vanish and the result is exactly `Sigma_rel`, irrespective of rho.
#'
#' @param Sigma_rel Features-by-features relative log-covariance matrix.
#' @param sigma_L,sigma_U Optional scale-SD bounds, `0 <= sigma_L <= sigma_U`.
#' @param rho_L,rho_U Optional scalar or feature-length correlation bounds
#'   in `[-1, 1]`.
#' @param rho_witness Optional feature-length correlation vector known to
#'   be feasible; supplying it can certify a genuine rho interval as non-empty
#'   without a separate feasibility solve. Exact support optimization still
#'   requires CVXR.
#' @param solver Optional CVXR solver name. `NULL` lets CVXR select its default.
#'   Used for exact support optimization under a genuine rho interval and for
#'   rank-deficient feasibility checks when no witness resolves feasibility.
#' @param eig_tol Relative eigenvalue tolerance used to detect rank
#'   deficiency in `Sigma_rel`.
#'
#' @return A `prism_cov_bounds` object: a list with features-by-features
#'   `lower`/`upper` matrices, `regime`, `sharp`, and the resolved `inputs`.
#' @export
cov_bounds <- function(Sigma_rel, sigma_L = NULL, sigma_U = NULL,
                        rho_L = NULL, rho_U = NULL, rho_witness = NULL,
                        solver = NULL, eig_tol = 1e-10) {
  A <- .validate_relative_covariance(Sigma_rel)
  D <- nrow(A)
  cfg <- .resolve_bound_regime(D, sigma_L, sigma_U, rho_L, rho_U)
  s <- sqrt(diag(A))

  sharp_off_diagonal <- TRUE
  kappa <- sqrt(pmax(outer(diag(A), diag(A), "+") + 2 * A, 0))

  if (is.finite(cfg$sigma_U) && cfg$sigma_U == 0) {
    lower <- upper <- A
    dimnames(lower) <- dimnames(upper) <- dimnames(A)
    return(.new_prism_cov_bounds(
      lower, upper,
      regime = cfg$regime, sharp_off_diagonal = TRUE,
      inputs = list(sigma_L = cfg$sigma_L, sigma_U = cfg$sigma_U, rho_L = cfg$rho_L, rho_U = cfg$rho_U)
    ))
  }

  if (cfg$regime == "bounded_scale_and_correlation") {
    fixed <- all(cfg$rho_L == cfg$rho_U)
    if (fixed) {
      .check_fixed_rho_psd(A, cfg$rho_L)
      t_point <- outer(s * cfg$rho_L, s * cfg$rho_L, "+")
      t_min <- t_max <- t_point
    } else {
      if (!is.null(rho_witness)) rho_witness <- .validate_rho_vec(rho_witness, D, "rho_witness")
      if (!.rho_box_feasible(A, cfg$rho_L, cfg$rho_U, witness = rho_witness, eig_tol = eig_tol, solver = solver)) {
        stop(
          "Empty E intersection C: rho_L/rho_U are not jointly compatible ",
          "with Sigma_rel (no single correlation vector in the box ",
          "satisfies the positive-semidefiniteness constraint).",
          call. = FALSE
        )
      }
      support <- .rho_box_support(A, cfg$rho_L, cfg$rho_U, eig_tol = eig_tol, solver = solver)
      t_min <- support$m
      t_max <- support$M
    }
  } else {
    t_min <- -kappa
    t_max <- kappa
  }

  sigma_star <- pmin(pmax(-0.5 * t_min, cfg$sigma_L), cfg$sigma_U)
  lower <- A + sigma_star^2 + t_min * sigma_star
  if (is.infinite(cfg$sigma_U)) {
    upper <- matrix(Inf, D, D)
  } else {
    upper <- pmax(A + cfg$sigma_L^2 + t_max * cfg$sigma_L, A + cfg$sigma_U^2 + t_max * cfg$sigma_U)
  }
  diag(lower) <- pmax(diag(lower), 0)
  lower <- (lower + t(lower)) / 2
  upper <- (upper + t(upper)) / 2
  dimnames(lower) <- dimnames(upper) <- dimnames(A)

  .new_prism_cov_bounds(
    lower, upper,
    regime = cfg$regime, sharp_off_diagonal = sharp_off_diagonal,
    inputs = list(sigma_L = cfg$sigma_L, sigma_U = cfg$sigma_U, rho_L = cfg$rho_L, rho_U = cfg$rho_U)
  )
}

.new_prism_cov_bounds <- function(lower, upper, regime, sharp_off_diagonal, inputs) {
  structure(
    list(
      lower = lower, upper = upper, regime = regime,
      sharp = list(diagonal = TRUE, off_diagonal = sharp_off_diagonal),
      inputs = inputs
    ),
    class = "prism_cov_bounds"
  )
}

#' @export
print.prism_cov_bounds <- function(x, ...) {
  cat("PRISM covariance bounds\n")
  cat("  regime: ", x$regime, "\n", sep = "")
  cat("  dimensions: ", nrow(x$lower), " x ", nrow(x$lower), "\n", sep = "")
  cat(
    "  off-diagonal bounds: ",
    if (x$sharp$off_diagonal) "sharp" else "conservative (not sharp)", "\n",
    sep = ""
  )
  invisible(x)
}

#' @keywords internal
.validate_relative_covariance <- function(Sigma_rel) {
  A <- as.matrix(Sigma_rel)
  if (!is.numeric(A) || nrow(A) != ncol(A) || nrow(A) < 1L || any(!is.finite(A))) {
    stop("Sigma_rel must be a finite, square numeric matrix.", call. = FALSE)
  }
  scale <- max(1, max(abs(A)))
  if (max(abs(A - t(A))) > 100 * .Machine$double.eps * scale) {
    stop("Sigma_rel must be symmetric.", call. = FALSE)
  }
  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  psd_tol <- 1e-8 * max(1, max(abs(eig)))
  if (min(eig) < -psd_tol) {
    stop("Sigma_rel must be positive semidefinite.", call. = FALSE)
  }
  variance_scale <- max(1, max(abs(diag(A))))
  if (any(!is.finite(diag(A))) || any(diag(A) <= .Machine$double.eps * variance_scale)) {
    stop(
      "Sigma_rel must have strictly positive marginal variances; ",
      "zero or numerically zero marginal variance is not supported.",
      call. = FALSE
    )
  }
  diag(A) <- pmax(diag(A), 0)
  A
}

#' @keywords internal
.resolve_bound_regime <- function(D, sigma_L, sigma_U, rho_L, rho_U) {
  has_sigma <- !is.null(sigma_L) || !is.null(sigma_U)
  has_rho <- !is.null(rho_L) || !is.null(rho_U)
  if (has_rho && !has_sigma) {
    stop("Correlation bounds require scale-SD bounds (sigma_L, sigma_U).", call. = FALSE)
  }
  if (!has_sigma) {
    return(list(regime = "unbounded_scale", sigma_L = 0, sigma_U = Inf, rho_L = NULL, rho_U = NULL))
  }
  if (is.null(sigma_L) || is.null(sigma_U)) {
    stop("sigma_L and sigma_U must be supplied together.", call. = FALSE)
  }
  sigma_L <- .assert_scalar_finite(sigma_L, "sigma_L", lower = 0)
  sigma_U <- .assert_scalar_finite(sigma_U, "sigma_U", lower = 0)
  if (sigma_L > sigma_U) stop("sigma_L must be <= sigma_U.", call. = FALSE)
  if (!has_rho) {
    return(list(regime = "bounded_scale", sigma_L = sigma_L, sigma_U = sigma_U, rho_L = NULL, rho_U = NULL))
  }
  if (is.null(rho_L) || is.null(rho_U)) {
    stop("rho_L and rho_U must be supplied together.", call. = FALSE)
  }
  rho_L <- .validate_rho_vec(rho_L, D, "rho_L")
  rho_U <- .validate_rho_vec(rho_U, D, "rho_U")
  if (any(rho_L > rho_U)) stop("rho_L must be <= rho_U elementwise.", call. = FALSE)
  list(regime = "bounded_scale_and_correlation", sigma_L = sigma_L, sigma_U = sigma_U, rho_L = rho_L, rho_U = rho_U)
}

.validate_rho_vec <- function(x, D, name) {
  ok <- is.numeric(x) && length(x) %in% c(1L, D) && all(is.finite(x)) && all(x >= -1) && all(x <= 1)
  if (!ok) stop(name, " must be finite, length 1 or D, and within [-1, 1].", call. = FALSE)
  if (length(x) == 1L) rep(as.numeric(x), D) else as.numeric(x)
}

#' PSD compatibility of a fixed correlation vector
#'
#' Checks that the augmented matrix `[[Sigma_rel, u], [u', 1]]`, with
#' `u = sqrt(diag(Sigma_rel)) * rho`, stays positive semidefinite -- the
#' exact condition for `rho` to be a realizable joint correlation between
#' log total scale and each feature's log composition.
#' @keywords internal
.check_fixed_rho_psd <- function(Sigma_rel, rho) {
  u <- sqrt(diag(Sigma_rel)) * rho
  aug <- rbind(cbind(Sigma_rel, u), c(u, 1))
  eig <- eigen(aug, symmetric = TRUE, only.values = TRUE)$values
  tol <- 1e-8 * max(1, max(abs(eig)))
  if (min(eig) < -tol) {
    stop(
      "Fixed correlations are incompatible with Sigma_rel: the augmented ",
      "covariance matrix is not positive semidefinite.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.require_cvxr_solver <- function(solver, purpose) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop(purpose, " requires the optional CVXR package.", call. = FALSE)
  }
  if (is.null(solver)) return(invisible(TRUE))
  if (!is.character(solver) || length(solver) != 1L || is.na(solver) || !nzchar(solver)) {
    stop("solver must be NULL or a single non-empty string.", call. = FALSE)
  }
  available <- CVXR::installed_solvers()
  if (!(solver %in% available)) {
    stop(
      "Solver '", solver, "' is not installed for CVXR. Installed solvers: ",
      paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.cvxr_psolve <- function(problem, solver = NULL) {
  result <- if (is.null(solver)) {
    CVXR::psolve(problem)
  } else {
    CVXR::psolve(problem, solver = solver)
  }
  if (is.list(result) && !is.null(result$status)) {
    return(list(status = result$status, value = result$value))
  }
  exports <- getNamespaceExports("CVXR")
  accessor <- intersect(c("status", "problem_status"), exports)
  if (!length(accessor)) {
    stop("CVXR did not expose a problem-status accessor.", call. = FALSE)
  }
  status <- getExportedValue("CVXR", accessor[1L])(problem)
  list(status = status, value = result)
}

#' Exact support of the SPSD-feasible rho set intersected with a box
#'
#' Computes, for every pair `(i, j)`, the exact extrema of
#' `s_i * rho_i + s_j * rho_j` over all `rho_L <= rho <= rho_U` satisfying
#' the singular-safe SPSD conditions `u in Range(Sigma_rel)` and
#' `u' Sigma_rel^+ u <= 1`, where `u = s * rho`.
#'
#' Two exact shortcuts avoid most conic solves:
#' \enumerate{
#'   \item Closed-form certificate. Over the ellipsoid alone, the maximum of
#'     `c'u` (with `c = e_i + e_j`) is `kappa = sqrt(c' Sigma_rel c)`,
#'     attained at `u* = Sigma_rel c / kappa`. Because the box-constrained
#'     feasible set is a subset of the ellipsoid, if `u*` also satisfies the
#'     box then it is optimal there too and `M_ij = kappa` exactly. The
#'     minimum is handled the same way at `-u*`. Every coordinate of `u*` is
#'     checked, not only `i` and `j`.
#'   \item Central symmetry. When the box is symmetric about zero
#'     (`rho_L == -rho_U`), the feasible set is centrally symmetric, so
#'     `m_ij = -M_ij` and only one solve is needed per pair.
#' }
#' CVXR is only required if at least one pair fails the certificate.
#'
#' @param shortcuts If `FALSE`, skip both shortcuts and solve every direction
#'   with CVXR (reference implementation, used for testing).
#' @return `list(m, M, n_solves)`, where `n_solves` counts CVXR calls.
#' @keywords internal
.rho_box_support <- function(Sigma_rel, rho_L, rho_U,
                             eig_tol = 1e-10, solver = NULL,
                             shortcuts = TRUE) {
  D <- nrow(Sigma_rel)
  s <- sqrt(diag(Sigma_rel))
  u_L <- s * rho_L
  u_U <- s * rho_U

  eig <- eigen(Sigma_rel, symmetric = TRUE)
  eig_scale <- max(1, max(abs(eig$values)))
  keep <- eig$values > eig_tol * eig_scale
  if (!any(keep)) {
    stop("Sigma_rel has no positive eigendirections at the requested eig_tol.", call. = FALSE)
  }
  V_r <- eig$vectors[, keep, drop = FALSE]
  lambda_r <- eig$values[keep]
  factor <- V_r * rep(sqrt(lambda_r), each = D)
  # Range-truncated covariance: the ellipsoid's shape matrix, consistent with
  # the CVXR parameterisation u = factor %*% y, ||y|| <= 1.
  A_r <- tcrossprod(factor)

  box_tol <- 1e-12 * max(1, max(s))
  symmetric_box <- max(abs(u_L + u_U)) <= box_tol

  cvxr_ready <- FALSE
  n_solves <- 0L
  y <- NULL
  u <- NULL
  constraints <- NULL
  solve_direction <- function(direction) {
    if (!cvxr_ready) {
      .require_cvxr_solver(solver, "Exact bounds for a non-fixed rho interval")
      y <<- CVXR::Variable(length(lambda_r))
      u <<- factor %*% y
      constraints <<- list(CVXR::p_norm(y, 2) <= 1, u >= u_L, u <= u_U)
      cvxr_ready <<- TRUE
    }
    n_solves <<- n_solves + 1L
    objective <- CVXR::sum_entries(as.numeric(direction) * u)
    problem <- CVXR::Problem(CVXR::Maximize(objective), constraints)
    solved <- .cvxr_psolve(problem, solver)
    status <- tolower(as.character(solved$status)[1L])
    value <- as.numeric(solved$value)[1L]
    if (!(status %in% c("optimal", "optimal_inaccurate")) || !is.finite(value)) {
      stop("Exact rho-support optimization failed with CVXR status '", status, "'.", call. = FALSE)
    }
    value
  }

  m <- M <- matrix(NA_real_, D, D, dimnames = dimnames(Sigma_rel))
  a_diag <- diag(A_r)
  for (i in seq_len(D)) {
    js <- i:D
    # Column k of `num` is Sigma_r (e_i + e_j) for j = js[k].
    num <- A_r[, js, drop = FALSE] + A_r[, i]
    kappa <- sqrt(pmax(a_diag[i] + a_diag[js] + 2 * A_r[i, js], 0))
    if (shortcuts) {
      degenerate <- kappa <= sqrt(.Machine$double.eps) * sqrt(eig_scale)
      u_star <- num / rep(ifelse(degenerate, 1, kappa), each = D)
      max_ok <- colSums(u_star < u_L - box_tol | u_star > u_U + box_tol) == 0
      min_ok <- colSums(-u_star < u_L - box_tol | -u_star > u_U + box_tol) == 0
    } else {
      degenerate <- max_ok <- min_ok <- rep(FALSE, length(js))
    }
    for (k in seq_along(js)) {
      j <- js[k]
      direction <- numeric(D)
      direction[i] <- direction[i] + 1
      direction[j] <- direction[j] + 1
      if (degenerate[k]) {
        # c is orthogonal to Range(Sigma_rel): c'u = 0 on the whole feasible set.
        M_ij <- m_ij <- 0
      } else {
        M_ij <- if (max_ok[k]) kappa[k] else solve_direction(direction)
        m_ij <- if (min_ok[k]) {
          -kappa[k]
        } else if (shortcuts && symmetric_box) {
          -M_ij
        } else {
          -solve_direction(-direction)
        }
      }
      m[i, j] <- m[j, i] <- m_ij
      M[i, j] <- M[j, i] <- M_ij
    }
  }
  scale <- max(1, max(abs(c(m, M))))
  if (any(m > M + 1e-6 * scale)) {
    stop("CVXR returned inconsistent lower and upper rho-support values.", call. = FALSE)
  }
  crossed <- m > M
  if (any(crossed)) {
    midpoint <- (m[crossed] + M[crossed]) / 2
    m[crossed] <- midpoint
    M[crossed] <- midpoint
  }
  list(m = m, M = M, n_solves = n_solves)
}

#' Is a rho box compatible with Sigma_rel?
#'
#' Certifies that the box `C` (all rho with `rho_L <= rho <= rho_U`)
#' contains at least one point whose `u = sqrt(diag(Sigma_rel)) * rho`
#' keeps the augmented covariance matrix positive semidefinite, i.e. that
#' `u` lies in the ellipsoid `E` (`u` in `Range(Sigma_rel)` with
#' `u' Sigma_rel^+ u <= 1`).
#'
#' @keywords internal
.rho_box_feasible <- function(Sigma_rel, rho_L, rho_U, witness = NULL,
                               eig_tol = 1e-10, solver = NULL) {
  s <- sqrt(diag(Sigma_rel))
  u_L <- s * rho_L
  u_U <- s * rho_U
  if (all(u_L <= 1e-12) && all(u_U >= -1e-12)) return(TRUE)

  eig <- eigen(Sigma_rel, symmetric = TRUE)
  eig_scale <- max(1, max(abs(eig$values)))
  keep <- eig$values > eig_tol * eig_scale
  rank_A <- sum(keep)
  V_r <- eig$vectors[, keep, drop = FALSE]
  lambda_r <- eig$values[keep]

  in_ellipsoid <- function(u) {
    if (any(u < u_L - 1e-8) || any(u > u_U + 1e-8)) return(FALSE)
    if (rank_A == 0L) return(all(abs(u) < 1e-8))
    coord <- as.numeric(crossprod(V_r, u))
    resid <- u - V_r %*% coord
    if (max(abs(resid)) > 1e-6 * max(1, max(abs(u)))) return(FALSE)
    sum(coord^2 / lambda_r) <= 1 + 1e-8
  }

  if (!is.null(witness) && in_ellipsoid(s * witness)) return(TRUE)
  if (rank_A == 0L) return(FALSE)

  if (rank_A == length(u_L)) {
    Sigma_inv <- solve(Sigma_rel)
    u0 <- pmin(pmax(0, u_L), u_U)
    fit <- stats::optim(
      u0,
      fn = function(u) as.numeric(t(u) %*% Sigma_inv %*% u),
      gr = function(u) as.numeric(2 * Sigma_inv %*% u),
      method = "L-BFGS-B", lower = u_L, upper = u_U,
      control = list(factr = 1e7, pgtol = 1e-14, maxit = 2000)
    )
    if (fit$convergence != 0) {
      stop("Rho-box feasibility check did not converge.", call. = FALSE)
    }
    return(fit$value <= 1 + 1e-6)
  }

  .require_cvxr_solver(
    solver,
    "Certifying a rank-deficient rho box without a feasible witness"
  )
  M <- V_r * rep(sqrt(lambda_r), each = nrow(V_r))
  y <- CVXR::Variable(rank_A)
  u <- M %*% y
  problem <- CVXR::Problem(
    CVXR::Minimize(0),
    list(CVXR::p_norm(y, 2) <= 1, u >= u_L, u <= u_U)
  )
  result <- .cvxr_psolve(problem, solver)
  isTRUE(tolower(as.character(result$status)[1L]) %in% c("optimal", "optimal_inaccurate"))
}
