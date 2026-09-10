#' Identification bounds for absolute log-covariance
#'
#' For each pair of features `(i, j)`,
#' `Cov(log W_i, log W_j) = Sigma_rel[i,j] + sigma^2 + sigma*(s_i*rho_i + s_j*rho_j)`,
#' where `sigma = SD(log total scale)`, `s_d = sqrt(Sigma_rel[d,d])`, and
#' `rho_d = Cor(log P_d, log total scale)`. `sigma` and `rho` are generally
#' not point-identified from compositional data alone, so this returns the
#' tightest closed-form bound on the covariance implied by whatever scale
#' information is supplied:
#'
#' \describe{
#'   \item{no `sigma_L`/`sigma_U`}{unbounded scale: finite lower bound,
#'     infinite upper bound}
#'   \item{`sigma_L`, `sigma_U`}{bounded scale: sharp closed-form bounds}
#'   \item{`sigma_L`, `sigma_U`, and a fixed `rho_L == rho_U`}{sharp
#'     closed-form bounds, after checking the fixed correlation is
#'     compatible with `Sigma_rel`}
#'   \item{`sigma_L`, `sigma_U`, and a genuine `rho_L < rho_U` interval}{
#'     conservative (not sharp) closed-form bounds, after checking the
#'     correlation box is compatible with `Sigma_rel`}
#' }
#'
#' Supplying `rho_L`/`rho_U` without `sigma_L`/`sigma_U` is an error, since
#' correlation information alone does not bound `sigma`.
#'
#' @param Sigma_rel Features-by-features relative log-covariance matrix.
#' @param sigma_L,sigma_U Optional scale-SD bounds, `0 <= sigma_L <= sigma_U`.
#' @param rho_L,rho_U Optional scalar or feature-length correlation bounds
#'   in `[-1, 1]`.
#' @param rho_witness Optional feature-length correlation vector known to
#'   be feasible; supplying it lets a genuine rho interval be certified
#'   non-empty without an optimizer.
#' @param solver CVXR solver used only when `Sigma_rel` is rank-deficient
#'   and no witness is supplied or works; see Details.
#' @param eig_tol Relative eigenvalue tolerance used to detect rank
#'   deficiency in `Sigma_rel`.
#'
#' @return A `prism_cov_bounds` object: a list with features-by-features
#'   `lower`/`upper` matrices, `regime`, `sharp` (whether the diagonal and
#'   off-diagonal bounds are known to be the tightest possible, versus
#'   conservative), and the resolved `inputs`.
#' @export
cov_bounds <- function(Sigma_rel, sigma_L = NULL, sigma_U = NULL,
                        rho_L = NULL, rho_U = NULL, rho_witness = NULL,
                        solver = "ECOS", eig_tol = 1e-10) {
  A <- .validate_relative_covariance(Sigma_rel)
  D <- nrow(A)
  cfg <- .resolve_bound_regime(D, sigma_L, sigma_U, rho_L, rho_U)
  s <- sqrt(diag(A))

  sharp_off_diagonal <- TRUE
  kappa <- sqrt(pmax(outer(diag(A), diag(A), "+") + 2 * A, 0))

  if (cfg$regime == "bounded_scale_and_correlation") {
    zero_var <- s == 0
    if (any(zero_var & (cfg$rho_L != 0 | cfg$rho_U != 0))) {
      stop("rho bounds for a zero-variance feature must be exactly 0.", call. = FALSE)
    }
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
      sharp_off_diagonal <- FALSE
      box_pos <- outer(s * cfg$rho_U, s * cfg$rho_U, "+")
      box_neg <- -outer(s * cfg$rho_L, s * cfg$rho_L, "+")
      t_max <- pmin(kappa, box_pos)
      t_min <- -pmin(kappa, box_neg)
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
                               eig_tol = 1e-10, solver = "ECOS") {
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

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop(
      "Sigma_rel is rank-deficient; certifying this rho box requires the ",
      "CVXR package, or a feasible rho_witness supplied directly.",
      call. = FALSE
    )
  }
  if (!(solver %in% CVXR::installed_solvers())) {
    stop("Solver '", solver, "' is not installed for CVXR.", call. = FALSE)
  }
  M <- V_r * rep(sqrt(lambda_r), each = nrow(V_r))
  y <- CVXR::Variable(rank_A)
  u <- M %*% y
  problem <- CVXR::Problem(
    CVXR::Minimize(0),
    list(CVXR::p_norm(y, 2) <= 1, u >= u_L, u <= u_U)
  )
  result <- CVXR::psolve(problem, solver = solver)
  isTRUE(result$status %in% c("optimal", "optimal_inaccurate"))
}
