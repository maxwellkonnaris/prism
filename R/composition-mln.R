.default_mln_options <- function() {
  list(
    design = NULL,
    tune_draws = 200L,
    center_scale_grid = c(1, 1.5, 2, 2.5, 3),
    m_grid = 10L,
    prior_mode = "tuned",
    fixed_Gamma = NULL,
    fixed_Xi = NULL,
    fixed_upsilon = NULL,
    prior_label = NULL,
    omega_structure_grid = c("ALR_I", "G_I_Gt"),
    gamma_grid = 10^seq(-2, 2, by = 0.5),
    max_hessian_gb = Inf,
    return_grid_pairwise = FALSE,
    return_grid_fits = FALSE,
    return_fit = FALSE
  )
}

.resolve_mln_options <- function(user, defaults) {
  unknown <- setdiff(names(user), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown MLN option(s): ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  utils::modifyList(defaults, user)
}

.fit_composition_mln <- function(Y, n_draws, pseudocount, ci_alpha, options, verbose) {
  if (!requireNamespace("fido", quietly = TRUE)) {
    stop("The MLN composition estimator requires the 'fido' package.", call. = FALSE)
  }
  D <- nrow(Y)
  N <- ncol(Y)
  props <- sweep(Y, 2, colSums(Y), "/")
  median_dominance <- stats::median(apply(props, 2, max))
  if (median_dominance > 0.99) {
    stop(
      sprintf(
        "MLN fitting is ill-conditioned: median maximum sample proportion is %.3f (> 0.99).",
        median_dominance
      ),
      call. = FALSE
    )
  }
  hessian_gb <- ((D - 1L) * N)^2 * 8 / 1e9
  if (hessian_gb > options$max_hessian_gb) {
    stop(
      sprintf(
        "Estimated MLN Hessian storage is %.2f GB, above max_hessian_gb=%.2f.",
        hessian_gb,
        options$max_hessian_gb
      ),
      call. = FALSE
    )
  }

  otu_table <- Y + pseudocount
  design <- .build_mln_design(options$design, N)
  Q <- nrow(design)
  theta <- matrix(0, nrow = D - 1L, ncol = Q)
  tuning <- NULL

  if (options$prior_mode == "fixed") {
    required <- c("fixed_Gamma", "fixed_Xi", "fixed_upsilon")
    if (any(vapply(options[required], is.null, logical(1)))) {
      stop("Fixed MLN priors require fixed_Gamma, fixed_Xi, and fixed_upsilon.", call. = FALSE)
    }
    Gamma <- options$fixed_Gamma
    Xi <- options$fixed_Xi
    upsilon <- as.integer(options$fixed_upsilon)
  } else {
    tuning <- .tune_mln_gamma_centered_prior(
      otu_table = otu_table,
      X = design,
      D = D,
      Theta = theta,
      S = options$tune_draws,
      gamma_grid = options$gamma_grid,
      center_scale_grid = options$center_scale_grid,
      m_grid = options$m_grid,
      omega_structure_grid = options$omega_structure_grid,
      ci_alpha = ci_alpha,
      return_grid_pairwise = options$return_grid_pairwise,
      return_grid_fits = options$return_grid_fits,
      verbose = verbose
    )
    Gamma <- tuning$best_Gamma
    Xi <- tuning$best_Xi
    upsilon <- tuning$best_upsilon
  }

  priors <- fido::pibble(
    otu_table,
    design,
    upsilon,
    theta,
    Gamma,
    Xi,
    n_samples = n_draws
  )
  fit <- fido::refit(priors, optim_method = "lbfgs")
  draws <- fido::to_proportions(fit)$Eta
  if (!is.numeric(draws) || length(dim(draws)) != 3L || dim(draws)[3] < 2L ||
      any(!is.finite(draws)) || any(draws <= 0)) {
    stop("MLN fitting did not return valid positive D x N x S composition draws.", call. = FALSE)
  }

  diagnostics <- list(
    design_matrix = design,
    estimated_hessian_gb = hessian_gb,
    median_maximum_proportion = median_dominance,
    prior_mode = options$prior_mode,
    prior_label = options$prior_label,
    selected_Gamma = Gamma,
    selected_Xi = Xi,
    upsilon = upsilon,
    grid_results = if (is.null(tuning)) NULL else tuning$grid_results,
    grid_fits = if (is.null(tuning)) NULL else tuning$grid_fits
  )
  list(draws = draws, fit = fit, fit_summary = summary(fit), diagnostics = diagnostics)
}

.build_mln_design <- function(X, N) {
    # If no design matrix is supplied, use an intercept-only model.
    if (is.null(X)) {
        X_out <- matrix(1, nrow = 1L, ncol = N)
        rownames(X_out) <- "(Intercept)"
        return(X_out)
    }

    X_out <- as.matrix(X)

    if (!is.numeric(X_out) || any(!is.finite(X_out))) {
        stop("mln_design must be a finite numeric matrix.")
    }

    # Accept either Q x N or N x Q.
    if (ncol(X_out) == N) {
        # already Q x N
    } else if (nrow(X_out) == N) {
        X_out <- t(X_out)
    } else {
        stop("mln_design must have either ncol = N or nrow = N.")
    }

    if (nrow(X_out) < 1L) {
        stop("mln_design must contain at least one row.")
    }

    X_out
}

.coerce_mln_gamma_diag <- function(gamma_diag, Q) {
    if (is.null(gamma_diag)) {
        return(NULL)
    }

    gamma_diag <- as.numeric(gamma_diag)

    if (length(gamma_diag) == 1L) {
        gamma_diag <- rep(gamma_diag, Q)
    }

    if (length(gamma_diag) != Q) {
        stop("mln_gamma_diag must have length 1 or Q = ", Q, ".")
    }

    if (any(!is.finite(gamma_diag)) || any(gamma_diag <= 0)) {
        stop("mln_gamma_diag must contain positive finite values.")
    }

    gamma_diag
}

.make_mln_G <- function(D, ref = D) {
    nonref <- setdiff(seq_len(D), ref)
    G <- matrix(0, nrow = D - 1L, ncol = D)

    for (k in seq_along(nonref)) {
        G[k, nonref[k]] <- 1
        G[k, ref] <- -1
    }

    G
}

.make_mln_A <- function(D, structure) {
    G <- .make_mln_G(D)
    K <- D - 1L

    A <- switch(structure,
        "ALR_I" = diag(K),
        "G_I_Gt" = {
            G %*% diag(D) %*% t(G)
        },
        "G_diagTau_Gt" = {
            tau2 <- exp(seq(log(0.5^2), log(2^2), length.out = D))
            tau2 <- tau2 / mean(tau2)
            G %*% diag(tau2, D) %*% t(G)
        },
        stop("Unknown omega_structure: ", structure)
    )

    A / mean(diag(A))
}

.summarise_mln_pairwise_rel_from_proportions <- function(Pi_draws, ci_alpha = 0.05) {
    if (is.null(Pi_draws) || !is.numeric(Pi_draws) || length(dim(Pi_draws)) != 3L) {
        return(NULL)
    }
    D <- dim(Pi_draws)[1]
    S_draws <- dim(Pi_draws)[3]
    if (D < 2L || S_draws < 1L) return(NULL)

    idx <- which(upper.tri(matrix(0, D, D), diag = FALSE), arr.ind = TRUE)
    rel_draws <- matrix(NA_real_, nrow = S_draws, ncol = nrow(idx))
    for (s in seq_len(S_draws)) {
        Pi_s <- Pi_draws[, , s]
        if (any(!is.finite(Pi_s)) || any(Pi_s <= 0)) next
        Sigma_rel_s <- stats::cov(t(log(Pi_s)))
        rel_draws[s, ] <- Sigma_rel_s[idx]
    }

    q_lo <- ci_alpha / 2
    q_hi <- 1 - ci_alpha / 2
    safe_stat <- function(x, fun) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) return(NA_real_)
        fun(x)
    }
    rel_hat <- apply(rel_draws, 2, safe_stat, fun = stats::median)
    rel_ci_lower <- apply(rel_draws, 2, safe_stat, fun = function(x) {
        stats::quantile(x, probs = q_lo, names = FALSE, type = 7)
    })
    rel_ci_upper <- apply(rel_draws, 2, safe_stat, fun = function(x) {
        stats::quantile(x, probs = q_hi, names = FALSE, type = 7)
    })

    data.frame(
        i = idx[, 1],
        j = idx[, 2],
        rel_hat = rel_hat,
        rel_ci_lower = rel_ci_lower,
        rel_ci_upper = rel_ci_upper,
        rel_ci_width = rel_ci_upper - rel_ci_lower
    )
}

.tune_mln_gamma_centered_prior <- function(
  otu_table,
  X,
  D,
  Theta,
  S,
  gamma_grid,
  center_scale_grid,
  m_grid,
  omega_structure_grid,
  ci_alpha = 0.05,
  return_grid_pairwise = FALSE,
  return_grid_fits = FALSE,
  verbose = TRUE
) {
    gamma_grid <- sort(unique(as.numeric(gamma_grid)))
    center_scale_grid <- sort(unique(as.numeric(center_scale_grid)))
    m_grid <- sort(unique(as.integer(m_grid)))
    omega_structure_grid <- unique(as.character(omega_structure_grid))

    if (length(gamma_grid) < 1L || any(!is.finite(gamma_grid)) || any(gamma_grid <= 0)) {
        stop("mln_gamma_grid must contain positive finite values.")
    }

    if (length(center_scale_grid) < 1L || any(!is.finite(center_scale_grid)) || any(center_scale_grid <= 0)) {
        stop("mln_center_scale_grid must contain positive finite values.")
    }

    if (length(m_grid) < 1L || any(!is.finite(m_grid)) || any(m_grid <= 0)) {
        stop("mln_m_grid must contain positive integer values.")
    }

    if (length(omega_structure_grid) < 1L) {
        stop("mln_omega_structure_grid must contain at least one structure.")
    }

    Q <- nrow(X)
    K <- nrow(Theta)

    if (K != D - 1L) {
        stop("Theta has ", K, " rows, but expected D - 1 = ", D - 1L, ".")
    }

    grid_results <- expand.grid(
        gamma_scale = gamma_grid,
        center_scale = center_scale_grid,
        m = m_grid,
        omega_structure = omega_structure_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    grid_results$upsilon <- D + as.integer(grid_results$m)
    grid_results$log_marginal_likelihood <- NA_real_
    grid_results$success <- FALSE
    grid_results$error <- NA_character_
    if (isTRUE(return_grid_pairwise)) {
        grid_results$pairwise_rel <- rep(list(NULL), nrow(grid_results))
        grid_results$refit_success <- NA
        grid_results$refit_error <- NA_character_
    }

    fit_list <- vector("list", nrow(grid_results))
    A_list <- vector("list", nrow(grid_results))
    Xi_list <- vector("list", nrow(grid_results))
    Gamma_list <- vector("list", nrow(grid_results))
    grid_fit_payloads <- if (isTRUE(return_grid_fits)) {
        vector("list", nrow(grid_results))
    } else {
        NULL
    }

    for (ii in seq_len(nrow(grid_results))) {
        g <- grid_results$gamma_scale[ii]
        center_scale <- grid_results$center_scale[ii]
        m <- as.integer(grid_results$m[ii])
        omega_structure <- grid_results$omega_structure[ii]
        upsilon_i <- D + m

        A_i <- .make_mln_A(D, omega_structure)

        if (!all(dim(A_i) == c(K, K))) {
            stop("A matrix has wrong dimensions for structure ", omega_structure, ".")
        }

        Gamma_i <- g * diag(Q)
        Xi_i <- m * center_scale * A_i

        if (isTRUE(verbose)) {
            message(sprintf(
                paste0(
                    "[prism MLN] centered prior tuning | candidate %d/%d | ",
                    "gamma = %.6g | tau^2 = %.6g | m = %d | upsilon = %d | A = %s"
                ),
                ii, nrow(grid_results), g, center_scale, m, upsilon_i, omega_structure
            ))
        }

        fit_i <- tryCatch(
            fido::pibble(
                otu_table,
                X,
                upsilon_i,
                Theta,
                Gamma_i,
                Xi_i,
                n_samples = S
            ),
            error = function(e) e
        )

        if (inherits(fit_i, "error")) {
            grid_results$error[ii] <- conditionMessage(fit_i)
            next
        }

        lml_i <- fit_i$logMarginalLikelihood

        if (is.null(lml_i) || !is.finite(lml_i)) {
            grid_results$error[ii] <- "logMarginalLikelihood missing or non-finite"
            next
        }

        fit_list[[ii]] <- fit_i
        A_list[[ii]] <- A_i
        Xi_list[[ii]] <- Xi_i
        Gamma_list[[ii]] <- Gamma_i

        grid_results$log_marginal_likelihood[ii] <- as.numeric(lml_i)
        grid_results$success[ii] <- TRUE

        if (isTRUE(return_grid_pairwise) || isTRUE(return_grid_fits)) {
            refit_i <- tryCatch(
                fido::refit(fit_i, optim_method = "lbfgs"),
                error = function(e) e
            )
            pairwise_i <- NULL
            refit_error <- NULL
            if (inherits(refit_i, "error")) {
                refit_error <- conditionMessage(refit_i)
            } else {
                Pi_i <- tryCatch(
                    fido::to_proportions(refit_i)$Eta,
                    error = function(e) e
                )
                if (inherits(Pi_i, "error")) {
                    refit_error <- conditionMessage(Pi_i)
                } else {
                    pairwise_i <- .summarise_mln_pairwise_rel_from_proportions(
                        Pi_i,
                        ci_alpha = ci_alpha
                    )
                }
            }
            if (isTRUE(return_grid_pairwise)) {
                grid_results$pairwise_rel[[ii]] <- pairwise_i
                grid_results$refit_success[ii] <- is.null(refit_error)
                grid_results$refit_error[ii] <- if (is.null(refit_error)) {
                    NA_character_
                } else {
                    refit_error
                }
            }
            if (isTRUE(return_grid_fits)) {
                grid_fit_payloads[[ii]] <- list(
                    prior_fit_object = fit_i,
                    mln_fit_object = if (inherits(refit_i, "error")) NULL else refit_i,
                    pairwise_rel = pairwise_i,
                    error = refit_error
                )
            }
        }
    }

    ok <- which(is.finite(grid_results$log_marginal_likelihood))

    if (length(ok) == 0L) {
        stop(
            "All centered-prior MLN candidates failed.\n",
            paste(
                sprintf(
                    "(gamma=%s, tau^2=%s, m=%s, A=%s) -> %s",
                    grid_results$gamma_scale,
                    grid_results$center_scale,
                    grid_results$m,
                    grid_results$omega_structure,
                    grid_results$error
                ),
                collapse = "\n"
            )
        )
    }

    best_idx <- ok[which.max(grid_results$log_marginal_likelihood[ok])]

    list(
        best_index = best_idx,
        best_gamma_scale = grid_results$gamma_scale[best_idx],
        best_center_scale = grid_results$center_scale[best_idx],
        best_m = as.integer(grid_results$m[best_idx]),
        best_upsilon = as.integer(grid_results$upsilon[best_idx]),
        best_omega_structure = grid_results$omega_structure[best_idx],
        best_Gamma = Gamma_list[[best_idx]],
        best_Xi = Xi_list[[best_idx]],
        best_A = A_list[[best_idx]],
        best_prior_fit = fit_list[[best_idx]],
        best_log_marginal_likelihood = grid_results$log_marginal_likelihood[best_idx],
        gamma_grid = gamma_grid,
        center_scale_grid = center_scale_grid,
        m_grid = m_grid,
        omega_structure_grid = omega_structure_grid,
        grid_results = grid_results,
        grid_fits = grid_fit_payloads
    )
}

.tune_mln_gamma_xi_identity <- function(
  otu_table,
  X,
  upsilon,
  Theta,
  S,
  gamma_grid,
  xi_grid,
  verbose = TRUE
) {
    gamma_grid <- sort(unique(as.numeric(gamma_grid)))
    xi_grid <- sort(unique(as.numeric(xi_grid)))

    if (length(gamma_grid) < 1L || any(!is.finite(gamma_grid)) || any(gamma_grid <= 0)) {
        stop("mln_gamma_grid must contain positive finite values.")
    }
    if (length(xi_grid) < 1L || any(!is.finite(xi_grid)) || any(xi_grid <= 0)) {
        stop("mln_xi_grid must contain positive finite values.")
    }

    Q <- nrow(X)
    K <- nrow(Theta)

    grid_results <- expand.grid(
        gamma_scale = gamma_grid,
        xi_scale = xi_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    grid_results$log_marginal_likelihood <- NA_real_
    grid_results$success <- FALSE
    grid_results$error <- NA_character_

    fit_list <- vector("list", nrow(grid_results))

    for (ii in seq_len(nrow(grid_results))) {
        g <- grid_results$gamma_scale[ii]
        x <- grid_results$xi_scale[ii]

        Gamma_i <- g * diag(Q)
        Xi_i <- x * diag(K)

        if (isTRUE(verbose)) {
            message(sprintf(
                paste0(
                    "[prism MLN] joint tuning Gamma and Xi | ",
                    "candidate %d/%d | gamma = %.6g | xi = %.6g"
                ),
                ii, nrow(grid_results), g, x
            ))
        }

        fit_i <- tryCatch(
            fido::pibble(
                otu_table,
                X,
                upsilon,
                Theta,
                Gamma_i,
                Xi_i,
                n_samples = S
            ),
            error = function(e) e
        )

        if (inherits(fit_i, "error")) {
            grid_results$error[ii] <- conditionMessage(fit_i)
            next
        }

        lml_i <- fit_i$logMarginalLikelihood
        if (is.null(lml_i) || !is.finite(lml_i)) {
            grid_results$error[ii] <- "logMarginalLikelihood missing or non-finite"
            next
        }

        fit_list[[ii]] <- fit_i
        grid_results$log_marginal_likelihood[ii] <- as.numeric(lml_i)
        grid_results$success[ii] <- TRUE
    }

    ok <- which(is.finite(grid_results$log_marginal_likelihood))
    if (length(ok) == 0L) {
        stop(
            "All (Gamma, Xi) candidates failed in MLN tuning.\n",
            paste(
                sprintf(
                    "(gamma=%s, xi=%s) -> %s",
                    grid_results$gamma_scale,
                    grid_results$xi_scale,
                    grid_results$error
                ),
                collapse = "\n"
            )
        )
    }

    best_idx <- ok[which.max(grid_results$log_marginal_likelihood[ok])]
    best_g <- grid_results$gamma_scale[best_idx]
    best_x <- grid_results$xi_scale[best_idx]

    list(
        best_index = best_idx,
        best_gamma_scale = best_g,
        best_xi_scale = best_x,
        best_Gamma = best_g * diag(Q),
        best_Xi = best_x * diag(K),
        best_prior_fit = fit_list[[best_idx]],
        best_log_marginal_likelihood = grid_results$log_marginal_likelihood[best_idx],
        gamma_grid = gamma_grid,
        xi_grid = xi_grid,
        grid_results = grid_results
    )
}


#' Compute true relative log-covariance from known abundances
#'
#' @description

