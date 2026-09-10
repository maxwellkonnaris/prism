.resolve_prism_bootstrap_mode <- function(bootstrap) {
    modes <- c("none", "both")
    if (is.logical(bootstrap) && length(bootstrap) == 1L && !is.na(bootstrap)) {
        return(if (isTRUE(bootstrap)) "both" else "none")
    }
    if (!is.character(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap)) {
        stop(
            "'bootstrap' must be TRUE, FALSE, 'none', or 'both'.",
            call. = FALSE
        )
    }
    if (bootstrap %in% c("composition", "scale")) {
        stop(
            paste0(
                "bootstrap='", bootstrap, "' is deprecated and unsupported. ",
                "Finite-sampling uncertainty must resample complete paired sample units; ",
                "use bootstrap='both' (or TRUE), or bootstrap='none' (or FALSE)."
            ),
            call. = FALSE
        )
    }
    match.arg(bootstrap, modes)
}

.prism_engine <- function(
  counts,
  counts_input = c("counts", "proportions"),
  composition = prism_composition_dirichlet(),
  bootstrap = TRUE,
  S = 1000L,
  max_attempts_per_draw = 100L,
  prevalence = 0,
  verbose = FALSE,
  scale = NULL,
  scale_log = NULL,
  scale_log_ci_level = 0,
  scale_log_rho = TRUE,
  scale_log_lower_zero = FALSE,
  ci_alpha = 0.05,
  sigma_L = NULL,
  sigma_U = NULL,
  rho_L = NULL,
  rho_U = NULL,
  delta = 0.1,
  p_adjust_method = "BH",
  seed = 1L,
  scale_sd_fallback = 0.1,
  scale_uncertainty_draws = 200L,
  return_rel_draws = FALSE,
  return_bound_draws = FALSE,
  stream = NULL,
  workers = 1L,
  high_resolution = FALSE,
  high_resolution_policy = c("borderline", "all", "never"),
  high_resolution_eps = 0.05,
  high_resolution_pair_index = NULL,
  solver = "ECOS",
  high_resolution_eig_tol = 1e-10
) {
    counts_input <- match.arg(counts_input)
    if (!is.logical(scale_log_rho) || length(scale_log_rho) != 1L || is.na(scale_log_rho)) {
        stop("scale_log_rho must be a single TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(scale_log_lower_zero) || length(scale_log_lower_zero) != 1L || is.na(scale_log_lower_zero)) {
        stop("scale_log_lower_zero must be a single TRUE or FALSE.", call. = FALSE)
    }
    composition_estimator <- .resolve_composition_estimator(
        composition,
        alpha = NULL,
        concentration = NULL
    )
    composition <- switch(
        composition_estimator$type,
        dirichlet_multinomial = "Dirichlet",
        fixed = "none",
        mln = "MLN",
        custom = "custom",
        stop("Unsupported composition estimator type.", call. = FALSE)
    )
    alpha <- composition_estimator$options$pseudocount %||% 1
    dirichlet_concentration <- composition_estimator$options$concentration %||% 1
    high_resolution_config <- .resolve_high_resolution(
        high_resolution,
        high_resolution_policy,
        high_resolution_eps
    )

    bootstrap_mode <- .resolve_prism_bootstrap_mode(bootstrap)
    bootstrap_active <- bootstrap_mode != "none"

    if (is.null(counts)) stop("counts is NULL.")
    Y <- as.matrix(counts)
    if (length(dim(Y)) != 2) stop("counts must be a 2D object (matrix-like).")
    if (!is.numeric(Y)) stop("counts must be numeric.")
    if (any(!is.finite(Y))) stop("counts contains NA/NaN/Inf.")
    if (any(Y < 0)) stop("counts must be nonnegative.")

    tol_int <- 1e-8
    if (identical(counts_input, "proportions")) {
        if (!identical(composition, "none")) {
            stop("counts_input='proportions' is only valid with composition='none'.")
        }
        if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha != 0) {
            stop("counts_input='proportions' requires alpha=0 so observed proportions are not altered.")
        }
        cs <- colSums(Y)
        if (any(!is.finite(cs)) || any(cs <= 0)) {
            stop("counts_input='proportions' requires positive finite column sums.")
        }
        Y <- sweep(Y, 2, cs, "/")
    } else if (any(abs(Y - round(Y)) > tol_int)) {
        stop(
            "counts must contain integer-valued raw counts; use counts_input='proportions' only with a fixed composition estimator.",
            call. = FALSE
        )
    }

    if (length(dirichlet_concentration) != 1L || !is.numeric(dirichlet_concentration) || !is.finite(dirichlet_concentration) || dirichlet_concentration <= 0) {
        stop("dirichlet_concentration must be a single positive numeric.", call. = FALSE)
    }
    dirichlet_concentration <- as.numeric(dirichlet_concentration)

    storage.mode(Y) <- "double"

    mln_fit_object <- NULL
    mln_options <- NULL

    if (length(prevalence) != 1 || !is.numeric(prevalence) || prevalence < 0 || prevalence > 1) {
        stop("prevalence must be a single numeric value in [0, 1].")
    }

    if (prevalence > 0) {
        freq_present <- rowMeans(Y > 0)
        keep_taxa <- freq_present >= prevalence
        if (sum(keep_taxa) < 2) {
            stop(sprintf(
                "After prevalence filtering (threshold %g), fewer than 2 taxa remain.",
                prevalence
            ))
        }
        if (sum(!keep_taxa) > 0) {
            if (verbose) {
                message(sprintf(
                    "[prism] Dropped %d taxa below prevalence threshold %g",
                    sum(!keep_taxa), prevalence
                ))
            }
            Y <- Y[keep_taxa, , drop = FALSE]
        }
    }

    D <- nrow(Y)
    N <- ncol(Y)
    if (D < 1) stop("Need at least 1 feature (nrow(counts) >= 1).")
    if (N < 2) stop("Need at least 2 samples (ncol(counts) >= 2).")

    sample_names <- colnames(Y)
    if (is.null(sample_names) || any(sample_names == "")) {
        sample_names <- paste0("col_", seq_len(N))
        colnames(Y) <- sample_names
    }

    # Drop bad columns/rows
    col_totals <- colSums(Y)
    bad_cols <- which(!is.finite(col_totals) | col_totals <= 0)
    if (length(bad_cols) > 0) {
        bad_labels <- paste0(sample_names[bad_cols], " (index ", bad_cols, ")")
        warning(
            "Dropping ", length(bad_cols), " sample(s) with column sum <= 0: ",
            paste(bad_labels, collapse = ", ")
        )
        Y <- Y[, -bad_cols, drop = FALSE]
        if (!is.null(scale) && is.vector(scale)) {
            scale <- scale[-bad_cols]
        } else if (!is.null(scale)) {
            scale_dim <- dim(as.matrix(scale))
            if (scale_dim[1] == length(sample_names)) {
                scale <- as.matrix(scale)[-bad_cols, , drop = FALSE]
            } else if (scale_dim[2] == length(sample_names)) {
                scale <- as.matrix(scale)[, -bad_cols, drop = FALSE]
            }
        }
        if (!is.null(scale_log)) scale_log <- scale_log[-bad_cols]
        N <- ncol(Y)
        if (N < 2) stop("After dropping samples with column sum <= 0, need >= 2 samples.")
    }

    row_totals <- rowSums(Y)
    bad_rows <- which(!is.finite(row_totals) | row_totals <= 0)
    if (length(bad_rows) > 0) {
        feat_names <- rownames(Y)
        if (is.null(feat_names) || any(feat_names == "")) {
            feat_names <- paste0("row_", seq_len(nrow(Y)))
        }
        warning(
            "Dropping ", length(bad_rows), " feature(s) with row sum <= 0: ",
            paste(paste0(feat_names[bad_rows], " (index ", bad_rows, ")"), collapse = ", ")
        )
        Y <- Y[-bad_rows, , drop = FALSE]
        D <- nrow(Y)
        if (D < 2) stop("After dropping zero-sum features, need >= 2 features.")
    }

    high_resolution_pair_index <- .validate_high_resolution_pair_index(
        high_resolution_pair_index,
        D
    )
    if (!is.null(high_resolution_pair_index) && !high_resolution_config$enabled) {
        stop(
            "high_resolution_pair_index requires high_resolution=TRUE and a policy other than 'never'.",
            call. = FALSE
        )
    }

    # ---- dataset diagnostics (on filtered Y, printed before composition step) ---
    {
        lib_sizes <- colSums(Y)
        taxon_prevalence <- rowMeans(Y > 0) # fraction of samples where each taxon is present
        sparsity <- mean(Y == 0)

        nonzero_vals <- Y[Y > 0]
        mean_nz <- if (length(nonzero_vals) > 0) mean(nonzero_vals) else NA_real_
        median_nz <- if (length(nonzero_vals) > 0) stats::median(nonzero_vals) else NA_real_

        lib_mean <- mean(lib_sizes)
        lib_sd <- stats::sd(lib_sizes)
        lib_cv <- if (is.finite(lib_mean) && lib_mean > 0 && is.finite(lib_sd) && lib_sd > 0) {
            lib_sd / lib_mean
        } else {
            0
        }

        n_singletons <- sum(rowSums(Y > 0) == 1) # taxa present in exactly 1 sample

        data_diagnostics <- list(
            D                    = D,
            N                    = N,
            sparsity             = sparsity,
            D_over_N             = D / N,
            lib_size_min         = min(lib_sizes),
            lib_size_med         = stats::median(lib_sizes),
            lib_size_max         = max(lib_sizes),
            lib_size_cv          = lib_cv,
            prevalence_min       = min(taxon_prevalence),
            prevalence_med       = stats::median(taxon_prevalence),
            prevalence_max       = max(taxon_prevalence),
            mean_nonzero_count   = mean_nz,
            median_nonzero_count = median_nz,
            n_singleton_taxa     = n_singletons
        )

        msg <- sprintf(
            paste0(
                "[prism] Diagnostics: D=%d, N=%d, sparsity=%.1f%%, D/N=%.2f | ",
                "lib_size min/med/max=[%g/%g/%g] CV=%.2f | ",
                "taxon prevalence min/med/max=[%.2f/%.2f/%.2f] | ",
                "mean_nz=%.1f, singleton_taxa=%d"
            ),
            D, N, 100 * sparsity, D / N,
            min(lib_sizes), stats::median(lib_sizes), max(lib_sizes), lib_cv,
            min(taxon_prevalence), stats::median(taxon_prevalence), max(taxon_prevalence),
            mean_nz, n_singletons
        )
        if (isTRUE(verbose)) message(msg)

        if (verbose) {
            .fmt_row <- function(label, value, note = "") {
                sprintf("  %-28s %12s   %s", label, value, note)
            }

            msg_lines <- c(
                "[prism] Dataset summary table:",
                "  --------------------------------------------------------",
                .fmt_row("Metric", "Value", "Note"),
                "  --------------------------------------------------------",
                .fmt_row(
                    "Taxa (D)", formatC(D, format = "d"),
                    if (D / N > 1) "! D > N (MLN may struggle)" else ""
                ),
                .fmt_row("Samples (N)", formatC(N, format = "d"), ""),
                .fmt_row(
                    "D / N ratio", sprintf("%.2f", D / N),
                    if (D / N > 2) {
                        "!! high"
                    } else if (D / N > 1) {
                        "! moderate"
                    } else {
                        "ok"
                    }
                ),
                .fmt_row(
                    "Sparsity", sprintf("%.1f%%", 100 * sparsity),
                    if (sparsity > 0.9) {
                        "!! very high"
                    } else if (sparsity > 0.7) {
                        "! high"
                    } else {
                        ""
                    }
                ),
                .fmt_row(
                    "Singleton taxa", formatC(n_singletons, format = "d"),
                    if (n_singletons > 0) "present in exactly 1 sample" else "none"
                ),
                .fmt_row("Lib size (min)", formatC(min(lib_sizes), format = "fg", digits = 4), ""),
                .fmt_row("Lib size (median)", formatC(stats::median(lib_sizes), format = "fg", digits = 4), ""),
                .fmt_row("Lib size (max)", formatC(max(lib_sizes), format = "fg", digits = 4), ""),
                .fmt_row(
                    "Lib size CV", sprintf("%.2f", lib_cv),
                    if (lib_cv > 1) "! high variability in depth" else ""
                ),
                .fmt_row("Taxon prev (min)", sprintf("%.3f", min(taxon_prevalence)), ""),
                .fmt_row("Taxon prev (median)", sprintf("%.3f", stats::median(taxon_prevalence)), ""),
                .fmt_row("Taxon prev (max)", sprintf("%.3f", max(taxon_prevalence)), ""),
                .fmt_row(
                    "Mean non-zero count", sprintf("%.1f", mean_nz),
                    if (!is.na(mean_nz) && mean_nz < 5) "! very low counts" else ""
                ),
                .fmt_row(
                    "Median non-zero count", sprintf("%.1f", median_nz),
                    if (!is.na(median_nz) && median_nz < 3) "! very low typical non-zero count" else ""
                ),
                "  --------------------------------------------------------"
            )

            message(paste(msg_lines, collapse = "\n"))
        }
    }

    # ---- validate scalars ------------------------------------------------------
  if (length(S) != 1 || !is.finite(S) || S < 1) stop("S must be a single number >= 1.")
  S <- as.integer(S)
  if (length(max_attempts_per_draw) != 1L ||
      !is.numeric(max_attempts_per_draw) ||
      !is.finite(max_attempts_per_draw) ||
      max_attempts_per_draw < 1) {
    stop("max_attempts_per_draw must be one positive integer.", call. = FALSE)
  }
  max_attempts_per_draw <- as.integer(max_attempts_per_draw)
    if (length(alpha) != 1 || !is.finite(alpha) || alpha < 0) {
        stop("The composition estimator pseudocount must be >= 0.")
    }
    if (alpha == 0 && composition %in% c("none", "Dirichlet") && any(Y == 0)) {
        stop(
            "A zero pseudocount is invalid when zero counts would be logged or used as zero Dirichlet parameters.",
            call. = FALSE
        )
    }
    if (length(ci_alpha) != 1 || !is.finite(ci_alpha) || !(ci_alpha > 0 && ci_alpha < 1)) {
        stop("ci_alpha must be a single number in (0, 1).")
    }
    if (length(delta) != 1 || !is.finite(delta) || delta < 0) stop("delta must be a single finite number >= 0.")
    if (!is.character(p_adjust_method) || length(p_adjust_method) != 1) stop("p_adjust_method must be a single string.")
    if (!is.null(seed)) {
        if (length(seed) != 1 || !is.finite(seed)) stop("seed must be NULL or a single finite number.")
        set.seed(as.integer(seed))
    }
    if (length(scale_sd_fallback) != 1 || !is.numeric(scale_sd_fallback) ||
        !is.finite(scale_sd_fallback) || scale_sd_fallback < 0) {
        stop("scale_sd_fallback must be a single finite number >= 0.")
    }
    if (!is.numeric(scale_uncertainty_draws) || length(scale_uncertainty_draws) != 1L ||
        !is.finite(scale_uncertainty_draws) || scale_uncertainty_draws < 2) {
        stop("scale_uncertainty_draws must be one integer >= 2.", call. = FALSE)
    }
    scale_uncertainty_draws <- as.integer(scale_uncertainty_draws)
    for (flag_name in c("return_rel_draws", "return_bound_draws")) {
        flag_value <- get(flag_name, inherits = FALSE)
        if (!is.logical(flag_value) || length(flag_value) != 1L || is.na(flag_value)) {
            stop(flag_name, " must be a single TRUE or FALSE.", call. = FALSE)
        }
    }

    # ---- validate sigma / rho bounds -------------------------------------------
    if (!is.null(sigma_L)) {
        if (length(sigma_L) != 1 || !is.numeric(sigma_L) || !is.finite(sigma_L) || sigma_L < 0) {
            stop("sigma_L must be a single finite number >= 0 (or NULL).")
        }
    }
    if (!is.null(sigma_U)) {
        if (length(sigma_U) != 1 || !is.numeric(sigma_U) || !is.finite(sigma_U) || sigma_U < 0) {
            stop("sigma_U must be a single finite number >= 0 (or NULL).")
        }
    }
    if (!is.null(sigma_L) && !is.null(sigma_U) && sigma_L > sigma_U) {
        stop("Need sigma_L <= sigma_U when both are supplied.")
    }
    if (xor(is.null(sigma_L), is.null(sigma_U))) {
        stop("Provide both sigma_L and sigma_U, or neither.", call. = FALSE)
    }
    check_rho_vec <- function(x, name) {
        if (is.null(x)) {
            return(invisible(NULL))
        }
        if (!is.numeric(x) || any(!is.finite(x))) stop(paste0(name, " must be numeric and finite (or NULL)."))
        if (!(length(x) %in% c(1, D))) stop(paste0(name, " must have length 1 or D = ", D, " (or NULL)."))
        if (any(x < -1 | x > 1)) stop(paste0(name, " entries must be in [-1, 1]."))
        invisible(NULL)
    }
    check_rho_vec(rho_L, "rho_L")
    check_rho_vec(rho_U, "rho_U")
    if (xor(is.null(rho_L), is.null(rho_U))) {
        stop("Provide both rho_L and rho_U, or neither.", call. = FALSE)
    }
    if (!is.null(rho_L) && !is.null(rho_U)) {
        rL_chk <- if (length(rho_L) == 1) rep(rho_L, D) else as.numeric(rho_L)
        rU_chk <- if (length(rho_U) == 1) rep(rho_U, D) else as.numeric(rho_U)
        if (any(rL_chk > rU_chk)) stop("Need rho_L[d] <= rho_U[d] for all features d.")
    }
    if ((!is.null(rho_L) || !is.null(rho_U)) &&
        is.null(sigma_L) && is.null(sigma_U) &&
        is.null(scale) && is.null(scale_log)) {
        stop(
            "Correlation bounds require scale-SD bounds: provide sigma_L and sigma_U with rho_L and rho_U.",
            call. = FALSE
        )
    }

    # ---- validate scale --------------------------------------------------------
    has_scale <- !is.null(scale) && length(scale) > 0
    scale_mat_full <- NULL
    if (has_scale) {
        if (is.vector(scale)) {
            if (!is.numeric(scale) || any(!is.finite(scale))) stop("scale must be numeric with finite values.")
            if (length(scale) != N) stop(paste0("scale must have length N = ", N, "."))
            if (any(scale <= 0)) stop("scale must be strictly positive.")
            scale_mat_full <- matrix(as.numeric(scale), nrow = N, ncol = 1)
        } else {
            sm0 <- as.matrix(scale)
            if (!is.numeric(sm0)) stop("scale matrix must be numeric.")
            if (nrow(sm0) == N) {
                scale_mat_full <- sm0
            } else if (ncol(sm0) == N) {
                scale_mat_full <- t(sm0)
            } else {
                stop("scale matrix must have nrow(scale)=N or ncol(scale)=N.")
            }
            if (all(!is.finite(scale_mat_full))) stop("scale has no finite measurements.")
            if (any(is.finite(scale_mat_full) & scale_mat_full <= 0)) stop("scale must be > 0 wherever finite.")
        }
        if (!is.null(sigma_L) || !is.null(sigma_U) || !is.null(rho_L) || !is.null(rho_U)) {
            warning("scale was supplied, so sigma_L/sigma_U and rho_L/rho_U are ignored.")
        }
        sigma_L <- NULL
        sigma_U <- NULL
        rho_L <- NULL
        rho_U <- NULL
    }

    # ---- validate scale_log --------------------------------------------------
    # scale_log takes priority over sigma_L/sigma_U/rho_L/rho_U but yields to
    # scale (raw measurements). When scale is supplied, scale_log is ignored.
    has_scale_log <- !is.null(scale_log) && !has_scale
    if (!is.null(scale_log) && has_scale) {
        warning("Both 'scale' and 'scale_log' supplied; 'scale_log' is ignored (scale takes priority).")
        scale_log <- NULL
        has_scale_log <- FALSE
    }
    if (has_scale_log) {
        if (!is.numeric(scale_log) || length(scale_log) != N) {
            stop("scale_log must be a numeric vector of length N = ", N, ".")
        }
        if (any(!is.finite(scale_log))) {
            stop("scale_log must be finite (no NA/Inf). Log-transform the raw load values before passing.")
        }
        if (!is.null(sigma_L) || !is.null(sigma_U) || !is.null(rho_L) || !is.null(rho_U)) {
            warning("scale_log supplied; sigma_L/sigma_U/rho_L/rho_U are ignored.")
        }
        sigma_L <- NULL
        sigma_U <- NULL
        rho_L <- NULL
        rho_U <- NULL
        if (length(scale_log_ci_level) != 1 || !is.finite(scale_log_ci_level) ||
            scale_log_ci_level < 0 || scale_log_ci_level >= 1) {
            stop("scale_log_ci_level must be a single number in [0, 1).")
        }
    }
    corr_safe <- function(x, u) {
        xc <- x - mean(x)
        uc <- u - mean(u)
        num <- sum(xc * uc)
        den <- sqrt(sum(xc * xc) * sum(uc * uc))
        if (!is.finite(den) || den <= 0) {
            return(0)
        }
        r <- num / den
        if (!is.finite(r)) r <- 0
        max(-1, min(1, r))
    }

    .resolve_prism_stream <- function(
      stream,
      D,
      S,
      default_block_size = 50000L,
      auto_pair_draw_threshold = 2e6
    ) {
        n_off <- D * (D - 1) / 2

        if (is.null(stream)) {
            use_stream <- (n_off * S) > auto_pair_draw_threshold
            block_size <- as.integer(default_block_size)
            mode <- "auto"
        } else if (is.logical(stream) && length(stream) == 1L && !is.na(stream)) {
            use_stream <- isTRUE(stream)
            block_size <- as.integer(default_block_size)
            mode <- if (use_stream) "forced" else "full"
        } else if (is.numeric(stream) && length(stream) == 1L && is.finite(stream) && stream > 0) {
            use_stream <- TRUE
            block_size <- as.integer(max(1L, floor(stream)))
            mode <- "forced_numeric"
        } else {
            stop(
                "stream must be NULL, TRUE, FALSE, or a positive numeric block size.",
                call. = FALSE
            )
        }

        list(
            use_stream = use_stream,
            mode = mode,
            block_size = block_size,
            n_offdiag_pairs = as.integer(n_off),
            n_diag_pairs = as.integer(D),
            auto_pair_draw_threshold = auto_pair_draw_threshold
        )
    }

    .make_prism_stream_index <- function(D) {
        off <- which(upper.tri(matrix(0, D, D), diag = FALSE), arr.ind = TRUE)

        # Diagonal entries are always included in streamed mode because
        # pairwise_rel needs i == i rows for relative variances/SDs/correlations.
        diag_idx <- cbind(seq_len(D), seq_len(D))

        list(
            off = off,
            diag = diag_idx,
            K_off = nrow(off),
            K_diag = nrow(diag_idx)
        )
    }

    alpha_num <- as.numeric(alpha)
    q_lo <- ci_alpha / 2
    q_hi <- 1 - ci_alpha / 2

    # Pair indices (upper triangle)
    idx_pairs_rel <- which(upper.tri(matrix(0, D, D), diag = TRUE), arr.ind = TRUE)
    i_rel <- idx_pairs_rel[, 1]
    j_rel <- idx_pairs_rel[, 2]
    n_pairs_rel <- nrow(idx_pairs_rel)
    col_idx_rel <- (j_rel - 1L) * D + i_rel

    make_pairwise_bound_draws <- function(draw_id, i, j, lower, upper) {
        if (!isTRUE(return_bound_draws)) return(NULL)
        pair_id <- paste(i, j, sep = "__")
        dimnames(lower) <- list(as.character(draw_id), pair_id)
        dimnames(upper) <- list(as.character(draw_id), pair_id)
        feature_names <- rownames(Y)
        pair_index <- data.frame(
            pair_id = pair_id,
            i = as.integer(i),
            j = as.integer(j),
            taxon_i = if (is.null(feature_names)) NA_character_ else feature_names[i],
            taxon_j = if (is.null(feature_names)) NA_character_ else feature_names[j],
            stringsAsFactors = FALSE
        )
        list(
            draw_id = as.integer(draw_id),
            pair_index = pair_index,
            lower = lower,
            upper = upper
        )
    }

    stream_cfg <- .resolve_prism_stream(stream = stream, D = D, S = S)

    # ===========================================================================
    # FAST PATH: no composition noise AND no bootstrap => single deterministic draw
    # ===========================================================================
    deterministic <- (composition == "none" && bootstrap_mode == "none")

    if (deterministic) {
        Y0 <- Y + alpha_num
        P0 <- sweep(Y0, 2, colSums(Y0), "/")
        if (any(!is.finite(P0)) || any(P0 <= 0)) stop("Deterministic mode: failed to form strictly-positive compositions.")
        logP0 <- log(P0)
        Sigma_rel <- stats::cov(t(logP0))
        rho_witness_use <- NULL

        # Scale correction (use as-is if scale data supplied, otherwise take supplied bounds)
        if (has_scale) {
            u_full <- log(scale_mat_full)
            mu_u_full <- apply(u_full, 1, function(v) mean(v, na.rm = TRUE))
            sigma_hat_sd <- stats::sd(mu_u_full)
            if (!is.finite(sigma_hat_sd)) sigma_hat_sd <- 0
            rho_hat <- numeric(D)
            for (d in seq_len(D)) rho_hat[d] <- corr_safe(logP0[d, ], mu_u_full)
            sL_use <- sigma_hat_sd
            sU_use <- sigma_hat_sd
            rL_use <- rho_hat
            rU_use <- rho_hat
        } else if (has_scale_log) {
            scale_bounds <- .estimate_scale_log_bounds(
                log_proportions = logP0,
                scale_log = scale_log,
                ci_level = scale_log_ci_level
            )
            sL_use <- scale_bounds$sigma_L
            sU_use <- scale_bounds$sigma_U
            if (isTRUE(scale_log_rho)) {
                rL_use <- scale_bounds$rho_L
                rU_use <- scale_bounds$rho_U
                rho_witness_use <- scale_bounds$rho_witness
            } else {
                rL_use <- NULL
                rU_use <- NULL
                rho_witness_use <- NULL
            }
        } else {
            sL_use <- sigma_L
            sU_use <- sigma_U
            rL_use <- rho_L
            rU_use <- rho_U
        }

        point <- cov_bounds(
            Sigma_rel,
            sigma_L = sL_use,
            sigma_U = sU_use,
            rho_L = rL_use,
            rho_U = rU_use,
            rho_witness = rho_witness_use,
            solver = solver,
            eig_tol = high_resolution_eig_tol
        )
        ci_lower <- point$lower
        ci_upper <- point$upper
        idx_pairs <- which(upper.tri(ci_lower), arr.ind = TRUE)
        i_vec <- idx_pairs[, 1]
        j_vec <- idx_pairs[, 2]

        pairwise <- data.frame(
            i = i_vec, j = j_vec,
            ci_lower = ci_lower[idx_pairs],
            ci_upper = ci_upper[idx_pairs],
            p_value = NA_real_,
            q_value = NA_real_,
            delta = delta
        )
        high_resolution_diagnostics <- NULL
        if (high_resolution_config$enabled) {
            if (point$regime != "bounded_scale_and_correlation") {
                warning("high_resolution=TRUE has no effect outside bounded rho interval mode.", call. = FALSE)
            } else if (all(point$inputs$rho_L == point$inputs$rho_U)) {
                warning("high_resolution=TRUE is unnecessary when rho_L == rho_U.", call. = FALSE)
            } else {
                refine <- .select_requested_high_resolution_pairs(
                    cbind(pairwise$i, pairwise$j),
                    high_resolution_pair_index,
                    pairwise$ci_lower,
                    pairwise$ci_upper,
                    delta,
                    high_resolution_config$policy,
                    high_resolution_config$eps
                )
                if (any(refine)) {
                    standard_lower <- pairwise$ci_lower[refine]
                    standard_upper <- pairwise$ci_upper[refine]
                    sharp <- cov_bounds_pairs(
                        Sigma_rel,
                        pair_index = cbind(pairwise$i[refine], pairwise$j[refine]),
                        sigma_L = sL_use,
                        sigma_U = sU_use,
                        rho_L = rL_use,
                        rho_U = rU_use,
                        rho_witness = rho_witness_use,
                        high_resolution = TRUE,
                        solver = solver,
                        eig_tol = high_resolution_eig_tol
                    )
                    tightened <- .validate_high_resolution_tightening(
                        standard_lower,
                        standard_upper,
                        sharp$lower,
                        sharp$upper
                    )
                    pairwise$ci_lower[refine] <- tightened$lower
                    pairwise$ci_upper[refine] <- tightened$upper
                    for (k in which(refine)) {
                        ii <- pairwise$i[k]
                        jj <- pairwise$j[k]
                        ci_lower[ii, jj] <- ci_lower[jj, ii] <- pairwise$ci_lower[k]
                        ci_upper[ii, jj] <- ci_upper[jj, ii] <- pairwise$ci_upper[k]
                    }
                    high_resolution_diagnostics <- .high_resolution_ci_diagnostics(
                        pairwise$i[refine],
                        pairwise$j[refine],
                        standard_lower,
                        standard_upper,
                        pairwise$ci_lower[refine],
                        pairwise$ci_upper[refine]
                    )
                }
            }
        }
        rel_hat <- Sigma_rel[col_idx_rel]
        pairwise_rel <- data.frame(
            i = i_rel, j = j_rel,
            rel_hat = rel_hat,
            rel_ci_lower = rel_hat,
            rel_ci_upper = rel_hat,
            rel_ci_width = rep(0, length(rel_hat))
        )
        pairwise_bound_draws <- make_pairwise_bound_draws(
            draw_id = 1L,
            i = i_vec,
            j = j_vec,
            lower = matrix(pairwise$ci_lower, nrow = 1L),
            upper = matrix(pairwise$ci_upper, nrow = 1L)
        )

        return(.finalize_prism_result(list(
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            point_lower = ci_lower,
            point_upper = ci_upper,
            pairwise = pairwise,
            pairwise_rel = pairwise_rel,
            high_resolution_diagnostics = high_resolution_diagnostics,
            rel_draws = NULL,
            pairwise_bound_draws = pairwise_bound_draws,
            mln_fit_summary = NULL,
            mln_diagnostics = NULL,
            data_diagnostics = data_diagnostics,
            sampling_diagnostics = list(
                bootstrap_mode = bootstrap_mode,
                requested_draws = 1L,
                accepted_draws = 1L,
                attempted_draws = 1L,
                rejected_attempts = 0L,
                acceptance_fraction = 1,
                rejection_fraction = 0,
                draws_requiring_replacement = 0L,
                attempts_per_requested_draw = c(min = 1, median = 1, max = 1),
                rejection_reasons = integer(),
                rejection_reason_fractions = numeric(),
                max_attempts_per_draw = 1L
            ),
            params = list(
                composition = composition,
                bootstrap = bootstrap_active,
                bootstrap_mode = bootstrap_mode,
                S = 1L, alpha = alpha, delta = delta,
                dirichlet_concentration = dirichlet_concentration,
                composition_storage = "fixed_matrix",
                stream = stream,
                workers = workers,
                return_bound_draws = return_bound_draws,
                stream_active = FALSE,
                scale_summary = if (has_scale || has_scale_log) list(
                    sigma_L = sL_use,
                    sigma_U = sU_use,
                    rho_L = rL_use,
                    rho_U = rU_use
                ) else NULL,
                high_resolution = high_resolution_config$enabled,
                high_resolution_policy = high_resolution_config$policy,
                high_resolution_eps = high_resolution_config$eps,
                high_resolution_selection = if (is.null(high_resolution_pair_index)) "policy" else "requested_pairs",
                high_resolution_requested_pairs = if (is.null(high_resolution_pair_index)) 0L else nrow(high_resolution_pair_index),
                solver = solver
            )
        ), composition_estimator))
    }

    # ===========================================================================
    # STOCHASTIC PATH: S draws, mix-and-match composition noise + bootstrap
    # ===========================================================================

    mln_fit_summary <- NULL
    mln_diagnostics <- NULL
    composition_rng_seed <- rep(NA_integer_, S)

    # -- pre-compute composition draws if needed ---------------------------------
    if (composition == "none") {
        Y0 <- Y + alpha_num
        P_fixed <- sweep(Y0, 2, colSums(Y0), "/")
        if (any(!is.finite(P_fixed)) || any(P_fixed <= 0)) stop("composition='none': failed to form strictly-positive compositions.")
        logP_fixed <- log(P_fixed)
        Pi <- NULL
    } else if (composition == "custom") {
        composition_fit_object <- composition_estimator$options$fit(
            counts = Y,
            n_draws = S,
            seed = seed,
            verbose = verbose
        )
        Pi <- NULL
        mln_diagnostics <- list(estimator = composition_estimator$name)
    } else if (composition == "MLN") {
        mln_options <- .resolve_mln_options(
            composition_estimator$options[setdiff(names(composition_estimator$options), "pseudocount")],
            .default_mln_options()
        )
        mln_fit <- .fit_composition_mln(
            Y = Y,
            n_draws = S,
            pseudocount = alpha_num,
            ci_alpha = ci_alpha,
            options = mln_options,
            verbose = verbose
        )
        Pi <- mln_fit$draws
        S <- dim(Pi)[3]
        mln_fit_object <- mln_fit$fit
        mln_fit_summary <- mln_fit$fit_summary
        mln_diagnostics <- mln_fit$diagnostics
    } else {
        # Draw on demand so streamed runs do not retain a D x N x S array.
        composition_rng_seed <- sample.int(.Machine$integer.max, S)
        Pi <- NULL
    }

    composition_storage <- switch(
        composition,
        Dirichlet = "draw_on_demand",
        none = "fixed_matrix",
        custom = "plugin_managed",
        MLN = "posterior_draw_array"
    )

    # Re-resolve stream after composition fitting because MLN paths may adjust S.
    stream_cfg <- .resolve_prism_stream(stream = stream, D = D, S = S)

    if (length(workers) != 1L || !is.numeric(workers) || !is.finite(workers) || workers < 1) {
        stop("workers must be a single positive integer.", call. = FALSE)
    }
    workers <- as.integer(floor(workers))
    workers <- max(1L, min(workers, S))

    if (workers > 1L && !isTRUE(stream_cfg$use_stream)) {
        warning(
            "workers > 1 currently applies only when streamed mode is active. ",
            "Using workers = 1 because stream resolved to full-matrix mode. ",
            "Set stream = TRUE or stream = a positive block size to enable internal parallel draw processing.",
            call. = FALSE
        )
        workers <- 1L
    }

    # -- allocate draw storage after S is final ----------------------------------
    if (isTRUE(stream_cfg$use_stream)) {
        if (!exists("cov_bounds_pairs", mode = "function")) {
            stop(
                "stream=TRUE requires cov_bounds_pairs(). ",
                "Source the patched bounds.R before source(prism.R).",
                call. = FALSE
            )
        }

        stream_index <- .make_prism_stream_index(D = D)

        i_off <- stream_index$off[, 1]
        j_off <- stream_index$off[, 2]
        i_diag <- stream_index$diag[, 1]
        j_diag <- stream_index$diag[, 2]

        K_off <- stream_index$K_off
        K_diag <- stream_index$K_diag

        lower_arr <- NULL
        upper_arr <- NULL
        rel_draws <- NULL
    } else {
        lower_arr <- array(NA_real_, dim = c(D, D, S))
        upper_arr <- array(NA_real_, dim = c(D, D, S))
        rel_draws <- matrix(NA_real_, nrow = S, ncol = n_pairs_rel)
    }

    valid_draw <- rep(FALSE, S)
    skip_reason <- rep(NA_character_, S)
    draw_attempts <- integer(S)
    draw_rejection_reasons <- vector("list", S)

    # tracking vectors for scale_log per-draw summaries
    sigma_L_track <- rep(NA_real_, S)
    sigma_U_track <- rep(NA_real_, S)
    rho_L_track <- rep(NA_real_, S)
    rho_U_track <- rep(NA_real_, S)
    Sigma_rel_high_resolution <- if (high_resolution_config$enabled) vector("list", S) else NULL
    sigma_L_high_resolution <- if (high_resolution_config$enabled) rep(NA_real_, S) else NULL
    sigma_U_high_resolution <- if (high_resolution_config$enabled) rep(NA_real_, S) else NULL
    rho_L_high_resolution <- if (high_resolution_config$enabled) vector("list", S) else NULL
    rho_U_high_resolution <- if (high_resolution_config$enabled) vector("list", S) else NULL
    rho_witness_high_resolution <- if (high_resolution_config$enabled) vector("list", S) else NULL

    # Precompute retry seeds so serial and parallel streamed runs use the same
    # bootstrap samples. Invalid hybrid draws may be replaced, up to this
    # finite cap, without changing the requested number of accepted draws.
  bootstrap_rng_seed <- if (bootstrap_active) {
        matrix(
            sample.int(.Machine$integer.max, size = S * max_attempts_per_draw),
            nrow = S,
            ncol = max_attempts_per_draw
        )
    } else {
        matrix(NA_integer_, nrow = S, ncol = 1L)
    }

    # The scale replicate branch simulates auxiliary load draws inside each Monte
    # Carlo draw. Precompute one seed per draw so workers are schedule-invariant.
    scale_rng_seed <- if (has_scale) {
        sample.int(.Machine$integer.max, size = S)
    } else {
        rep(NA_integer_, S)
    }

    .bad_bounds_vec <- function(lower, upper) {
        if (is.null(lower) || is.null(upper)) {
            return(TRUE)
        }
        if (length(lower) != length(upper)) {
            return(TRUE)
        }

        # Lower bounds should be finite. Upper bounds may be +Inf in the
        # unbounded-scale case, but should not be NA, NaN, or -Inf.
        bad_lower <- any(is.na(lower) | is.nan(lower) | is.infinite(lower))
        bad_upper <- any(is.na(upper) | is.nan(upper) | (is.infinite(upper) & upper < 0))

        bad_lower || bad_upper
    }

    .draw_fail <- function(s, reason) {
        list(
            s = s,
            valid = FALSE,
            skip_reason = reason,
            sigma_L_track = NA_real_,
            sigma_U_track = NA_real_,
            rho_L_track = NA_real_,
            rho_U_track = NA_real_,
            rel_off = NULL,
            lower_off = NULL,
            upper_off = NULL,
            rel_diag = NULL,
            lower_diag = NULL,
            upper_diag = NULL,
            rel_full = NULL,
            lower_full = NULL,
            upper_full = NULL
        )
    }

    .process_one_attempt <- function(s, attempt, stream_mode, high_resolution_pairs = NULL) {
        idx_boot <- seq_len(N)
        if (bootstrap_active) {
            set.seed(bootstrap_rng_seed[s, attempt])
            idx_boot <- sample.int(N, size = N, replace = TRUE)
        }
        idx_composition <- idx_boot
        idx_scale <- idx_boot

        # 1. COMPOSITION: get log-proportions for draw s after bootstrap resampling.
        if (composition == "none") {
            logP_base <- logP_fixed
        } else {
            Pi_s <- if (composition == "custom") {
                draw <- composition_estimator$options$draw(
                    fit = composition_fit_object,
                    draw = s
                )
                .validate_composition_draw(
                    draw,
                    D = D,
                    N = N,
                    estimator_name = composition_estimator$name
                )
            } else if (composition == "Dirichlet") {
                set.seed(composition_rng_seed[s])
                draw <- matrix(NA_real_, nrow = D, ncol = N)
                for (n in seq_len(N)) {
                    shape <- dirichlet_concentration * (Y[, n] + alpha_num)
                    gamma <- stats::rgamma(D, shape = shape, rate = 1)
                    draw[, n] <- gamma / sum(gamma)
                }
                draw
            } else {
                Pi[, , s]
            }
            Pi_boot <- Pi_s

            if (any(!is.finite(Pi_boot))) {
                return(.draw_fail(s, "non_finite_Pi_boot"))
            }

            Pi_boot <- pmax(Pi_boot, .Machine$double.xmin)
            cs <- colSums(Pi_boot)
            if (any(!is.finite(cs)) || any(cs <= 0)) {
                return(.draw_fail(s, "bad_column_sums_after_clip"))
            }
            Pi_boot <- sweep(Pi_boot, 2, cs, "/")

            if (any(!is.finite(Pi_boot)) || any(Pi_boot <= 0)) {
                return(.draw_fail(s, "bad_Pi_boot_after_renorm"))
            }

            logP_base <- log(Pi_boot)
        }

        logP_s <- logP_base[, idx_composition, drop = FALSE]

        if (any(!is.finite(logP_s))) {
            return(.draw_fail(s, "non_finite_logP"))
        }

    Sigma_rel_s <- tryCatch(stats::cov(t(logP_s)), error = function(e) NULL)
    if (is.null(Sigma_rel_s) || any(!is.finite(Sigma_rel_s))) {
      return(.draw_fail(s, "bad_Sigma_rel"))
    }
    variance_scale <- max(1, max(abs(diag(Sigma_rel_s))))
    variance_tol <- 100 * .Machine$double.eps * variance_scale
    if (any(diag(Sigma_rel_s) <= variance_tol)) {
      return(.draw_fail(s, "zero_variance_relative_feature"))
    }

        if (isTRUE(stream_mode)) {
            rel_off_s <- if (K_off > 0L) Sigma_rel_s[cbind(i_off, j_off)] else numeric(0)
            rel_diag_s <- if (K_diag > 0L) Sigma_rel_s[cbind(i_diag, j_diag)] else numeric(0)
            rel_full_s <- NULL
        } else {
            rel_off_s <- NULL
            rel_diag_s <- NULL
            rel_full_s <- Sigma_rel_s[col_idx_rel]
        }

        # 2. SCALE: determine sigma and rho for this draw.
        sigma_L_s <- sigma_L
        sigma_U_s <- sigma_U
        rho_L_s <- rho_L
        rho_U_s <- rho_U
        rho_witness_s <- NULL

        sigma_L_track_s <- NA_real_
        sigma_U_track_s <- NA_real_
        rho_L_track_s <- NA_real_
        rho_U_track_s <- NA_real_

        # 2a. scale_log branch: bootstrap-aware per-draw sigma and rho estimation.
        if (has_scale_log) {
            u_sigma <- scale_log[idx_scale]
            n_boot <- length(u_sigma)

            s_boot <- stats::sd(u_sigma)
            if (!is.finite(s_boot) || s_boot < 0) s_boot <- 0

            if (scale_log_ci_level == 0 || n_boot < 2) {
                sigma_L_s <- s_boot
                sigma_U_s <- s_boot
            } else {
                df_boot <- n_boot - 1L
                alpha_boot <- 1 - scale_log_ci_level
                sigma_L_s <- sqrt(df_boot * s_boot^2 / stats::qchisq(1 - alpha_boot / 2, df = df_boot))
                sigma_U_s <- sqrt(df_boot * s_boot^2 / stats::qchisq(alpha_boot / 2, df = df_boot))
            }
            # Asymmetric scale-CI variant: sigma is a standard deviation, so it
            # cannot be negative -- the two-sided chi-square lower quantile
            # above asserts a positive floor that isn't a real constraint on
            # sigma itself, only an artifact of splitting alpha_boot evenly
            # across both tails. Forcing it to 0 keeps sigma_U_s (and
            # everything downstream that depends on it) identical to the
            # ordinary two-sided regime at the same scale_log_ci_level, while
            # removing that lower-tail assertion entirely.
            if (isTRUE(scale_log_lower_zero)) sigma_L_s <- 0

            if (isTRUE(scale_log_rho)) {
            u_rho <- scale_log[idx_scale]
            logP_rho <- logP_base[, idx_scale, drop = FALSE]
            u_c <- u_rho - mean(u_rho)
            sdu <- sqrt(sum(u_c^2) / (n_boot - 1))
            Xc_s <- logP_rho - rowMeans(logP_rho, na.rm = TRUE)

            if (!is.finite(sdu) || sdu <= 0 || n_boot < 3) {
                r_boot <- rep(0, D)
            } else {
                sd_x_s <- sqrt(rowSums(Xc_s^2, na.rm = TRUE) / (n_boot - 1))
                sd_x_s[!is.finite(sd_x_s) | sd_x_s <= 0] <- NA_real_
                cov_xu <- rowSums(Xc_s * rep(u_c, each = D), na.rm = TRUE) / (n_boot - 1)
                r_boot <- cov_xu / (sd_x_s * sdu)
                r_boot[!is.finite(r_boot)] <- 0
                r_boot <- pmax(-1, pmin(1, r_boot))
            }

            if (scale_log_ci_level == 0 || n_boot < 4) {
                rho_L_s <- r_boot
                rho_U_s <- r_boot
            } else {
                alpha_boot <- 1 - scale_log_ci_level
                z_boot <- atanh(r_boot)
                z_boot[is.infinite(z_boot)] <- sign(z_boot[is.infinite(z_boot)]) * 5
                se_z_boot <- 1 / sqrt(n_boot - 3)
                me_z <- stats::qnorm(1 - alpha_boot / 2) * se_z_boot
                rho_L_s <- tanh(z_boot - me_z)
                rho_U_s <- tanh(z_boot + me_z)
                rho_L_s <- pmax(-1, pmin(1, rho_L_s))
                rho_U_s <- pmax(-1, pmin(1, rho_U_s))
                rho_L_s <- pmin(rho_L_s, r_boot)
                rho_U_s <- pmax(rho_U_s, r_boot)
            }
            rho_witness_s <- r_boot
            } else {
                rho_L_s <- NULL
                rho_U_s <- NULL
                rho_witness_s <- NULL
            }

            sigma_L_track_s <- sigma_L_s
            sigma_U_track_s <- sigma_U_s
            if (!is.null(rho_L_s)) rho_L_track_s <- mean(rho_L_s, na.rm = TRUE)
            if (!is.null(rho_U_s)) rho_U_track_s <- mean(rho_U_s, na.rm = TRUE)
        }

        # 2b. scale replicate branch.
        if (has_scale) {
            if (is.finite(scale_rng_seed[s])) {
                set.seed(scale_rng_seed[s])
            }

            scale_mat <- scale_mat_full[idx_scale, , drop = FALSE]
            if (all(!is.finite(scale_mat))) stop("scale has no finite measurements (after resampling).")
            if (any(is.finite(scale_mat) & scale_mat <= 0)) stop("scale must be > 0 wherever finite.")
            u_mat <- log(scale_mat)
            nrep <- rowSums(is.finite(u_mat))
            if (any(nrep == 0)) stop("Some resampled samples have no finite scale measurements.")
            mu_u <- apply(u_mat, 1, function(v) mean(v, na.rm = TRUE))
            sd_within <- apply(u_mat, 1, function(v) {
                vv <- v[is.finite(v)]
                if (length(vv) < 2) NA_real_ else stats::sd(vv)
            })

            idx_rep <- which(nrep >= 2 & is.finite(sd_within))
            if (length(idx_rep) == 0) {
                if (s == 1L) warning("No samples have >= 2 scale replicates; using scale_sd_fallback = ", scale_sd_fallback)
                pooled_within_sd <- scale_sd_fallback
            } else {
                pooled_within_sd <- sqrt(stats::median(sd_within[idx_rep]^2))
            }
            sd_within_filled <- sd_within
            sd_within_filled[!is.finite(sd_within_filled)] <- pooled_within_sd
            sd_mean <- sd_within_filled / sqrt(nrep)
            mu_u <- mu_u - mean(mu_u)

            u_draws <- matrix(
                stats::rnorm(
                    N * scale_uncertainty_draws,
                    mean = rep(mu_u, times = scale_uncertainty_draws),
                    sd = rep(sd_mean, times = scale_uncertainty_draws)
                ),
                nrow = N, ncol = scale_uncertainty_draws
            )
            Uc <- sweep(u_draws, 2, colMeans(u_draws), "-")
            sigma_draw <- sqrt(colSums(Uc^2) / (N - 1))
            sigma_L_s <- as.numeric(stats::quantile(sigma_draw, probs = q_lo, names = FALSE, type = 7))
            sigma_U_s <- as.numeric(stats::quantile(sigma_draw, probs = q_hi, names = FALSE, type = 7))

            logP_scale <- logP_base[, idx_scale, drop = FALSE]
            Xc <- logP_scale - rowMeans(logP_scale, na.rm = TRUE)
            sd_x <- sqrt(rowSums(Xc^2, na.rm = TRUE) / (N - 1))
            sd_x[!is.finite(sd_x) | sd_x <= 0] <- NA_real_
            rho_L_s <- numeric(D)
            rho_U_s <- numeric(D)
            joint_witness_draw <- rep(TRUE, scale_uncertainty_draws)

            chunk <- 2000L
            for (start in seq.int(1L, D, by = chunk)) {
                end <- min(D, start + chunk - 1L)
                idx_feat <- start:end
                rho_blk <- matrix(0, nrow = length(idx_feat), ncol = scale_uncertainty_draws)
                Xc_blk <- Xc[idx_feat, , drop = FALSE]
                sd_x_blk <- sd_x[idx_feat]
                for (ss in seq_len(scale_uncertainty_draws)) {
                    u <- u_draws[, ss]
                    uc <- u - mean(u)
                    sdu <- sqrt(sum(uc * uc) / (N - 1))
                    if (!is.finite(sdu) || sdu <= 0) {
                        rho_blk[, ss] <- 0
                    } else {
                        cov_xu <- rowSums(Xc_blk * rep(uc, each = nrow(Xc_blk))) / (N - 1)
                        denom <- sd_x_blk * sdu
                        r <- cov_xu / denom
                        r[!is.finite(r)] <- 0
                        r[r > 1] <- 1
                        r[r < -1] <- -1
                        rho_blk[, ss] <- r
                    }
                }
                rho_L_s[idx_feat] <- apply(rho_blk, 1, stats::quantile, probs = q_lo, names = FALSE, type = 7)
                rho_U_s[idx_feat] <- apply(rho_blk, 1, stats::quantile, probs = q_hi, names = FALSE, type = 7)
                joint_witness_draw <- joint_witness_draw & colSums(
                    rho_blk >= rho_L_s[idx_feat] & rho_blk <= rho_U_s[idx_feat]
                ) == length(idx_feat)
            }
            if (any(joint_witness_draw)) {
                witness_draw <- which(joint_witness_draw)[1L]
                u_witness <- u_draws[, witness_draw]
                u_witness <- u_witness - mean(u_witness)
                sd_u_witness <- sqrt(sum(u_witness^2) / (N - 1))
                if (is.finite(sd_u_witness) && sd_u_witness > 0) {
                    rho_witness_s <- rowSums(
                        Xc * rep(u_witness, each = D)
                    ) / ((N - 1) * sd_x * sd_u_witness)
                    rho_witness_s[!is.finite(rho_witness_s)] <- 0
                    rho_witness_s <- pmax(-1, pmin(1, rho_witness_s))
                }
            }
            sigma_L_track_s <- sigma_L_s
            sigma_U_track_s <- sigma_U_s
            rho_L_track_s <- mean(rho_L_s, na.rm = TRUE)
            rho_U_track_s <- mean(rho_U_s, na.rm = TRUE)
        }

        if (!is.null(high_resolution_pairs)) {
            if (is.null(sigma_L_s) || is.null(sigma_U_s) ||
                is.null(rho_L_s) || is.null(rho_U_s) ||
                all(rho_L_s == rho_U_s)) {
                return(list(s = s, valid = TRUE, sharp_applied = FALSE))
            }
            sharp <- cov_bounds_pairs(
                Sigma_rel = Sigma_rel_s,
                pair_index = high_resolution_pairs,
                sigma_L = sigma_L_s,
                sigma_U = sigma_U_s,
                rho_L = rho_L_s,
                rho_U = rho_U_s,
                rho_witness = rho_witness_s,
                high_resolution = TRUE,
                solver = solver,
                eig_tol = high_resolution_eig_tol
            )
            if (nrow(sharp) != nrow(high_resolution_pairs) ||
                .bad_bounds_vec(sharp$lower, sharp$upper)) {
                return(.draw_fail(s, "bad_high_resolution_bounds"))
            }
            return(list(
                s = s,
                valid = TRUE,
                sharp_applied = TRUE,
                lower = sharp$lower,
                upper = sharp$upper
            ))
        }

        # 3. COV_BOUNDS.
        if (isTRUE(stream_mode)) {
            pair_index_all <- rbind(
                cbind(i_off, j_off),
                cbind(i_diag, j_diag)
            )
            point_all <- cov_bounds_pairs(
                Sigma_rel = Sigma_rel_s,
                pair_index = pair_index_all,
                sigma_L = sigma_L_s,
                sigma_U = sigma_U_s,
                rho_L = rho_L_s,
                rho_U = rho_U_s,
                rho_witness = rho_witness_s,
                solver = solver,
                eig_tol = high_resolution_eig_tol
            )

            if (nrow(point_all) != K_off + K_diag ||
                .bad_bounds_vec(point_all$lower, point_all$upper)) {
                return(.draw_fail(s, "bad_cov_bounds"))
            }

            off_rows <- seq_len(K_off)
            diag_rows <- K_off + seq_len(K_diag)
            lower_off_s <- point_all$lower[off_rows]
            upper_off_s <- point_all$upper[off_rows]
            lower_diag_s <- point_all$lower[diag_rows]
            upper_diag_s <- point_all$upper[diag_rows]

            return(list(
                s = s,
                valid = TRUE,
                skip_reason = NA_character_,
                sigma_L_track = sigma_L_track_s,
                sigma_U_track = sigma_U_track_s,
                rho_L_track = rho_L_track_s,
                rho_U_track = rho_U_track_s,
                rel_off = rel_off_s,
                lower_off = lower_off_s,
                upper_off = upper_off_s,
                rel_diag = rel_diag_s,
                lower_diag = lower_diag_s,
                upper_diag = upper_diag_s,
                rel_full = NULL,
                lower_full = NULL,
                upper_full = NULL
            ))
        }

        point <- cov_bounds(
            Sigma_rel_s,
            sigma_L = sigma_L_s,
            sigma_U = sigma_U_s,
            rho_L = rho_L_s,
            rho_U = rho_U_s,
            rho_witness = rho_witness_s,
            solver = solver,
            eig_tol = high_resolution_eig_tol
        )

        if (is.null(point$lower) || is.null(point$upper) ||
            .bad_bounds_vec(as.vector(point$lower), as.vector(point$upper))) {
            return(.draw_fail(s, "bad_cov_bounds"))
        }

        list(
            s = s,
            valid = TRUE,
            skip_reason = NA_character_,
            sigma_L_track = sigma_L_track_s,
            sigma_U_track = sigma_U_track_s,
            rho_L_track = rho_L_track_s,
            rho_U_track = rho_U_track_s,
            rel_off = NULL,
            lower_off = NULL,
            upper_off = NULL,
            rel_diag = NULL,
            lower_diag = NULL,
            upper_diag = NULL,
            rel_full = rel_full_s,
            lower_full = point$lower,
            upper_full = point$upper,
            Sigma_rel = Sigma_rel_s,
            sigma_L = point$inputs$sigma_L,
            sigma_U = point$inputs$sigma_U,
            rho_L = point$inputs$rho_L,
            rho_U = point$inputs$rho_U,
            rho_witness = rho_witness_s
        )
    }

    .is_retriable_draw_error <- function(message) {
        patterns <- c(
            "fixed correlations are incompatible",
            "Empty E intersection C",
            "rho bounds are incompatible",
            "relative-covariance geometry"
        )
        any(vapply(patterns, grepl, logical(1), x = message, fixed = TRUE))
    }

    .process_one_draw <- function(s, stream_mode, high_resolution_pairs = NULL) {
        rejected <- character()
        for (attempt in seq_len(if (bootstrap_active) max_attempts_per_draw else 1L)) {
            result <- tryCatch(
                .process_one_attempt(s, attempt, stream_mode, high_resolution_pairs),
                error = function(error) {
                    message <- conditionMessage(error)
                    if (.is_retriable_draw_error(message)) {
                        return(.draw_fail(s, "incompatible_covariance_parameters"))
                    }
                    stop(error)
                }
            )
            if (isTRUE(result$valid)) {
                result$attempts <- attempt
                result$rejection_reasons <- rejected
                return(result)
            }
            rejected <- c(rejected, result$skip_reason %||% "unknown_draw_failure")
        }
        result$attempts <- length(rejected)
        result$rejection_reasons <- rejected
        result
    }

    .process_stream_chunk <- function(draw_ids) {
        n_local <- length(draw_ids)
        fields <- c(
            "lower_off", "upper_off", "rel_off",
            "lower_diag", "upper_diag", "rel_diag"
        )
        paths <- stats::setNames(
            vapply(
                fields,
                function(field) tempfile(paste0("prism_", field, "_"), fileext = ".bin"),
                character(1)
            ),
            fields
        )
        connections <- lapply(paths, file, open = "wb")
        on.exit(lapply(connections, close), add = TRUE)

        out <- list(
            draw_ids = draw_ids,
            valid = rep(FALSE, n_local),
            skip_reason = rep(NA_character_, n_local),
            sigma_L_track = rep(NA_real_, n_local),
            sigma_U_track = rep(NA_real_, n_local),
            rho_L_track = rep(NA_real_, n_local),
            rho_U_track = rep(NA_real_, n_local),
            attempts = integer(n_local),
            rejection_reasons = vector("list", n_local),
            paths = paths
        )

        write_values <- function(field, values, expected_length) {
            if (is.null(values)) values <- rep(NA_real_, expected_length)
            if (length(values) != expected_length) {
                stop("Internal streamed draw length mismatch for ", field, ".", call. = FALSE)
            }
            writeBin(
                as.double(values),
                connections[[field]],
                size = 8L,
                endian = "little"
            )
        }

        for (ii in seq_along(draw_ids)) {
            res <- .process_one_draw(draw_ids[ii], stream_mode = TRUE)

            out$valid[ii] <- isTRUE(res$valid)
            out$skip_reason[ii] <- res$skip_reason
            out$sigma_L_track[ii] <- res$sigma_L_track
            out$sigma_U_track[ii] <- res$sigma_U_track
            out$rho_L_track[ii] <- res$rho_L_track
            out$rho_U_track[ii] <- res$rho_U_track
            out$attempts[ii] <- res$attempts
            out$rejection_reasons[[ii]] <- res$rejection_reasons

            write_values("lower_off", res$lower_off, K_off)
            write_values("upper_off", res$upper_off, K_off)
            write_values("rel_off", res$rel_off, K_off)
            write_values("lower_diag", res$lower_diag, K_diag)
            write_values("upper_diag", res$upper_diag, K_diag)
            write_values("rel_diag", res$rel_diag, K_diag)
        }

        out
    }

    # -- S-draw processing -------------------------------------------------------
    if (isTRUE(stream_cfg$use_stream)) {
        if (workers > 1L) {
            if (!requireNamespace("future.apply", quietly = TRUE)) {
                stop(
                    "workers > 1 requires the 'future.apply' package. Install with: install.packages('future.apply')",
                    call. = FALSE
                )
            }
            if (!requireNamespace("future", quietly = TRUE)) {
                stop(
                    "workers > 1 requires the 'future' package. Install with: install.packages('future')",
                    call. = FALSE
                )
            }

            chunk_size <- ceiling(S / workers)
            draw_chunks <- split(seq_len(S), ceiling(seq_along(seq_len(S)) / chunk_size))

            old_future_plan <- future::plan()
            on.exit(future::plan(old_future_plan), add = TRUE)
            if (future::supportsMulticore()) {
                future::plan(future::multicore, workers = workers)
            } else {
                future::plan(future::multisession, workers = workers)
            }

            stream_results <- future.apply::future_lapply(
                draw_chunks,
                .process_stream_chunk,
                # Every stochastic operation uses a precomputed draw-specific
                # integer seed. Do not let future.apply switch workers to an
                # L'Ecuyer stream, which would change set.seed(seed) results
                # relative to serial Mersenne-Twister execution.
                future.seed = NULL
            )
        } else {
            stream_results <- list(.process_stream_chunk(seq_len(S)))
        }

        for (res in stream_results) {
            rows <- res$draw_ids

            valid_draw[rows] <- res$valid
            skip_reason[rows] <- res$skip_reason
            sigma_L_track[rows] <- res$sigma_L_track
            sigma_U_track[rows] <- res$sigma_U_track
            rho_L_track[rows] <- res$rho_L_track
            rho_U_track[rows] <- res$rho_U_track
            draw_attempts[rows] <- res$attempts
            draw_rejection_reasons[rows] <- res$rejection_reasons
        }
        stream_temp_files <- unlist(lapply(stream_results, `[[`, "paths"), use.names = FALSE)
        on.exit(unlink(stream_temp_files), add = TRUE)
    } else {
        for (s in seq_len(S)) {
            res <- .process_one_draw(s, stream_mode = FALSE)

            valid_draw[s] <- isTRUE(res$valid)
            skip_reason[s] <- res$skip_reason
            sigma_L_track[s] <- res$sigma_L_track
            sigma_U_track[s] <- res$sigma_U_track
            rho_L_track[s] <- res$rho_L_track
            rho_U_track[s] <- res$rho_U_track
            draw_attempts[s] <- res$attempts
            draw_rejection_reasons[[s]] <- res$rejection_reasons

            if (isTRUE(res$valid)) {
                rel_draws[s, ] <- res$rel_full
                lower_arr[, , s] <- res$lower_full
                upper_arr[, , s] <- res$upper_full
                if (high_resolution_config$enabled) {
                    Sigma_rel_high_resolution[s] <- list(res$Sigma_rel)
                    sigma_L_high_resolution[s] <- res$sigma_L
                    sigma_U_high_resolution[s] <- res$sigma_U
                    rho_L_high_resolution[s] <- list(res$rho_L)
                    rho_U_high_resolution[s] <- list(res$rho_U)
                    rho_witness_high_resolution[s] <- list(res$rho_witness)
                }
            }
        }
    }

    # ---- aggregate draws -------------------------------------------------------
    n_valid <- sum(valid_draw)
    rejection_reasons <- unlist(draw_rejection_reasons, use.names = FALSE)
    n_rejected_attempts <- length(rejection_reasons)
    n_attempts <- n_valid + n_rejected_attempts
    sampling_diagnostics <- list(
        bootstrap_mode = bootstrap_mode,
        requested_draws = S,
        accepted_draws = n_valid,
        attempted_draws = n_attempts,
        rejected_attempts = n_rejected_attempts,
        acceptance_fraction = if (n_attempts > 0L) n_valid / n_attempts else 0,
        rejection_fraction = if (n_attempts > 0L) n_rejected_attempts / n_attempts else 0,
        draws_requiring_replacement = sum(valid_draw & draw_attempts > 1L),
        attempts_per_requested_draw = c(
            min = min(draw_attempts),
            median = stats::median(draw_attempts),
            max = max(draw_attempts)
        ),
        rejection_reasons = sort(table(rejection_reasons), decreasing = TRUE),
        rejection_reason_fractions = if (n_attempts > 0L) {
            sort(table(rejection_reasons), decreasing = TRUE) / n_attempts
        } else {
            numeric()
        },
        max_attempts_per_draw = max_attempts_per_draw
    )

    if (bootstrap_active && n_valid < S) {
        tab <- sort(table(rejection_reasons), decreasing = TRUE)
        stop(
            paste0(
                "prism(): replacement bootstrap exhausted before obtaining S valid draws. ",
                "Valid draws = ", n_valid, " out of ", S, ".\n",
                "Rejected-attempt reasons:\n",
                paste(sprintf("  %s: %d", names(tab), as.integer(tab)), collapse = "\n")
            ),
            call. = FALSE
        )
    }

    if (n_valid < 2L) {
        tab <- sort(table(skip_reason[!is.na(skip_reason)]), decreasing = TRUE)
        stop(
            paste0(
                "prism(): too few valid stochastic draws after filtering. ",
                "Valid draws = ", n_valid, " out of ", S, ".\n",
                "Skip reasons:\n",
                paste(sprintf("  %s: %d", names(tab), as.integer(tab)), collapse = "\n")
            )
        )
    }

    if (verbose && n_rejected_attempts > 0L) {
        tab <- sort(table(rejection_reasons), decreasing = TRUE)
        message(
            paste0(
                "[prism] Replaced ", n_rejected_attempts,
                " rejected stochastic attempt(s) while obtaining ", S,
                " accepted draw(s). Reasons:\n",
                paste(sprintf("  %s: %d", names(tab), as.integer(tab)), collapse = "\n")
            )
        )
    }

    if (isTRUE(stream_cfg$use_stream)) {
        read_stream_block <- function(field, n_cols, start, end) {
            width <- end - start + 1L
            pieces <- lapply(stream_results, function(chunk) {
                n_local <- length(chunk$draw_ids)
                values <- matrix(NA_real_, nrow = n_local, ncol = width)
                con <- file(chunk$paths[[field]], open = "rb")
                on.exit(close(con), add = TRUE)
                for (row in seq_len(n_local)) {
                    offset <- ((row - 1L) * n_cols + (start - 1L)) * 8
                    seek(con, where = offset, origin = "start")
                    value <- readBin(
                        con,
                        what = "double",
                        n = width,
                        size = 8L,
                        endian = "little"
                    )
                    if (length(value) != width) {
                        stop("Unexpected end of a PRISM stream file.", call. = FALSE)
                    }
                    values[row, ] <- value
                }
                values[chunk$valid, , drop = FALSE]
            })
            do.call(rbind, pieces)
        }

        read_stream_columns <- function(field, n_cols, columns) {
            do.call(
                cbind,
                lapply(
                    columns,
                    function(column) read_stream_block(
                        field,
                        n_cols,
                        column,
                        column
                    )
                )
            )
        }

        pair_block_size <- max(
            1L,
            as.integer(floor(stream_cfg$block_size / max(1L, n_valid)))
        )

        aggregate_stream_group <- function(prefix, K, calculate_inference) {
            out <- list(
                ci_lower = numeric(K),
                ci_upper = numeric(K),
                id_width = numeric(K),
                lower_iqr = numeric(K),
                upper_iqr = numeric(K),
                rel_hat = numeric(K),
                rel_ci_lower = numeric(K),
                rel_ci_upper = numeric(K),
                p_value = rep(NA_real_, K),
                rel_draws = if (isTRUE(return_rel_draws)) {
                    matrix(NA_real_, nrow = n_valid, ncol = K)
                } else {
                    NULL
                },
                lower_draws = if (isTRUE(return_bound_draws)) {
                    matrix(NA_real_, nrow = n_valid, ncol = K)
                } else {
                    NULL
                },
                upper_draws = if (isTRUE(return_bound_draws)) {
                    matrix(NA_real_, nrow = n_valid, ncol = K)
                } else {
                    NULL
                }
            )
            if (K == 0L) return(out)

            for (start in seq.int(1L, K, by = pair_block_size)) {
                end <- min(K, start + pair_block_size - 1L)
                cols <- start:end
                lower <- read_stream_block(paste0("lower_", prefix), K, start, end)
                upper <- read_stream_block(paste0("upper_", prefix), K, start, end)
                relative <- read_stream_block(paste0("rel_", prefix), K, start, end)

                out$ci_lower[cols] <- apply(
                    lower,
                    2,
                    stats::quantile,
                    probs = q_lo,
                    names = FALSE,
                    type = 7,
                    na.rm = TRUE
                )
                out$ci_upper[cols] <- apply(
                    upper,
                    2,
                    stats::quantile,
                    probs = q_hi,
                    names = FALSE,
                    type = 7,
                    na.rm = TRUE
                )
                out$id_width[cols] <- apply(
                    upper - lower,
                    2,
                    stats::median,
                    na.rm = TRUE
                )
                out$lower_iqr[cols] <- apply(lower, 2, stats::IQR, na.rm = TRUE)
                out$upper_iqr[cols] <- apply(upper, 2, stats::IQR, na.rm = TRUE)
                out$rel_hat[cols] <- apply(relative, 2, stats::median, na.rm = TRUE)
                out$rel_ci_lower[cols] <- apply(
                    relative,
                    2,
                    stats::quantile,
                    probs = q_lo,
                    names = FALSE,
                    type = 7,
                    na.rm = TRUE
                )
                out$rel_ci_upper[cols] <- apply(
                    relative,
                    2,
                    stats::quantile,
                    probs = q_hi,
                    names = FALSE,
                    type = 7,
                    na.rm = TRUE
                )
                if (isTRUE(return_rel_draws)) out$rel_draws[, cols] <- relative
                if (isTRUE(return_bound_draws)) {
                    out$lower_draws[, cols] <- lower
                    out$upper_draws[, cols] <- upper
                }

                if (isTRUE(calculate_inference)) {
                    ok <- !is.na(lower) & !is.na(upper)
                    valid_count <- colSums(ok)
                    pos_count <- colSums(ok & lower <= delta)
                    neg_count <- colSums(ok & upper >= -delta)
                    has_valid <- valid_count > 0L
                    p_block <- rep(NA_real_, length(cols))
                    denom <- 1L + valid_count[has_valid]
                    p_plus <- (1 + pos_count[has_valid]) / denom
                    p_minus <- (1 + neg_count[has_valid]) / denom
                    p_block[has_valid] <- pmin(1, 2 * pmin(p_plus, p_minus))
                    out$p_value[cols] <- p_block
                }
            }
            out
        }

        off_summary <- aggregate_stream_group("off", K_off, TRUE)
        diag_summary <- aggregate_stream_group("diag", K_diag, FALSE)
        high_resolution_diagnostics <- NULL
        if (high_resolution_config$enabled) {
            refine <- .select_requested_high_resolution_pairs(
                cbind(i_off, j_off),
                high_resolution_pair_index,
                off_summary$ci_lower,
                off_summary$ci_upper,
                delta,
                high_resolution_config$policy,
                high_resolution_config$eps
            )
            if (any(refine)) {
                selected <- which(refine)
                standard_selected_lower <- off_summary$ci_lower[selected]
                standard_selected_upper <- off_summary$ci_upper[selected]
                valid_ids <- which(valid_draw)
                any_sharp_draw <- FALSE

                for (block_start in seq.int(1L, length(selected), by = pair_block_size)) {
                    block_end <- min(
                        length(selected),
                        block_start + pair_block_size - 1L
                    )
                    selected_block <- selected[block_start:block_end]
                    pair_block <- cbind(i_off[selected_block], j_off[selected_block])
                    standard_lower_draw <- read_stream_columns(
                        "lower_off",
                        K_off,
                        selected_block
                    )
                    standard_upper_draw <- read_stream_columns(
                        "upper_off",
                        K_off,
                        selected_block
                    )
                    sharp_lower_draw <- standard_lower_draw
                    sharp_upper_draw <- standard_upper_draw

                    for (row in seq_along(valid_ids)) {
                        sharp <- .process_one_draw(
                            valid_ids[row],
                            stream_mode = TRUE,
                            high_resolution_pairs = pair_block
                        )
                        if (!isTRUE(sharp$valid)) {
                            stop(
                                "High-resolution streamed recomputation failed for draw ",
                                valid_ids[row],
                                ": ",
                                sharp$skip_reason,
                                call. = FALSE
                            )
                        }
                        if (!isTRUE(sharp$sharp_applied)) next
                        tightened <- .validate_high_resolution_tightening(
                            standard_lower_draw[row, ],
                            standard_upper_draw[row, ],
                            sharp$lower,
                            sharp$upper
                        )
                        sharp_lower_draw[row, ] <- tightened$lower
                        sharp_upper_draw[row, ] <- tightened$upper
                        any_sharp_draw <- TRUE
                    }

                    off_summary$ci_lower[selected_block] <- apply(
                        sharp_lower_draw,
                        2,
                        stats::quantile,
                        probs = q_lo,
                        names = FALSE,
                        type = 7,
                        na.rm = TRUE
                    )
                    off_summary$ci_upper[selected_block] <- apply(
                        sharp_upper_draw,
                        2,
                        stats::quantile,
                        probs = q_hi,
                        names = FALSE,
                        type = 7,
                        na.rm = TRUE
                    )
                    off_summary$id_width[selected_block] <- apply(
                        sharp_upper_draw - sharp_lower_draw,
                        2,
                        stats::median,
                        na.rm = TRUE
                    )
                    off_summary$lower_iqr[selected_block] <- apply(
                        sharp_lower_draw,
                        2,
                        stats::IQR,
                        na.rm = TRUE
                    )
                    off_summary$upper_iqr[selected_block] <- apply(
                        sharp_upper_draw,
                        2,
                        stats::IQR,
                        na.rm = TRUE
                    )
                    if (isTRUE(return_bound_draws)) {
                        off_summary$lower_draws[, selected_block] <- sharp_lower_draw
                        off_summary$upper_draws[, selected_block] <- sharp_upper_draw
                    }

                    ok <- !is.na(sharp_lower_draw) & !is.na(sharp_upper_draw)
                    valid_count <- colSums(ok)
                    pos_count <- colSums(ok & sharp_lower_draw <= delta)
                    neg_count <- colSums(ok & sharp_upper_draw >= -delta)
                    has_valid <- valid_count > 0L
                    p_block <- rep(NA_real_, length(selected_block))
                    denom <- 1L + valid_count[has_valid]
                    p_plus <- (1 + pos_count[has_valid]) / denom
                    p_minus <- (1 + neg_count[has_valid]) / denom
                    p_block[has_valid] <- pmin(1, 2 * pmin(p_plus, p_minus))
                    off_summary$p_value[selected_block] <- p_block
                }

                if (any_sharp_draw) {
                    high_resolution_diagnostics <- .high_resolution_ci_diagnostics(
                        i_off[selected],
                        j_off[selected],
                        standard_selected_lower,
                        standard_selected_upper,
                        off_summary$ci_lower[selected],
                        off_summary$ci_upper[selected]
                    )
                } else {
                    warning(
                        "high_resolution=TRUE is unnecessary because no valid draw has nonzero-width rho intervals.",
                        call. = FALSE
                    )
                }
            }
        }
        ci_lower_off <- off_summary$ci_lower
        ci_upper_off <- off_summary$ci_upper
        p_value <- off_summary$p_value
        q_value <- rep(NA_real_, K_off)
        has_valid <- is.finite(p_value)
        q_value[has_valid] <- stats::p.adjust(p_value[has_valid], method = p_adjust_method)

        pairwise <- data.frame(
            i = i_off,
            j = j_off,
            ci_lower = ci_lower_off,
            ci_upper = ci_upper_off,
            p_value = p_value,
            q_value = q_value,
            delta = delta
        )

        ci_lower_diag <- diag_summary$ci_lower
        ci_upper_diag <- diag_summary$ci_upper

        ci_lower <- matrix(NA_real_, D, D)
        ci_upper <- matrix(NA_real_, D, D)
        id_region_width_median <- matrix(NA_real_, D, D)
        lower_iqr <- matrix(NA_real_, D, D)
        upper_iqr <- matrix(NA_real_, D, D)

        ci_lower[cbind(i_off, j_off)] <- ci_lower_off
        ci_lower[cbind(j_off, i_off)] <- ci_lower_off
        ci_upper[cbind(i_off, j_off)] <- ci_upper_off
        ci_upper[cbind(j_off, i_off)] <- ci_upper_off

        id_width_off <- off_summary$id_width
        lower_iqr_off <- off_summary$lower_iqr
        upper_iqr_off <- off_summary$upper_iqr
        id_region_width_median[cbind(i_off, j_off)] <- id_width_off
        id_region_width_median[cbind(j_off, i_off)] <- id_width_off
        lower_iqr[cbind(i_off, j_off)] <- lower_iqr_off
        lower_iqr[cbind(j_off, i_off)] <- lower_iqr_off
        upper_iqr[cbind(i_off, j_off)] <- upper_iqr_off
        upper_iqr[cbind(j_off, i_off)] <- upper_iqr_off

        if (K_diag > 0L) {
            diag(ci_lower) <- ci_lower_diag
            diag(ci_upper) <- ci_upper_diag
            diag(id_region_width_median) <- diag_summary$id_width
            diag(lower_iqr) <- diag_summary$lower_iqr
            diag(upper_iqr) <- diag_summary$upper_iqr
        }

        rel_i <- c(i_diag, i_off)
        rel_j <- c(j_diag, j_off)
        rel_order <- match(
            paste(i_rel, j_rel, sep = "_"),
            paste(rel_i, rel_j, sep = "_")
        )
        rel_hat <- c(diag_summary$rel_hat, off_summary$rel_hat)[rel_order]
        rel_ci_lower <- c(
            diag_summary$rel_ci_lower,
            off_summary$rel_ci_lower
        )[rel_order]
        rel_ci_upper <- c(
            diag_summary$rel_ci_upper,
            off_summary$rel_ci_upper
        )[rel_order]
        rel_combined <- if (isTRUE(return_rel_draws)) {
            cbind(diag_summary$rel_draws, off_summary$rel_draws)[, rel_order, drop = FALSE]
        } else {
            NULL
        }
        pairwise_rel <- data.frame(
            i = i_rel,
            j = j_rel,
            rel_hat = rel_hat,
            rel_ci_lower = rel_ci_lower,
            rel_ci_upper = rel_ci_upper,
            rel_ci_width = rel_ci_upper - rel_ci_lower
        )
        pairwise_bound_draws <- make_pairwise_bound_draws(
            draw_id = which(valid_draw),
            i = i_off,
            j = j_off,
            lower = off_summary$lower_draws,
            upper = off_summary$upper_draws
        )

        scale_summary <- if (has_scale_log || has_scale) {
            valid_s <- is.finite(sigma_L_track) & is.finite(sigma_U_track)
            list(
                sigma_L_med = stats::median(sigma_L_track[valid_s]),
                sigma_U_med = stats::median(sigma_U_track[valid_s]),
                rho_L_med = stats::median(rho_L_track[valid_s], na.rm = TRUE),
                rho_U_med = stats::median(rho_U_track[valid_s], na.rm = TRUE)
            )
        } else {
            NULL
        }

        return(.finalize_prism_result(list(
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            point_lower = NULL,
            point_upper = NULL,
            pairwise = pairwise,
            pairwise_rel = pairwise_rel,
            high_resolution_diagnostics = high_resolution_diagnostics,
            rel_draws = if (isTRUE(return_rel_draws)) rel_combined else NULL,
            pairwise_bound_draws = pairwise_bound_draws,
            mln_fit_summary = mln_fit_summary,
            mln_diagnostics = mln_diagnostics,
            mln_fit_object = if (composition == "MLN" && isTRUE(mln_options$return_fit)) mln_fit_object else NULL,
            data_diagnostics = data_diagnostics,
            sampling_diagnostics = sampling_diagnostics,
            id_region_width_median = id_region_width_median,
            lower_iqr = lower_iqr,
            upper_iqr = upper_iqr,
            params = list(
                composition = composition,
                bootstrap = bootstrap_active,
                bootstrap_mode = bootstrap_mode,
                S = S,
                alpha = alpha,
                dirichlet_concentration = dirichlet_concentration,
                composition_storage = composition_storage,
                delta = delta,
                stream = stream,
                workers = workers,
                return_bound_draws = return_bound_draws,
                stream_active = TRUE,
                stream_mode = stream_cfg$mode,
                stream_block_size = stream_cfg$block_size,
                stream_pair_block_size = pair_block_size,
                stream_storage = "temporary_binary_files",
                stream_n_offdiag_pairs = K_off,
                stream_n_diag_pairs = K_diag,
                stream_diagonal_policy = "always_include_i_equals_i_in_pairwise_rel",
                scale_summary = scale_summary,
                high_resolution = high_resolution_config$enabled,
                high_resolution_policy = high_resolution_config$policy,
                high_resolution_eps = high_resolution_config$eps,
                high_resolution_selection = if (is.null(high_resolution_pair_index)) "policy" else "requested_pairs",
                high_resolution_requested_pairs = if (is.null(high_resolution_pair_index)) 0L else nrow(high_resolution_pair_index),
                solver = solver,
                scale_uncertainty_draws = scale_uncertainty_draws
            )
        ), composition_estimator))
    }

    lower_arr_valid <- lower_arr[, , valid_draw, drop = FALSE]
    upper_arr_valid <- upper_arr[, , valid_draw, drop = FALSE]
    rel_draws_valid <- rel_draws[valid_draw, , drop = FALSE]

    ci_lower <- apply(lower_arr_valid, c(1, 2), stats::quantile, probs = q_lo, names = FALSE, type = 7, na.rm = TRUE)
    ci_upper <- apply(upper_arr_valid, c(1, 2), stats::quantile, probs = q_hi, names = FALSE, type = 7, na.rm = TRUE)

    idx_pairs <- which(upper.tri(ci_lower), arr.ind = TRUE)
    i_vec <- idx_pairs[, 1]
    j_vec <- idx_pairs[, 2]
    n_pairs <- nrow(idx_pairs)

    high_resolution_diagnostics <- NULL
    if (high_resolution_config$enabled) {
        valid_indices <- which(valid_draw)
        interval_draw <- vapply(valid_indices, function(s) {
            rL <- rho_L_high_resolution[[s]]
            rU <- rho_U_high_resolution[[s]]
            !is.null(rL) && !is.null(rU) && any(rU > rL)
        }, logical(1))
        if (!any(interval_draw)) {
            warning("high_resolution=TRUE is unnecessary because no valid draw has nonzero-width rho intervals.", call. = FALSE)
        } else {
            standard_pair_lower <- ci_lower[idx_pairs]
            standard_pair_upper <- ci_upper[idx_pairs]
            refine <- .select_requested_high_resolution_pairs(
                idx_pairs,
                high_resolution_pair_index,
                standard_pair_lower,
                standard_pair_upper,
                delta,
                high_resolution_config$policy,
                high_resolution_config$eps
            )
            if (any(refine)) {
                selected_pairs <- idx_pairs[refine, , drop = FALSE]
                for (s in valid_indices) {
                    rL <- rho_L_high_resolution[[s]]
                    rU <- rho_U_high_resolution[[s]]
                    if (is.null(rL) || is.null(rU) || all(rL == rU)) next
                    standard_lower_s <- lower_arr[cbind(selected_pairs[, 1], selected_pairs[, 2], s)]
                    standard_upper_s <- upper_arr[cbind(selected_pairs[, 1], selected_pairs[, 2], s)]
                    sharp_s <- cov_bounds_pairs(
                        Sigma_rel_high_resolution[[s]],
                        pair_index = selected_pairs,
                        sigma_L = sigma_L_high_resolution[s],
                        sigma_U = sigma_U_high_resolution[s],
                        rho_L = rL,
                        rho_U = rU,
                        rho_witness = rho_witness_high_resolution[[s]],
                        high_resolution = TRUE,
                        solver = solver,
                        eig_tol = high_resolution_eig_tol
                    )
                    tightened_s <- .validate_high_resolution_tightening(
                        standard_lower_s,
                        standard_upper_s,
                        sharp_s$lower,
                        sharp_s$upper
                    )
                    for (k in seq_len(nrow(selected_pairs))) {
                        ii <- selected_pairs[k, 1]
                        jj <- selected_pairs[k, 2]
                        lower_arr[ii, jj, s] <- lower_arr[jj, ii, s] <- tightened_s$lower[k]
                        upper_arr[ii, jj, s] <- upper_arr[jj, ii, s] <- tightened_s$upper[k]
                    }
                }
                lower_arr_valid <- lower_arr[, , valid_draw, drop = FALSE]
                upper_arr_valid <- upper_arr[, , valid_draw, drop = FALSE]
                ci_lower <- apply(lower_arr_valid, c(1, 2), stats::quantile, probs = q_lo, names = FALSE, type = 7, na.rm = TRUE)
                ci_upper <- apply(upper_arr_valid, c(1, 2), stats::quantile, probs = q_hi, names = FALSE, type = 7, na.rm = TRUE)
                high_resolution_diagnostics <- .high_resolution_ci_diagnostics(
                    idx_pairs[refine, 1],
                    idx_pairs[refine, 2],
                    standard_pair_lower[refine],
                    standard_pair_upper[refine],
                    ci_lower[idx_pairs][refine],
                    ci_upper[idx_pairs][refine]
                )
            }
        }
    }

    id_width_arr <- upper_arr_valid - lower_arr_valid
    id_region_width_median <- apply(id_width_arr, c(1, 2), stats::median, na.rm = TRUE)
    lower_iqr <- apply(lower_arr, c(1, 2), stats::IQR, na.rm = TRUE)
    upper_iqr <- apply(upper_arr, c(1, 2), stats::IQR, na.rm = TRUE)

    S_eff <- n_valid
    lower_flat <- matrix(aperm(lower_arr_valid, c(3, 1, 2)), nrow = S_eff, ncol = D * D)
    upper_flat <- matrix(aperm(upper_arr_valid, c(3, 1, 2)), nrow = S_eff, ncol = D * D)
    col_idx <- (j_vec - 1L) * D + i_vec

    chunk_pairs <- 5000L
    pos_count <- integer(n_pairs)
    neg_count <- integer(n_pairs)
    valid_count <- integer(n_pairs)

    for (start in seq.int(1L, n_pairs, by = chunk_pairs)) {
        end <- min(n_pairs, start + chunk_pairs - 1L)
        cols <- col_idx[start:end]
        Ls <- lower_flat[, cols, drop = FALSE]
        Us <- upper_flat[, cols, drop = FALSE]
        ok <- !is.na(Ls) & !is.na(Us)
        valid_count[start:end] <- colSums(ok)
        pos_count[start:end] <- colSums(ok & (Ls <= delta))
        neg_count[start:end] <- colSums(ok & (Us >= -delta))
    }

    p_value <- rep(NA_real_, n_pairs)
    has_valid <- valid_count > 0
    denom <- 1L + valid_count[has_valid]
    p_plus <- (1 + pos_count[has_valid]) / denom
    p_minus <- (1 + neg_count[has_valid]) / denom
    p_value[has_valid] <- pmin(1, 2 * pmin(p_plus, p_minus))
    q_value <- rep(NA_real_, n_pairs)
    q_value[has_valid] <- stats::p.adjust(p_value[has_valid], method = p_adjust_method)

    pairwise <- data.frame(
        i = i_vec, j = j_vec,
        ci_lower = ci_lower[idx_pairs],
        ci_upper = ci_upper[idx_pairs],
        p_value = p_value,
        q_value = q_value,
        delta = delta
    )

    rel_hat <- apply(rel_draws_valid, 2, stats::median, na.rm = TRUE)
    rel_ci_lower <- apply(rel_draws_valid, 2, stats::quantile, probs = q_lo, names = FALSE, type = 7, na.rm = TRUE)
    rel_ci_upper <- apply(rel_draws_valid, 2, stats::quantile, probs = q_hi, names = FALSE, type = 7, na.rm = TRUE)

    pairwise_rel <- data.frame(
        i = i_rel, j = j_rel,
        rel_hat = rel_hat,
        rel_ci_lower = rel_ci_lower,
        rel_ci_upper = rel_ci_upper,
        rel_ci_width = rel_ci_upper - rel_ci_lower
    )
    pairwise_bound_draws <- make_pairwise_bound_draws(
        draw_id = which(valid_draw),
        i = i_vec,
        j = j_vec,
        lower = lower_flat[, col_idx, drop = FALSE],
        upper = upper_flat[, col_idx, drop = FALSE]
    )

    if (!isTRUE(return_rel_draws)) rel_draws <- NULL

    # build scale_summary from tracking vectors if scale_log was used
    scale_summary <- if (has_scale_log || has_scale) {
        valid_s <- is.finite(sigma_L_track) & is.finite(sigma_U_track)
        list(
            sigma_L_med = stats::median(sigma_L_track[valid_s]),
            sigma_U_med = stats::median(sigma_U_track[valid_s]),
            rho_L_med   = stats::median(rho_L_track[valid_s], na.rm = TRUE),
            rho_U_med   = stats::median(rho_U_track[valid_s], na.rm = TRUE)
        )
    } else {
        NULL
    }

    .finalize_prism_result(list(
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        point_lower = NULL,
        point_upper = NULL,
        pairwise = pairwise,
        pairwise_rel = pairwise_rel,
        high_resolution_diagnostics = high_resolution_diagnostics,
        rel_draws = rel_draws,
        pairwise_bound_draws = pairwise_bound_draws,
        mln_fit_summary = mln_fit_summary,
        mln_diagnostics = mln_diagnostics,
        mln_fit_object = if (composition == "MLN" && isTRUE(mln_options$return_fit)) mln_fit_object else NULL,
        data_diagnostics = data_diagnostics,
        sampling_diagnostics = sampling_diagnostics,
        id_region_width_median = id_region_width_median,
        lower_iqr = lower_iqr,
        upper_iqr = upper_iqr,
        params = list(
            composition = composition,
            bootstrap = bootstrap_active,
            bootstrap_mode = bootstrap_mode,
            S = S,
            alpha = alpha,
            dirichlet_concentration = dirichlet_concentration,
            composition_storage = composition_storage,
            delta = delta,
            stream = stream,
            workers = workers,
            return_bound_draws = return_bound_draws,
            stream_active = FALSE,
            stream_mode = stream_cfg$mode,
            stream_block_size = stream_cfg$block_size,
            scale_summary = scale_summary,
            high_resolution = high_resolution_config$enabled,
            high_resolution_policy = high_resolution_config$policy,
            high_resolution_eps = high_resolution_config$eps,
            high_resolution_selection = if (is.null(high_resolution_pair_index)) "policy" else "requested_pairs",
            high_resolution_requested_pairs = if (is.null(high_resolution_pair_index)) 0L else nrow(high_resolution_pair_index),
            solver = solver,
            scale_uncertainty_draws = scale_uncertainty_draws
        )
    ), composition_estimator)
}
