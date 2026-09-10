#' Estimate absolute log-covariance bounds from compositional data
#'
#' @description
#' PRISM combines a sample-level compositional uncertainty estimator with exact
#' or bounded information about log scale. The default composition estimator is
#' Dirichlet-multinomial. The method returns identification bounds, Monte Carlo
#' confidence intervals, pairwise inference, and diagnostics without writing
#' files.
#'
#' Counts are features by samples. Scale information may be supplied as
#' sample-level raw measurements (scale), sample-level log measurements
#' (scale_log), or explicit sigma and rho bounds. Correlation bounds require
#' scale-SD bounds.
#'
#' @param counts Non-negative features-by-samples count matrix.
#' @param counts_input Either "counts" or "proportions". Proportions are
#'   accepted only with [prism_composition_fixed()] and pseudocount zero.
#' @param composition A prism composition estimator. Legacy strings
#'   "Dirichlet", "none", and "MLN" are also accepted.
#' @param bootstrap A logical value or one of `"none"` or `"both"`. `FALSE` is
#'   equivalent to `"none"`; `TRUE` is equivalent to `"both"`. Paired
#'   finite-sample resampling applies one sample index to the count composition
#'   and matched scale measurements, then re-estimates every empirical
#'   covariance parameter from that draw. The former `"composition"` and
#'   `"scale"` hybrid modes are deprecated and unsupported because independently
#'   resampled parameter blocks need not define a positive-semidefinite joint
#'   covariance. Sampling counts, fractions, reasons, and attempts-per-draw
#'   summaries are returned under `diagnostics$sampling`.
#' @param S Number of Monte Carlo draws.
#' @param max_attempts_per_draw Maximum bootstrap proposals allowed for each
#'   requested accepted draw. The default is 100. Increase only when explicit
#'   draw-rejection diagnostics show that valid proposals are rare.
#' @param prevalence Minimum feature prevalence retained, in [0, 1].
#' @param verbose Emit progress messages.
#' @param scale Optional positive sample-level scale vector or replicate matrix.
#' @param scale_log Optional finite sample-level log-scale vector.
#' @param scale_log_ci_level Coverage for chi-square scale-SD and Fisher-z
#'   correlation intervals. Zero uses point estimates within each draw.
#' @param scale_log_rho Estimate taxon--scale correlations from `scale_log`.
#' @param scale_log_lower_zero When TRUE, forces the lower end of the
#'   per-draw scale-SD interval to 0 instead of the chi-square lower
#'   quantile, leaving the upper end (and everything else) identical to the
#'   ordinary two-sided regime at the same `scale_log_ci_level`. Sigma
#'   cannot be negative, so the two-sided lower quantile asserts a positive
#'   floor that is an artifact of splitting alpha evenly across both tails,
#'   not a real constraint; this removes that assertion.
#'   Set to `FALSE` to retain paired scale-SD resampling while omitting rho.
#' @param ci_alpha Tail probability for final confidence intervals.
#' @param sigma_L,sigma_U Optional bounds satisfying 0 <= sigma_L <= sigma_U.
#' @param rho_L,rho_U Optional scalar or feature-level correlation bounds in
#'   [-1, 1]; these require sigma bounds. Genuine intervals must contain at
#'   least one joint correlation vector compatible with `Sigma_rel`.
#' @param delta Half-width of the pairwise null region.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param seed Reproducibility seed, or NULL.
#' @param scale_sd_fallback Within-sample log-scale SD used when no sample has
#'   replicate scale measurements.
#' @param scale_uncertainty_draws Inner Monte Carlo draws used to propagate raw
#'   scale replicate uncertainty. This budget is independent of S.
#' @param return_rel_draws Return the relative-covariance draw matrix. This can
#'   be large and necessarily defeats bounded-memory output for that component.
#' @param return_bound_draws Return accepted draw-level lower and upper
#'   off-diagonal covariance bounds under `pairwise_bound_draws`. The two
#'   matrices have draws in rows and the same pair order as `pairwise` in
#'   columns. This can be large and necessarily defeats bounded-memory output
#'   for those components.
#' @param stream NULL for automatic disk-backed streaming, a logical value,
#'   or a positive pair-draw cell budget.
#' @param workers Number of workers for streamed draw processing.
#' @param high_resolution Use sharp numerical bounds for eligible bounded-rho
#'   pairs.
#' @param high_resolution_policy One of "borderline", "all", or "never".
#' @param high_resolution_eps Borderline distance as a fraction of conservative
#'   interval width.
#' @param high_resolution_pair_index Optional two-column matrix or data frame of
#'   off-diagonal integer feature indices to refine. Indices refer to the
#'   retained feature order after PRISM filtering. When supplied, these pairs
#'   replace policy-based pair selection, while refinement still occurs inside
#'   every stochastic draw and multiplicity adjustment still uses all pairs.
#' @param solver CVXR solver used for high-resolution support problems and for
#'   rank-deficient rho-box feasibility checks when no empirical witness or
#'   zero-containing box resolves feasibility.
#' @param high_resolution_eig_tol Relative eigenvalue tolerance for numerical
#'   range-space calculations.
#' @param ... Deprecated PRISM 0.0.x arguments accepted for analysis-script
#'   compatibility. New code should configure composition uncertainty with a
#'   composition-estimator constructor.
#'
#' @return A prism_result containing covariance matrices, pairwise tables,
#'   diagnostics, and effective parameters.
#' @export

prism <- function(
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
  high_resolution_eig_tol = 1e-10,
  ...
) {
  legacy <- list(...)
  supported_legacy <- c(
    "alpha", "dirichlet_concentration", "counts_are_truth", "samplingmethod",
    "return_mln_fit", "return_mln_grid_pairwise", "return_mln_grid_fits",
    "plot_diagnostics", "plot_diagnostics_dir", "plot_diagnostics_pairs",
    "plot_diagnostics_n_pairs", "mln_design", "mln_tune_S",
    "mln_center_scale_grid", "mln_m_grid", "mln_prior_mode",
    "mln_fixed_Gamma", "mln_fixed_Xi", "mln_fixed_upsilon",
    "mln_prior_label", "mln_omega_structure_grid", "mln_gamma_grid",
    "mln_gamma_tune", "mln_xi_grid", "mln_cv_folds", "mu_prior_mean",
    "mu_prior_sd", "sigma_prior_scale", "lkj_eta", "mln_ref",
    "mln_chains", "mln_iter", "mln_warmup", "max_treedepth",
    "adapt_delta", "mln_max_hessian_gb"
  )
  unknown <- setdiff(names(legacy), supported_legacy)
  if (length(unknown)) {
    stop("Unused argument(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.null(legacy$samplingmethod)) composition <- legacy$samplingmethod
  if (isTRUE(legacy$counts_are_truth)) {
    composition <- "none"
    bootstrap <- FALSE
  }
  if (isTRUE(legacy$plot_diagnostics)) {
    warning(
      "plot_diagnostics is deprecated and no files were written; use the two exported high-resolution plot functions.",
      call. = FALSE
    )
  }

  alpha <- legacy$alpha
  concentration <- legacy$dirichlet_concentration
  mln_map <- c(
    mln_design = "design",
    mln_tune_S = "tune_draws",
    mln_center_scale_grid = "center_scale_grid",
    mln_m_grid = "m_grid",
    mln_prior_mode = "prior_mode",
    mln_fixed_Gamma = "fixed_Gamma",
    mln_fixed_Xi = "fixed_Xi",
    mln_fixed_upsilon = "fixed_upsilon",
    mln_prior_label = "prior_label",
    mln_omega_structure_grid = "omega_structure_grid",
    mln_gamma_grid = "gamma_grid",
    mln_max_hessian_gb = "max_hessian_gb",
    return_mln_grid_pairwise = "return_grid_pairwise",
    return_mln_grid_fits = "return_grid_fits",
    return_mln_fit = "return_fit"
  )
  mln_options <- list()
  for (old_name in names(mln_map)) {
    is_return_flag <- old_name %in% c(
      "return_mln_fit", "return_mln_grid_pairwise", "return_mln_grid_fits"
    )
    if (!is.null(legacy[[old_name]]) && (!is_return_flag || isTRUE(legacy[[old_name]]))) {
      mln_options[[mln_map[[old_name]]]] <- legacy[[old_name]]
    }
  }

  if (inherits(composition, "prism_composition_estimator")) {
    if (!is.null(alpha)) composition$options$pseudocount <- .validate_pseudocount(alpha)
    if (!is.null(concentration)) {
      if (composition$type != "dirichlet_multinomial") {
        stop("dirichlet_concentration applies only to the Dirichlet-multinomial estimator.", call. = FALSE)
      }
      composition$options$concentration <- .validate_positive_scalar(
        concentration,
        "dirichlet_concentration"
      )
    }
    if (length(mln_options)) {
      if (composition$type != "mln") {
        stop("MLN compatibility arguments require an MLN composition estimator.", call. = FALSE)
      }
      composition$options[names(mln_options)] <- mln_options
    }
  } else {
    choice <- match.arg(composition, c("Dirichlet", "none", "MLN"))
    if (choice == "Dirichlet" && (!is.null(alpha) || !is.null(concentration))) {
      composition <- prism_composition_dirichlet(
        pseudocount = alpha %||% 1,
        concentration = concentration %||% 1
      )
    } else if (choice == "none" && !is.null(alpha)) {
      composition <- prism_composition_fixed(alpha)
    } else if (choice == "MLN" && (!is.null(alpha) || length(mln_options))) {
      composition <- do.call(
        prism_composition_mln,
        c(list(pseudocount = alpha %||% 1), mln_options)
      )
    }
  }

  .prism_engine(
    counts = counts,
    counts_input = counts_input,
    composition = composition,
    bootstrap = bootstrap,
    S = S,
    max_attempts_per_draw = max_attempts_per_draw,
    prevalence = prevalence,
    verbose = verbose,
    scale = scale,
    scale_log = scale_log,
    scale_log_ci_level = scale_log_ci_level,
    scale_log_rho = scale_log_rho,
    scale_log_lower_zero = scale_log_lower_zero,
    ci_alpha = ci_alpha,
    sigma_L = sigma_L,
    sigma_U = sigma_U,
    rho_L = rho_L,
    rho_U = rho_U,
    delta = delta,
    p_adjust_method = p_adjust_method,
    seed = seed,
    scale_sd_fallback = scale_sd_fallback,
    scale_uncertainty_draws = scale_uncertainty_draws,
    return_rel_draws = return_rel_draws,
    return_bound_draws = return_bound_draws,
    stream = stream,
    workers = workers,
    high_resolution = high_resolution,
    high_resolution_policy = high_resolution_policy,
    high_resolution_eps = high_resolution_eps,
    high_resolution_pair_index = high_resolution_pair_index,
    solver = solver,
    high_resolution_eig_tol = high_resolution_eig_tol
  )
}
