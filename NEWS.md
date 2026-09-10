# prism 0.2.0

Ground-up rewrite. The statistical target is unchanged (identification
bounds and bootstrap confidence intervals for absolute log-covariance of
compositional data with optional scale information), but the
implementation is new and several features are gone rather than
preserved:

- Replaced the 2213-line single-function engine and the 1215-line bounds
  file with nine small, single-purpose files (`relative-covariance.R`,
  `composition.R`, `scale.R`, `bounds.R`, `bootstrap.R`, `aggregate.R`,
  `prism.R`, `result.R`, `utils.R`), each independently unit-tested.
- `cov_bounds()` is now one function (was `cov_bounds()` +
  `cov_bounds_pairs()` + a separate diagonal-only helper, three
  implementations of the same formula). Verified algebraically and by
  test that a single vertex-of-quadratic-in-sigma kernel reproduces the
  old unbounded/bounded-scale/diagonal formulas as special cases.
- Removed the CVXR-based "high resolution" sharp-bounds feature entirely
  (`high_resolution*` arguments, `plot_high_resolution_*()`). The
  rho-box feasibility check it depended on is kept, since it is a
  correctness gate, not part of the sharpening.
- Removed disk-backed streaming and `future.apply` parallel workers
  (`stream`, `workers` arguments gone). Draws are held in memory.
- Removed the MLN (`fido`-based) composition estimator: it was largely
  untested, had dead code and an unfinished function stub. The
  `prism_composition_estimator()` plug-in interface remains as the
  extension point for reintroducing it later.
- Removed all PRISM 0.0.x/0.1.x legacy argument compatibility (`...`,
  `alpha`, `samplingmethod`, `counts_are_truth`, `mln_*`,
  `plot_diagnostics`, `composition = "Dirichlet"/"none"/"MLN"` strings,
  `counts_input = "proportions"`).
- Raw scale-replicate uncertainty (`scale` as a replicate matrix) now
  resamples one replicate column per sample per bootstrap draw and
  reuses the same `scale_log` chi-square/Fisher-z CI machinery, instead
  of a separate bespoke inner Monte Carlo "joint witness" search.
- `scale` and `scale_log` are strictly mutually exclusive with each
  other and with `sigma_L`/`sigma_U`/`rho_L`/`rho_U` (an error now,
  rather than one silently overriding another with a warning).
- Test suite rewritten module-by-module; 98.3% line coverage (up from
  78.8%), `R CMD check` clean.
- Existing callers in `SRCoV/` and `stepbystep/code/` in the parent
  repo use the old API and will need updating separately.

# prism 0.1.0

- Rebuilt PRISM as a focused user-facing R package.
- Made Dirichlet-multinomial composition uncertainty the default.
- Added fixed, optional MLN, and custom composition-estimator specifications.
- Isolated MLN fitting and tuning from the PRISM engine.
- Added a documented `prism_result` object and print method.
- Centralized returned diagnostics and removed estimator plotting side effects.
- Retained only the two high-resolution diagnostic plot functions.
- Added bounded-memory disk streaming with full-output equivalence tests.
- Added draw-on-demand Dirichlet generation.
- Removed quadratic dependence on the outer draw count from raw-scale
  uncertainty propagation.
- Added draw-level high-resolution refinement and full-family BH tests.
- Added validated selected-pair high-resolution replay without changing the
  default policy-based interface.
- Added global ellipsoid-box feasibility checks for nonzero-width rho bounds.
- Added same-draw rho witnesses for sample-level scale data and solver-free
  rank-deficient certification when a witness is valid.
- Combined streamed diagonal and off-diagonal bounds into one feasibility check
  per draw.
- Added package, interface, scale-alignment, edge-case, streaming, and load tests.
- Restricted finite-sampling bootstrap to coherent paired sample resampling.
  The former component-only `"composition"` and `"scale"` modes now fail with
  a migration message because their hybrid covariance parameters need not be
  jointly positive semidefinite.
- Used fork-based streamed workers when supported, with socket-based workers as
  the cross-platform fallback.
- Moved analysis code, benchmarks, unrelated methods, data, and binary artifacts
  outside the package source without deleting them.
