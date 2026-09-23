# prism 0.5.0

- Exact rho-support problems now go directly to ECOS instead of rebuilding
  the same formulation through CVXR for every covariance entry. The original
  CVXR formulation is retained as a regression and benchmark reference.
- Bootstrap draws can now run concurrently with `workers > 1`. Per-draw seeds
  make results invariant to worker count, and results are collected in batches
  no larger than `workers` so only the parent process writes bounded-memory
  streams in draw order.
- Genuine non-fixed `rho_L`/`rho_U` intervals now produce sharp bounds by
  solving the exact ellipsoid-box support problems for every covariance entry
  with conic optimization, including the singular `Sigma_rel` range
  constraint. The previous conservative intersection of separate ellipsoid
  and box support bounds was removed.
- `sigma_U = 0` now returns `Sigma_rel` directly, because all scale terms
  vanish and rho is irrelevant. CVXR result handling supports both the legacy
  list-returning API and the current numeric-returning `psolve()` API.
- Removed bootstrap rejection sampling and replacement draws. Every requested
  draw is now generated exactly once, and invalid or infeasible draws stop with
  their original error. Paired scale estimates pass their known feasible
  `rho_witness` directly to `cov_bounds()`.
- P-values now use opposing tails of the bootstrapped per-draw identification
  intervals and are `NA` when paired subject bootstrapping is disabled.
  `pairwise` also reports `covers_zero` and explicit bootstrap counts and
  frequencies for intervals above zero, below zero, or covering zero.
- Removed the unsupported conversion of an identification bound into an
  implied symmetric Wald p-value.
- Sample-level scale data are now always assumed to be already log transformed;
  bare numeric vectors and matrices are accepted and PRISM never applies a
  logarithm. `prism_scale_log()` remains available for non-default estimation
  settings. Within-draw analytic scale intervals require
  `experimental_ci = TRUE`.
- Constant scale measurements now imply `sigma_L = sigma_U = 0`; rho is omitted
  because it is undefined and irrelevant when sigma is zero. Zero or
  numerically zero marginal variance in a retained log-composition feature
  still raises an error.
- Technical scale replicates are bootstrapped within subject and averaged on
  the log scale rather than selecting a single replicate. Invalid samples are
  removed before prevalence is calculated, and all exclusions warn.
- Removed the bootstrap mean `rel_hat` from `pairwise_rel`; relative-covariance
  uncertainty is reported through its bootstrap interval only.
- Changed the default Dirichlet pseudocount from 1 to 0.5.

# prism 0.4.0

- `prism()` no longer leaks its internal `set.seed(seed)` reproducibility
  reseed into the caller's subsequent random draws: it now saves the
  caller's `.Random.seed` on entry and restores it via `on.exit()`
  (verified to still restore correctly on every exit path -- an early
  validation error, a mid-bootstrap exhaustion error, a non-default RNG
  kind, and even a "rogue" custom composition estimator that calls
  `set.seed()` itself). `prism()`'s own output is unaffected -- it was
  always fully determined by `seed`, independent of the caller's prior
  RNG state; only the (undocumented, unintended) side effect on code
  that runs *after* `prism()` returns is gone. `seed = NULL` is
  unaffected (nothing to restore -- it was never reseeding).
- Added `return_draws = FALSE`: when `TRUE`, the result gains a `draws`
  component with every draw's lower/upper/relative-covariance
  values (S_effective x n_pairs matrices, aligned via `pair_index`) --
  the data `ci_lower`/`ci_upper`/`pairwise_rel` intervals are already
  computed from, but previously discarded after aggregation. Works
  identically whether or not `stream` was engaged.
- A second, harder adversarial pass (extreme Dirichlet parameters beyond
  the one already fixed, near-boundary PSD/feasibility cases cross-checked
  against an independent CVXR-only formulation, object aliasing/mutation
  safety for reused `scale` objects, feature-name edge cases, nested
  `prism()`-inside-a-custom-estimator calls, and a large `return_draws`
  request) found no further bugs.

# prism 0.3.0

Follow-up round on the 0.2.0 rewrite, based on review feedback:

- Consolidated `scale`, `scale_log`, `scale_log_ci_level`, `scale_log_rho`,
  `scale_log_lower_zero` (5 arguments) into one `scale` argument, mirroring
  the `composition = prism_composition_*()` pattern: `scale` now accepts
  `NULL`, a bare vector/matrix (raw or auto-detected log-scale), a
  [`prism_scale_log()`], or a [`prism_scale_bounds()`]. `sigma_L`/`sigma_U`/
  `rho_L`/`rho_U` remain available as `prism()` arguments as shorthand for
  `prism_scale_bounds()`, forwarded internally. Renamed
  `scale_log_ci_level`/`scale_log_rho`/`scale_log_lower_zero` to
  `scale_ci_level`/`scale_estimate_rho`/`scale_lower_zero`.
- Auto-detects raw vs. already-log-scale for a bare `scale` vector/matrix by
  sign (a raw measurement can't be non-positive); all-positive already-log
  data is genuinely ambiguous and must be wrapped in `prism_scale_log()`
  explicitly -- no magnitude-based guessing.
- `scale_lower_zero` (and the `prism_scale_log()`/`.estimate_scale_log_bounds()`
  default) is now `TRUE`: the lower scale-SD CI endpoint defaults to 0
  rather than the two-sided chi-square lower quantile, since sigma cannot
  be negative and that quantile was asserting a floor that isn't a real
  constraint.
- `S` default raised from 1000 to 2000. `delta` default changed from 0.1 to
  0 (tests significance against a point null rather than a null region, by
  default).
- Reintroduced streaming, simplified relative to the pre-0.2.0 version (no
  longer needs to support single-column high-resolution reads, since that
  feature is gone): `stream = FALSE` by default, but automatically enabled
  when the estimated in-memory cost of all draws exceeds
  `stream_memory_limit` (2 GB default) so a large `D`/`S` cannot silently
  exhaust memory; `stream = TRUE` forces it on. Verified numerically
  identical to non-streamed output for the same seed.
- Fixed a latent bug found while adding streaming support: `.pair_draws()`
  mishandled the `S == 1` (single deterministic draw) case -- `vapply()`
  collapses to a plain vector rather than a matrix when each call returns
  one value, so `t()` transposed the wrong way. This path was previously
  untested because the pre-streaming aggregation code happened to avoid
  calling it when `S == 1`; a refactor to share aggregation logic between
  the streamed and in-memory paths started exercising it and surfaced the
  bug immediately.

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
