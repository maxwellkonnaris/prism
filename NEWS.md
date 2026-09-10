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
