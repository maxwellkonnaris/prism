# PRISM package audit: 2026-07-11

## Scope

This audit treats PRISM as a reusable R package, not as an analysis directory.
The review covered statistical branch selection, compositional uncertainty,
scale uncertainty, high-resolution refinement, streaming, output stability,
side effects, dependencies, tests, documentation, and package installation.

## Baseline facts

Before restructuring, the package source included benchmark frameworks,
cross-cohort analysis, plotting pipelines, several unrelated covariance
estimators, executables, serialized objects, and roughly 99 MB analysis data.
The installed package was 95.7 MB. `R CMD check` reported two errors, six
warnings, and two notes. The core method tests passed, while cross-cohort tests
failed independently of PRISM.

All displaced material was preserved outside the package under:

```text
analysis/prism_package_archive_2026-07-11/
```

## Package boundary after restructuring

```text
R/
  bounds.R                    covariance equations and sharp refinement
  composition-estimators.R    public estimator specifications and plug-in checks
  composition-mln.R           optional MLN fit and tuning implementation
  prism.R                     short documented public interface
  engine.R                    validation, draw processing, streaming, inference
  result.R                    prism_result construction and printing
  utils.R                     compact numerical helpers
```

Only PRISM, covariance bounds, composition-estimator constructors, and the two
high-resolution plot functions are exported.

## Statistical decisions

- The default composition estimator is Dirichlet-multinomial.
- MLN is optional and is not part of the default estimator.
- Fixed composition is explicit and adds its pseudocount before closure/logging.
- A zero pseudocount is allowed only when no zero logarithm or zero Dirichlet
  parameter can occur.
- Rho bounds without sigma bounds are rejected.
- Sigma and rho domains are enforced.
- Fixed rho compatibility is checked through an augmented PSD matrix.
- Diagonal bounds use rho when rho information is supplied.
- High-resolution refinement occurs within selected draws before endpoint
  aggregation and full-family BH adjustment.

## Performance decisions

- Streaming reconstructs the same final matrices and tables as the full path.
- Stream files are temporary and removed on exit.
- Dirichlet compositions are generated on demand rather than held as a
  `D x N x S` array.
- Raw-scale uncertainty uses an independent inner draw budget rather than an
  `S`-sized inner simulation inside every one of `S` outer draws.
- Returning relative draws remains an explicit memory-expensive request.
- MLN posterior arrays remain estimator-managed because `fido` returns them as
  a fitted posterior object.

### Local profile

The reproducible command in `tools/benchmark-memory.R` used `D=40`, `N=80`,
`S=40` with Dirichlet-multinomial draws. On the audit machine it produced:

| mode | elapsed seconds | largest profiled allocation MB | result size MB |
|---|---:|---:|---:|
| full | 0.810 | 0.488 | 0.140 |
| streamed | 1.335 | 0.125 | 0.141 |

This is a small local profile, not a general benchmark. It supports the specific
claim that streaming reduced the largest allocation observed by `Rprofmem` in
this run, at the cost of temporary-file I/O and longer elapsed time. It does not
measure process peak resident memory.

## Remaining uncertainty

The exact frequentist calibration of the final intervals depends on the
composition estimator, bootstrap assumptions, and supplied scale uncertainty.
The high-resolution solver is tested for numerical tightening, but solver
behavior still depends on installed CVXR backends and conditioning. MLN runtime
and conditioning cannot be guaranteed from dimensions alone; its Hessian memory
guard is a preflight estimate, not a complete runtime predictor.

## Evidence that would change this audit

- a simulation showing systematic undercoverage under the intended data model;
- a full-versus-stream mismatch under a fixed seed;
- a high-resolution interval wider than its conservative input beyond tolerance;
- package-check failures after installation in a clean library;
- profiling that identifies a different dominant allocation or runtime kernel.
