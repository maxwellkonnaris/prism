# prism

`prism` estimates identification regions and confidence intervals for absolute
log-covariances from compositional observations and optional scale information.

## Statistical target

For feature `d`, write absolute log abundance as

```text
log W_d = log P_d + U,
```

where `P_d` is relative abundance and `U` is log total scale. For pair `(i,j)`,

```text
Cov(log W_i, log W_j)
  = Sigma_rel[i,j] + sigma^2
    + sigma * (s_i * rho_i + s_j * rho_j),
```

where:

- `Sigma_rel` is the covariance of sample log-compositions;
- `sigma = SD(U)`;
- `s_d = sqrt(Sigma_rel[d,d])`;
- `rho_d = Cor(log P_d, U)`.

The target is the population absolute log-covariance. `Sigma_rel` and any
sample-level scale summaries are estimated from observed samples. The package
reports identification bounds because the absolute covariance is generally not
point identified from compositions alone.

## Minimal use

```r
library(prism)

fit <- prism(
  counts = count_matrix,
  composition = prism_composition_dirichlet(pseudocount = 1, concentration = 1),
  bootstrap = TRUE,
  S = 2000,
  sigma_L = 0.2,
  sigma_U = 0.8,
  seed = 1
)

fit
fit$pairwise
fit$diagnostics
```

The count matrix must be features by samples, non-negative, and
integer-valued. `bootstrap` accepts `TRUE`/`FALSE` or `"both"`/`"none"`.
`"both"` is a paired sample bootstrap: one resampled sample index is applied
to the count composition and its matched scale measurement, and all empirical
covariance parameters are re-estimated from that draw. Independently
resampling only the composition or only the scale block is not supported,
since those blocks need not jointly define a positive-semidefinite covariance
unless resampled together. Attempted, accepted, and rejected counts;
acceptance and rejection fractions; and rejection-reason counts are returned
in `fit$diagnostics$sampling`.

Draws are held in memory by default, but streamed to temporary disk files
when the estimated memory cost of holding all of them exceeds
`stream_memory_limit` (2 GB by default) -- automatically, whether or not
`stream = TRUE` was requested -- so a large `D`/`S` does not silently exhaust
memory. Pass `stream = TRUE` to force it on for any size, or raise
`stream_memory_limit` to opt out of the safety margin.

## Composition estimators

The default is:

```r
prism_composition_dirichlet(pseudocount = 1, concentration = 1)
```

Independent sample-level Dirichlet distributions with parameters
`concentration * (counts[, n] + pseudocount)`.

Fixed observed composition (no composition uncertainty) is requested
explicitly:

```r
prism_composition_fixed(pseudocount = 0.5)
```

A custom estimator has two functions:

```r
my_estimator <- prism_composition_estimator(
  name = "my_model",
  fit = function(counts, n_draws, seed, verbose) {
    # Called once, on the full (unresampled) counts. Return fitted state.
  },
  draw = function(fit, draw) {
    # Return one finite, positive D x N composition matrix for draw index `draw`.
  }
)
```

Every draw is checked and reclosed before covariance estimation. A custom
estimator is not resampled at the sample level by PRISM's bootstrap (its
`fit()`/`draw()` own all compositional uncertainty); scale-side resampling,
if `scale` carries sample-level data, still applies.

## Scale information

`scale` is the single argument for "what do I know about total scale?" --
it accepts, in increasing order of how much it specifies:

| `scale` | Regime | Result |
|---|---|---|
| `NULL` (default), no `sigma_L`/etc. either | unbounded scale | finite lower bound, infinite upper bound |
| a bare vector (length N) or matrix (N rows, replicate columns) | estimated sigma and rho | point or bounded values according to `scale_ci_level`, with replicate columns resampled per draw when present |
| [`prism_scale_log()`] | estimated sigma and rho, values already logged | same as above, but the vector/matrix is not re-logged, and its own `ci_level`/`estimate_rho`/`lower_zero` apply |
| [`prism_scale_bounds()`], or the `sigma_L`/`sigma_U`/`rho_L`/`rho_U` shorthand | fixed bounds, no data | closed-form bounds (sharp if `rho_L == rho_U`, conservative otherwise), after a compatibility check |

A bare vector/matrix is assumed raw (un-logged), and PRISM takes the log
itself -- *unless* any value is non-positive, which a raw scale measurement
can never be, in which case the whole input is treated as already log-scale.
All-positive values that are already log-scale are genuinely ambiguous by
value alone and must be wrapped explicitly with `prism_scale_log()`; there is
no reliable magnitude-based heuristic for this and PRISM does not attempt
one. `scale_ci_level`/`scale_estimate_rho`/`scale_lower_zero` configure the
within-draw estimate for a bare vector/matrix; they are an error when
combined with `prism_scale_log()` (which carries its own settings) or with
`prism_scale_bounds()`/`sigma_L` etc. (fixed, nothing to estimate).

Supplying `rho_L`/`rho_U` without `sigma_L`/`sigma_U` is an error. Sigma
bounds are checked against `[0, Inf)` and rho bounds against `[-1, 1]`. A
fixed rho vector must keep the augmented covariance matrix
`[[Sigma_rel, u], [u', 1]]` (with `u = sqrt(diag(Sigma_rel)) * rho`) positive
semidefinite. A genuine rho interval must contain at least one such feasible
point in the box; this is certified once via a closed-form shortcut, a
supplied `rho_witness`, an L-BFGS-B search (`Sigma_rel` full rank), or CVXR
(`Sigma_rel` rank-deficient and no witness works -- an optional dependency,
only needed in that case).

The off-diagonal bound for a genuine (non-fixed) rho interval is a valid but
conservative closed-form outer bound, not the tightest possible one; check
`bounds$sharp$off_diagonal` (from [cov_bounds()]) to see when this applies.

## Result object

`prism()` returns a `prism_result` with:

- `ci_lower`, `ci_upper`: symmetric absolute-covariance CI endpoint matrices;
- `pairwise`: off-diagonal tested pairs with CI endpoints, p-value, and
  BH-adjusted q-value (`p_value`/`q_value` are `NA` when there is only one
  deterministic draw -- see below);
- `pairwise_rel`: relative log-covariance estimates and CIs, including the
  diagonal;
- `diagnostics`: data (dimensions, dropped samples/features) and sampling
  (bootstrap accept/reject bookkeeping);
- `parameters`: effective estimator and computational settings.

With no bootstrap resampling and a deterministic (non-random) composition
estimator, PRISM evaluates the bound once rather than `S` times; `ci_lower`/
`ci_upper` are then the identification interval itself and `p_value` is `NA`
(there is no bootstrap distribution to test against the null region).

## Principal failure modes

- zeros with a zero pseudocount;
- non-integer values supplied as counts;
- fewer than two valid samples or features after filtering;
- rho information without sigma information;
- infeasible fixed correlations or an empty rho-box/ellipsoid intersection;
- excessive rejection when paired bootstrap resampling produces degenerate
  (e.g. zero-variance) relative covariance under a fixed rho constraint;
- too few valid Monte Carlo draws;
- custom composition draws with invalid dimensions or values;
- unavailable CVXR solver for a rank-deficient rho-box feasibility check.

The conclusion should be reconsidered if composition draws are miscalibrated,
scale bounds are not scientifically defensible, or bootstrap samples are not
exchangeable.
