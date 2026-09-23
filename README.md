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
  composition = prism_composition_dirichlet(pseudocount = 0.5, concentration = 1),
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
unless resampled together. Requested and completed draw counts are returned
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
prism_composition_dirichlet(pseudocount = 0.5, concentration = 1)
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
| numeric vector or matrix | estimated sigma and rho from log-transformed measurements | all values are assumed to already be on the log scale; technical replicates are bootstrapped within subject and averaged |
| [`prism_scale_log()`] | the same log-scale data with non-default estimation settings | technical replicates are bootstrapped within subject and averaged |
| [`prism_scale_bounds()`], or the `sigma_L`/`sigma_U`/`rho_L`/`rho_U` shorthand | fixed bounds, no data | sharp bounds after a compatibility check; genuine rho intervals use exact CVXR optimization over the feasible ellipsoid-box intersection |

Scale measurements must be supplied on the log scale. Bare numeric vectors and
matrices are accepted under that contract; PRISM never applies a logarithm or
attempts to infer a transformation from the values. Bare input uses point
estimates of sigma and rho within each bootstrap draw by default. Use
`prism_scale_log()` only to change those settings. Setting `ci_level > 0` adds
chi-square/Fisher-z intervals inside every draw and is restricted to the
explicit `experimental_ci = TRUE` sensitivity-analysis mode.

A constant log-scale vector implies `sigma = 0`; PRISM returns
`sigma_L = sigma_U = 0` and does not estimate rho because rho is undefined and
irrelevant when multiplied by zero. PRISM still errors when any retained
log-composition feature has zero or numerically zero marginal variance.

Supplying `rho_L`/`rho_U` without `sigma_L`/`sigma_U` is an error. Sigma
bounds are checked against `[0, Inf)` and rho bounds against `[-1, 1]`. A
fixed rho vector must keep the augmented covariance matrix
`[[Sigma_rel, u], [u', 1]]` (with `u = sqrt(diag(Sigma_rel)) * rho`) positive
semidefinite. A genuine rho interval must contain at least one such feasible
point in the box; this is certified once via a closed-form shortcut, a
supplied `rho_witness`, an L-BFGS-B search (`Sigma_rel` full rank), or CVXR
(`Sigma_rel` rank-deficient and no witness works).

For a genuine rho interval, PRISM then uses CVXR to compute the exact support
values

```text
m_ij = min (s_i * rho_i + s_j * rho_j)
M_ij = max (s_i * rho_i + s_j * rho_j)
```

over the same globally feasible ellipsoid-box intersection. These values are
propagated through the bounded-sigma quadratic to obtain sharp entrywise
covariance endpoints. CVXR is therefore required for genuine rho intervals,
but not for fixed rho, bounded sigma without rho restrictions, or unbounded
scale. Exact interval-rho calculation solves two conic support problems for
each upper-triangle matrix entry and can be substantially slower for large
feature sets or many bootstrap draws.

If `sigma_U = 0`, all scale terms vanish and PRISM returns `Sigma_rel`
directly; rho feasibility is irrelevant in this degenerate-scale case.

## Result object

`prism()` returns a `prism_result` with:

- `ci_lower`, `ci_upper`: symmetric absolute-covariance CI endpoint matrices;
- `pairwise`: off-diagonal pairs with CI endpoints, `covers_zero`, bootstrap
  positive/negative/covers-zero counts and frequencies, and p/q-values when
  paired subject bootstrapping is active;
- `pairwise_rel`: relative log-covariance CIs, including the
  diagonal;
- `diagnostics`: data (dimensions, dropped samples/features) and sampling
  (bootstrap draw and streaming bookkeeping);
- `parameters`: effective estimator and computational settings.

Without paired subject bootstrapping, `p_value`/`q_value` and the bootstrap
sign frequencies are `NA`; `covers_zero` still reports directly whether the
aggregated identification interval contains zero. With bootstrapping, the
sign summaries count intervals wholly above zero, wholly below zero, or
covering zero. The two-sided p-value against `[-delta, delta]` is computed from
the opposing bootstrap tails of the per-draw identification intervals.

## Principal failure modes

- zeros with a zero pseudocount;
- non-integer values supplied as counts;
- fewer than two valid samples or features after filtering;
- rho information without sigma information;
- infeasible fixed correlations or an empty rho-box/ellipsoid intersection;
- zero or numerically zero marginal variance in a retained log-composition
  feature;
- custom composition draws with invalid dimensions or values;
- unavailable CVXR or conic solver for exact genuine-rho-interval bounds.

The conclusion should be reconsidered if composition draws are miscalibrated,
scale bounds are not scientifically defensible, or bootstrap samples are not
exchangeable.
