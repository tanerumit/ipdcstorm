# Compute return level CIs via parametric bootstrap

Generates bootstrap CIs by:

1.  Resampling the annual maxima

2.  Refitting the hurdle-GEV to each bootstrap sample

3.  Computing return levels from each fit

## Usage

``` r
bootstrap_return_level_ci(
  annual_max,
  return_periods = c(5, 10, 25, 50),
  n_boot = 500,
  xi_bounds = c(-0.3, 0.4),
  conf_level = 0.95
)
```

## Arguments

- annual_max:

  Numeric vector of annual maxima (including zeros).

- return_periods:

  Numeric vector of return periods.

- n_boot:

  Integer; number of bootstrap replicates.

- xi_bounds:

  Bounds on GEV shape.

- conf_level:

  Numeric; confidence level for CIs (default: 0.90).

## Value

Tibble with columns: `return_period`, `sim_median`, `sim_ci_lo`,
`sim_ci_hi`, `sim_lo_50`, `sim_hi_50`. The `ci_lo`/`ci_hi` columns
correspond to the specified `conf_level`.
