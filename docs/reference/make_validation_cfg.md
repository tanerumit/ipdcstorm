# Create a validation configuration

Creates a typed configuration object for
[`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md).
Follows the same pattern as
[`make_hazard_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md)
and `make_sst_cfg()`: common parameters up front, expert-only knobs in
`advanced`.

## Usage

``` r
make_validation_cfg(
  holdout_years = 10L,
  n_sim = NULL,
  return_periods = c(5, 10, 25, 50),
  conf_level = 0.95,
  tier1_pass_rate = 0.8,
  full_output = FALSE,
  seed = 42L,
  out_dir = "output/validation",
  save_plots = TRUE,
  save_tables = TRUE,
  advanced = NULL
)
```

## Arguments

- holdout_years:

  Integer; number of years to hold out from the end of the historical
  record for train/test split (default: 10).

- n_sim:

  Integer or `NULL`; number of synthetic years to simulate for hindcast
  comparison. If `NULL` (default),
  [`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md)
  inherits the simulation length from the hazard model output.

- return_periods:

  Numeric vector of return periods (years) to compare (default: 5, 10,
  25, 50).

- conf_level:

  Numeric; confidence level for return-level CIs (default: 0.90). Must
  be in `(0, 1)`. This controls the width of the confidence interval
  around GEV-based return level estimates, computed via the delta method
  or parametric bootstrap of the hurdle-GEV model. Common choices: 0.90
  (standard), 0.95 (stricter), 0.80 (exploratory). **Note**: this
  quantifies estimation uncertainty in return levels given the observed
  record length — it is not a prediction interval for future storm
  intensities.

- tier1_pass_rate:

  Numeric; minimum fraction of return-period/location comparisons that
  must fall within the CI for Tier 1 to show ✓ in the summary (default:
  0.80). A value of 1.0 restores the original strict all-must-pass
  behaviour. Must be in `(0, 1]`.

- full_output:

  Logical; when `FALSE` (default)
  [`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md)
  returns a concise list (`summary`, `comparison`, `rate_check`,
  `wind_field`, `plots`). When `TRUE`, the full internal objects are
  also returned (`hindcast`, `climate_response`, `artifacts`). The slim
  output is sufficient for most diagnostics; use `full_output = TRUE`
  when you need the per-island breakdown or to re-generate plots from
  the returned object.

- seed:

  Integer; random seed for reproducibility.

- out_dir:

  Character; output directory for saved plots and tables.

- save_plots:

  Logical; whether to save standard validation figures.

- save_tables:

  Logical; whether to save CSV + markdown tables.

- advanced:

  Optional named list of expert parameters. Most users should leave this
  as `NULL`. Supported names:

  `xi_bounds`

  :   Numeric vector of length 2; allowed range for GEV shape parameter
      (default: `c(-0.3, 0.4)`).

  `base_size`

  :   Numeric; base font size for ggplot themes (default: 11).

  `hindcast_use_raw_rates`

  :   Logical; whether Tier 1A hindcast uses raw fitted rates instead of
      adjusted rates (default: `TRUE`). This is the validated package
      default because adjusted lambda worsened hindcast bias at Saba and
      Statia.

  The default validation path is the current best validated overall
  tradeoff: Tier 1A hindcast uses raw lambda, the package runs the
  legacy wind-field path unless an internal option overrides it, and the
  hindcast sampler defaults to `legacy`. Experimental wind-field
  comparison modes remain internal and are never activated by
  `make_validation_cfg()` defaults.

## Value

A list with class `c("validation_cfg", "list")`.

## Examples

``` r
# Defaults - suitable for most users
val_cfg <- make_validation_cfg()

# Custom holdout and output location
val_cfg <- make_validation_cfg(holdout_years = 15, out_dir = "results/val")

# Expert tuning
val_cfg <- make_validation_cfg(
  n_sim = 10000,
  advanced = list(xi_bounds = c(-0.4, 0.5), base_size = 13)
)
```
