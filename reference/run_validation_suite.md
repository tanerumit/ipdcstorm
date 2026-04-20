# Run the three-tier hazard model validation suite

Evaluates a stationary or climate-adjusted hazard model against the
historical record using three complementary checks:

- Tier 1 — Hindcast (return levels):

  Withholds the most recent `holdout_years` of the record, re-fits the
  model on the training period, simulates `n_sim` synthetic years, and
  compares GEV-based return levels against those estimated from the full
  observed record. A location/return-period pair passes (\\\checkmark\\)
  when the simulated return level falls within the observed `conf_level`
  confidence interval.

- Tier 2 — Storm rate check:

  Compares calibrated annual TS and hurricane passage rates per location
  against published HURDAT2/IBTrACS reference bounds. Rates outside the
  reference range are flagged and may indicate a search-radius or
  record-length issue.

- Tier 3 — Wind field spot-checks:

  Compares Holland wind-profile site estimates against historical
  station observations for individual storm passages. Each comparison is
  graded against an observation-quality-dependent tolerance band (A:
  \\\pm\\15%, B: \\\pm\\25%, C: \\\pm\\35%).

Results are printed to the console and, when `cfg$save_plots` or
`cfg$save_tables` is `TRUE`, written to `cfg$out_dir`.

## Usage

``` r
run_validation_suite(out, cfg = make_validation_cfg())
```

## Arguments

- out:

  Named list returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).
  The following elements are used:

  `events`

  :   Tibble of historical storm events with columns `location`, `year`,
      `V_max_kt`, and storm-class identifiers. Required for all three
      tiers.

  `fit`

  :   Tibble of fitted model parameters (lambda, GEV shape, `beta_sst`,
      `gamma_intensity`). Used to replicate the model's climate
      adjustments in the hindcast simulation.

  `sim`

  :   Simulation output; its length sets `n_sim` when `cfg$n_sim` is
      `NULL`.

  `trackpoints`

  :   Named list of track-point tibbles per location, used by the
      wind-field tier. Can be `NULL` if wind-field validation is not
      needed.

  `run_metadata`

  :   List of run-level metadata (seed, data ID, etc.) used for
      provenance labelling in outputs.

- cfg:

  A `validation_cfg` object created by
  [`make_validation_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_validation_cfg.md).
  Controls holdout length, return periods, confidence level, simulation
  size, output paths, and expert GEV bounds. Defaults to
  [`make_validation_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_validation_cfg.md)
  when omitted.

## Value

A named list with the following elements:

- `hindcast`:

  Results from Tier 1. A list with:

  `comparison`

  :   Tibble with one row per `(location, return_period)` combination
      and columns: `location`, `return_period`, `obs_full_rl` (observed
      GEV return level in kt), `sim_rl` (simulated GEV return level in
      kt), `obs_ci_lo` / `obs_ci_hi` (lower/upper bounds of the
      `conf_level` CI around the observed return level), `bias_pct`
      (relative bias in %), `obs_in_ci` (logical pass/fail), `sim_ci_lo`
      / `sim_ci_hi` (CI bounds from the simulation-side bootstrap, for
      reference), and `obs_ci_method` (`"delta"` or `"bootstrap"`).

  `per_island`

  :   Named list of per-location diagnostic objects, each containing
      training/test year splits, observed and simulated GEV fits, and
      bias decomposition.

- `rate_check`:

  Tibble from
  [`validate_rates()`](https://tanerumit.github.io/ipdcstorm/reference/validate_rates.md)
  with one row per `(location, storm_class)` and columns: `location`,
  `storm_class`, `lambda_fitted`, `ref_lo`, `ref_hi`, `flag` (`"OK"` or
  `"FLAG"`).

- `wind_field`:

  Tibble from
  [`validate_wind_field()`](https://tanerumit.github.io/ipdcstorm/reference/validate_wind_field.md)
  with one row per `(storm, location)` spot-check and columns:
  `storm_name`, `location`, `obs_1min_equiv_kt`, `model_V_site_kt`,
  `bias_kt`, `bias_pct`, `obs_quality`, `bias_threshold_pct`, `bias_ok`.
  `NULL` when no matching observations are found.

- `summary`:

  Compact tibble with one row per location aggregating pass rates across
  return periods (Tier 1), rate-check flags (Tier 2), and wind-field
  tolerances (Tier 3).

- `climate_response`:

  List of climate-sensitivity diagnostics derived from the fitted model,
  including basin count ratios, hurricane-fraction shifts, and intensity
  redistribution metrics.

- `artifacts`:

  Named list with two sub-lists:

  `plots`

  :   Named character vector of PNG file paths written when
      `cfg$save_plots = TRUE`: `hindcast_return_levels`,
      `bias_decomposition`, `wind_field_scatter`, `qq_annual_max`,
      `exceedance_comparison`.

  `tables`

  :   Named character vector of CSV file paths written when
      `cfg$save_tables = TRUE`: `hindcast_csv`, `rate_check_csv`,
      `wind_field_csv`, `summary_csv`.

## See also

[`make_validation_cfg`](https://tanerumit.github.io/ipdcstorm/reference/make_validation_cfg.md)
for configuring the validation run,
[`run_hazard_model`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
for producing the required `out` object,
[`plot_hindcast_validation`](https://tanerumit.github.io/ipdcstorm/reference/plot_hindcast_validation.md)
for standalone hindcast plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run with default settings
val_out <- run_validation_suite(out_base)

# Custom configuration: 15-yr holdout, 95% CI, save outputs
val_cfg <- make_validation_cfg(
  holdout_years  = 15L,
  return_periods = c(5, 10, 25, 50),
  conf_level     = 0.95,
  save_plots     = TRUE,
  out_dir        = "output/validation"
)
val_out <- run_validation_suite(out_base, cfg = val_cfg)

# Inspect Tier 1 return-level comparison
val_out$hindcast$comparison

# Overall pass rate
mean(val_out$hindcast$comparison$obs_in_ci, na.rm = TRUE)
} # }
```
