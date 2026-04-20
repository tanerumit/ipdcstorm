# Build a climate configuration object for the hazard model

Creates a climate configuration object for
[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).
The climate workflow first calibrates historical baseline sensitivities
(`beta_0`, `gamma_0`) from historical annual counts and MDR SST
anomalies, then resolves scenario-specific simulation scalars from
`delta_sst`. By default the historical sensitivities are used unchanged.
Expert users may instead enable a linear sensitivity shift:

`beta_sst = beta_0 * (1 + k_beta * delta_sst)`

`gamma = gamma_0 * (1 + k_gamma * delta_sst)`

Optional storm perturbation shifts individual event intensity, duration,
and wind radii by `delta_sst`. It is disabled by default; pass
`perturb = TRUE` to enable scenario-appropriate defaults.

Climate can be specified in two equivalent ways:

1.  Scenario helper mode: provide `scenario` plus `target_year` to
    derive `delta_sst`.

2.  Direct mode: provide `delta_sst` explicitly.

After `delta_sst` is resolved, downstream hazard logic depends only on
the resolved `delta_sst` and fixed model settings. Scenario names and
timing are retained only as provenance metadata. If both direct
`delta_sst` and scenario-derived inputs are supplied, they must resolve
to the same `delta_sst` or an error is raised.

If `scenario = "stationary"` and `delta_sst` is not supplied, the
returned config represents the canonical baseline hazard-model
specification with `delta_sst = 0`. Climate is always resolved through
this configuration; there is no disabled climate mode.

## Usage

``` r
make_climate_cfg(
  scenario = "stationary",
  sst_source = c("builtin", "csv", "ersst_nc"),
  sst_path = NULL,
  baseline_years = 1991L:2020L,
  start_year = 2025L,
  delta_sst = NULL,
  target_year = NULL,
  sensitivity_mode = c("fixed", "linear_shifted"),
  k_beta = 0,
  k_gamma = 0,
  perturb = NULL
)
```

## Arguments

- scenario:

  Optional character; climate scenario name used for helper resolution
  and provenance metadata.

- sst_source:

  Character; one of "builtin", "csv", "ersst_nc".

- sst_path:

  Optional character; path to SST data file (CSV or NetCDF).

- baseline_years:

  Integer vector; years for climatological baseline.

- start_year:

  Integer; first year of the simulation scenario.

- delta_sst:

  Optional numeric scalar explicit SST shift in degC.

- target_year:

  Optional numeric scalar target year used to derive scenario-helper
  `delta_sst`.

- sensitivity_mode:

  Character scalar; one of `"fixed"` or `"linear_shifted"`.

- k_beta:

  Numeric scalar; linear sensitivity-shift coefficient applied to
  `beta_0` when `sensitivity_mode = "linear_shifted"`.

- k_gamma:

  Numeric scalar; linear sensitivity-shift coefficient applied to
  `gamma_0` when `sensitivity_mode = "linear_shifted"`.

- perturb:

  Controls storm-perturbation of individual event properties (intensity,
  duration, wind radii) by `delta_sst`.

  - `TRUE`: enable with scenario-appropriate defaults (recommended for
    KNMI scenarios).

  - `FALSE` or `NULL`: disable (default).

  - named list: enable with explicit parameter overrides (see
    [`default_cc_params()`](https://tanerumit.github.io/ipdcstorm/reference/default_cc_params.md)).

## Value

A list with class "climate_cfg" containing climate configuration
parameters.

## Examples

``` r
cfg <- make_climate_cfg(
  sst_source = "builtin",
  scenario = "ssp245",
  target_year = 2050
)
cfg_direct <- make_climate_cfg(delta_sst = 0.5)
cfg
#> Climate configuration
#>   Input mode        : scenario helper
#>   Scenario          : ssp245
#>   SST source        : builtin
#>   Baseline          : 1991-2020
#>   Target year       : 2050.0
#>   Climate mode      : future climate run
#>   Sensitivity mode  : fixed
#>   Shift coefficients: fixed historical sensitivities
#>   Storm perturbation: disabled
```
