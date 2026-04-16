# Resolve climate scalars for the hazard model

Calibrates historical baseline sensitivities from annual counts and MDR
SST anomalies, resolves one canonical `delta_sst`, and returns flat
simulation-ready climate scalars for
[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

## Usage

``` r
resolve_climate_inputs(
  climate_cfg,
  annual_counts = NULL,
  lambda_table = NULL,
  min_year = 1970L,
  verbose = TRUE
)
```

## Arguments

- climate_cfg:

  List from
  [`make_climate_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_climate_cfg.md).

- annual_counts:

  Optional tibble of annual counts for baseline sensitivity estimation.

- lambda_table:

  Optional tibble from
  [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)
  for historical `p_hur_base`.

- min_year:

  Integer; passed to estimation functions.

- verbose:

  Logical.

## Value

A list with:

- delta_sst:

  Resolved SST shift in degC.

- beta_sst:

  Final rate sensitivity used in simulation.

- gamma:

  Final intensity sensitivity used in simulation.

- p_hur_base:

  Historical baseline hurricane fraction.

- beta_0:

  Historical baseline rate sensitivity.

- gamma_0:

  Historical baseline intensity sensitivity.

- scenario:

  Resolved scenario name.

- input_mode:

  Climate input mode: scenario helper or direct delta_sst.

- target_year:

  Target year used to derive scenario-helper `delta_sst`, if any.

- sensitivity_mode:

  Resolved sensitivity mode.

- k_beta:

  Configured linear shift coefficient for `beta_0`.

- k_gamma:

  Configured linear shift coefficient for `gamma_0`.

- beta_sst_raw:

  Raw rate sensitivity after any sensitivity shift, before basin-rate
  damping is applied.

- rate_scale:

  Applied SST-driven count multiplier used by the run.

- raw_rate_scale:

  Raw SST-driven count multiplier before damping.

- annual_count_series:

  Basin-consistent annual total count series used for `beta_0`
  calibration, in storms/year.

- annual_count_source:

  Character label describing the provenance of `annual_count_series`.

- beta_guardrail:

  List describing any `beta_0` plausibility fallback.

- sst_scale_guardrail:

  List describing any multiplier guardrail applied to `beta_sst` for the
  resolved scenario.

- baseline_years:

  SST anomaly baseline years.

- perturb:

  Resolved storm-perturbation parameters (or `NULL`).

- perturb_state:

  Storm-perturbation state label.

- source:

  SST source used for calibration.

## Examples

``` r
cfg <- make_climate_cfg(sst_source = "builtin", scenario = "stationary")
climate <- resolve_climate_inputs(cfg, verbose = FALSE)
climate$beta_sst
#> [1] 0.6
```
