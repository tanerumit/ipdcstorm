## 1. Goal
Isolate which bundled change(s) worsened hindcast positive bias, restore the valid Tier 1A default to raw lambda, and keep only behavior-neutral R34 provenance diagnostics.

## 2. Scope
Included:
- Reverted `make_validation_cfg()` Tier 1A default to raw rates.
- Restored legacy wind-field fallback behavior as the package default.
- Restored legacy hindcast sampler behavior as the package default.
- Added internal workspace-only attribution helpers for wind/rate/sampler reruns.
- Updated targeted tests and generated docs.

Excluded:
- No public API renames or signature changes.
- No retuning of storm physics constants beyond restoring pre-bundled defaults.
- No pipeline/script edits outside the allowed task files.

## 3. Problem solved
Recent bundled changes mixed three levers at once, so hindcast deterioration at Saba and Statia could not be attributed cleanly. Tier 1A also defaulted to adjusted lambda, which invalidated comparison to the pre-change workspace path.

## 4. Summary
The package default path now matches the pre-bundled workspace behavior for Tier 1A rate mode, wind fallback behavior, and sampler behavior. `R34_source` is preserved as cheap provenance for diagnostics, and an internal helper `.run_hindcast_attribution_grid()` reruns the workspace model under controlled wind/rate/sampler modes with recorded seed, data identifier, and parameter identifier metadata.

Controlled reruns for `Saba`, `Statia`, and `St_Martin` showed that adjusted lambda was the dominant driver of worse positive RL bias at Saba and Statia. The diagnostic new wind-field fallback also worsened positive bias at both sites, but less than the rate-mode change. The bounded sampler effect was smaller and secondary.

## 5. Files changed
- `R/hazard_core.R` — 85 additions, 23 deletions. Restored legacy default wind behavior; kept diagnostic `R34_source`; added internal wind-mode toggle for attribution reruns.
- `R/hazard_validation.R` — 266 additions, 33 deletions. Restored raw-rate default; restored legacy sampler default; added internal attribution grid and metadata summaries.
- `man/make_validation_cfg.Rd` — regenerated to show the raw-rate default.
- `tests/testthat/test-hazard-validation-exported.R` — 83 additions. Updated default expectations; added bounded-sampler diagnostic test; added attribution-grid metadata test.
- `tests/testthat/test-hazard_outer_cutoff.R` — 13 additions, 8 deletions. Restored legacy default cutoff expectations.
- `tests/testthat/test-rmw_estimation.R` — 63 additions, 1 deletion. Updated wind-mode expectations so tighter fallback behavior is diagnostic-only.

## 6. Commands run
- `Rscript -e "parse(file='R/hazard_core.R'); parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='tests/testthat/test-hazard-validation-exported.R'); parse(file='tests/testthat/test-hazard_outer_cutoff.R'); parse(file='tests/testthat/test-rmw_estimation.R')"`
- `Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|outer_cutoff|rmw)|hindcast|wind')"`
- `Rscript -e "devtools::document()"`
- `Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|outer_cutoff|rmw)|hindcast|wind')"`
- `Rscript -e "pkgload::load_all('.'); targets <- tibble::tribble(~location, ~lat, ~lon, 'St_Martin', 18.0708, -63.0501, 'Saba', 17.6350, -63.2300, 'Statia', 17.4890, -62.9740); hazard_cfg <- make_hazard_cfg(data_path = 'inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv', search_radius_km = 800, historical_start_year = 1970L, simulation_years = 500L); validation_cfg <- make_validation_cfg(holdout_years = 10L, n_sim = 500L, return_periods = c(5,10,25,50), conf_level = 0.90, seed = 42L, save_plots = FALSE, save_tables = FALSE); res <- ipdcstorm:::.run_hindcast_attribution_grid(cfg = hazard_cfg, targets = targets, validation_cfg = validation_cfg, climate = make_climate_cfg(scenario = 'stationary'), locations = c('Saba','Statia','St_Martin'), model_seed = 42L); print(as.data.frame(res[['summary']][, c('location','wind_field_mode','annual_rate_mode','sampler_mode','obs_xi','sim_xi','rl_bias_rp5','rl_bias_rp10','rl_bias_rp25','rl_bias_rp50')]), row.names = FALSE)"`

## 7. Test results
- Parse checks: passed for touched R files and targeted tests.
- Targeted tests: passed, `69` assertions in the requested validation/hindcast/wind contexts, `0` failures, `0` warnings, `0` skips.
- Documentation regeneration: passed.
- Controlled attribution rerun: completed from `pkgload::load_all('.')` workspace code path.

## 8. Behavior changes
- `make_validation_cfg()` now defaults `advanced$hindcast_use_raw_rates = TRUE` again.
- Default wind-field execution path again uses the legacy fallback behavior.
- Default hindcast intensity sampling again uses the legacy sampler behavior.
- `R34_source` remains available for diagnostics/provenance.
- Internal helper `.run_hindcast_attribution_grid()` now produces comparable workspace reruns with:
  - `wind_field_mode`
  - `annual_rate_mode`
  - `sampler_mode`
  - `model_seed`
  - `validation_seed`
  - `data_id`
  - `parameter_id`
  - `lambda_scaler_id`

## 9. Follow-ups/risks
- The attribution helper is internal and intentionally not part of the public API; downstream scripts should call it via `ipdcstorm:::` only for diagnostics.
- `R34_source` is preserved, but only the legacy pathway is active by default; the tighter fallback path is diagnostic-only until stronger methodological justification exists.
- St. Martin uses a slightly different training split between legacy and diagnostic-new wind reruns because the wind mode changes which years retain TS+ events; that is visible in the controlled output and should be considered when interpreting smaller deltas there.

