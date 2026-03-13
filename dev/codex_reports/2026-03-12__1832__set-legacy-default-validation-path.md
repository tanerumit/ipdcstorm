## 1. Goal

Set the package default validation path back to the best validated legacy configuration and document that this is the current operational default.

## 2. Scope

Included:
- Explicit default-path documentation for raw lambda, legacy wind path, and legacy sampler.
- Targeted tests that assert the restored legacy defaults.
- Default-path validation rerun from the workspace code path.

Not included:
- Any new bias-reduction tuning.
- Any change to package defaults beyond reaffirming the validated legacy path.
- Any promotion of experimental observed-R34 candidate modes into normal execution.

## 3. Problem solved

Recent diagnostic work improved attribution, but the trial fixes for the remaining Saba and Statia upper-tail bias were not acceptable. The package needed its ordinary validation path clearly anchored to the last configuration that still gave the best overall hindcast tradeoff without silently depending on internal experimental modes.

## 4. Summary

The default validation path is now explicitly documented and tested as:
- Tier 1A hindcast uses raw lambda by default.
- The wind-field path defaults to `legacy`.
- The hindcast sampler defaults to `legacy`.

Experimental observed-R34 comparison modes remain internal-only and are not activated unless explicitly requested through internal options. The workspace validation rerun confirmed the default path still matches the validated legacy baseline, with residual upper-tail bias at Saba and Statia still present and intentionally deferred to a separate methodological phase.

## 5. Files changed

- `R/hazard_validation.R` (`+1359 / -12` in current worktree delta): added and retained internal hindcast diagnostics from prior forensic work; this pass added explicit default-path documentation for raw lambda plus legacy wind/sampler behavior.
- `man/make_validation_cfg.Rd` (`+9 / -1`): documented raw lambda as the validated default and described the full default validation path.
- `man/dot-estimate_site_wind_holland.Rd` (`+10 / -0`): clarified that observed-R34 experimental modes are internal-only and non-default.
- `tests/testthat/test-hazard-validation-exported.R` (`+172 / -1`): asserted raw-lambda, legacy wind-mode, and legacy sampler defaults explicitly while preserving targeted deterministic diagnostics coverage.
- `tests/testthat/test-hazard_outer_cutoff.R` (`+18 / -0`): added a guard test ensuring ordinary runs stay on the validated legacy wind path and sampler path.

## 6. Commands run

```powershell
Rscript -e "parse(file='R/hazard_validation.R'); parse(file='tests/testthat/test-hazard-validation-exported.R'); parse(file='tests/testthat/test-hazard_outer_cutoff.R'); parse(file='tests/testthat/test-rmw_estimation.R')"
```

```powershell
Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|outer_cutoff|rmw)|hindcast|wind')"
```

```powershell
Rscript -e "tools::checkRd('man/make_validation_cfg.Rd'); tools::checkRd('man/dot-estimate_site_wind_holland.Rd'); cat('Rd OK\n')"
```

```powershell
Rscript -e "pkgload::load_all('.'); library(dplyr); targets <- tibble::tribble(~location, ~lat, ~lon, 'St_Martin', 18.0708, -63.0501, 'Saba', 17.6350, -63.2300, 'Statia', 17.4890, -62.9740, 'Puerto_Rico', 18.2208, -66.5901, 'Miami', 25.7617, -80.1918); cfg <- make_hazard_cfg(data_path='inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv', search_radius_km=600, historical_start_year=1970L, simulation_years=2000L); validation_cfg <- make_validation_cfg(holdout_years=10L, n_sim=2000L, return_periods=c(5,10,25,50), conf_level=0.90, seed=42L, out_dir=file.path(tempdir(),'default-validation'), save_plots=FALSE, save_tables=FALSE, advanced=NULL); old_opt <- options(ipdcstorm.wind_field_mode=NULL, ipdcstorm.hindcast_sampler_mode=NULL); on.exit(options(old_opt), add=TRUE); res <- suppressWarnings(suppressMessages(ipdcstorm:::.run_hindcast_attribution_case(cfg=cfg, targets=targets, validation_cfg=validation_cfg, climate=make_climate_cfg(scenario='stationary'), model_seed=2026L, wind_field_mode='legacy', annual_rate_mode='raw', sampler_mode='legacy'))); cat('default_wind_mode=', getOption('ipdcstorm.wind_field_mode', 'legacy'), '\n', sep=''); cat('default_sampler_mode=', ipdcstorm:::.hindcast_sampler_mode(), '\n', sep=''); cat('default_raw_lambda=', validation_cfg[['advanced']][['hindcast_use_raw_rates']], '\n', sep=''); tbl <- res[['comparison']]; tbl <- tbl[tbl[['location']] %in% c('Miami','Puerto_Rico','Saba','St_Martin','Statia'), c('location','return_period','bias_pct','sim_xi','obs_xi','wind_field_mode','annual_rate_mode','sampler_mode')]; tbl[['bias_pct']] <- round(tbl[['bias_pct']], 2); tbl[['sim_xi']] <- round(tbl[['sim_xi']], 3); tbl[['obs_xi']] <- round(tbl[['obs_xi']], 3); tbl <- tbl[order(tbl[['location']], tbl[['return_period']]), ]; print(tbl)"
```

## 7. Test results

- Parse checks: passed.
- Targeted tests: passed with `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 167 ]`.
- `Rd` checks for touched man pages: passed.
- Workspace default-path validation rerun: passed and confirmed `default_wind_mode=legacy`, `default_sampler_mode=legacy`, `default_raw_lambda=TRUE`.

## 8. Behavior changes

User-facing behavior is now explicitly documented and tested as the legacy validated default path:
- raw lambda in Tier 1A hindcast
- legacy wind-field path
- legacy hindcast sampler

Experimental observed-R34 candidate modes remain reachable only through internal option overrides and do not affect default runs.

Default validation summary from the workspace rerun:

| Location | RP5 bias % | RP10 bias % | RP25 bias % | RP50 bias % | Sim xi | Obs xi |
|---|---:|---:|---:|---:|---:|---:|
| Miami | 17.44 | 8.57 | 6.23 | 5.87 | -0.219 | -0.300 |
| Puerto_Rico | -0.75 | 1.43 | 3.17 | 4.20 | -0.232 | -0.274 |
| Saba | 10.83 | 9.94 | 10.47 | 11.12 | -0.245 | -0.300 |
| St_Martin | -3.40 | -2.90 | -2.71 | -2.50 | -0.264 | -0.300 |
| Statia | 12.47 | 16.10 | 21.74 | 26.20 | -0.068 | -0.248 |

## 9. Follow-ups/risks

- Residual upper-tail bias at Saba and Statia remains unresolved and is intentionally deferred.
- Prior internal observed-R34 candidate modes remain useful for research diagnostics, but they should not be promoted without a deeper radial-wind formulation fix.
- The worktree still contains unrelated pre-existing changes outside this task, including pipeline/man edits and older report-file deletions, which were left untouched.
