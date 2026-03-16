## Goal
Fix the direct-`delta_sst` runtime error exposed by `inst/extcode/02-climate-change.R`:
`Error in if (is.finite(climate_info$target_year)) { : argument is of length zero`

## Scope
Included: narrow inspection of the script and the invoked verbose climate reporting path, a minimal guard fix, parse validation, and targeted runtime validation for direct-`delta_sst` and scenario/year runs.
Excluded: refactors, simulation logic changes, climate calibration changes, and edits to unrelated files.

## Skills loaded
- always : r-coding
- loaded : (none)
- skipped: design-doc-mermaid, flowchart-creator, obsidian-skills, pptx, r-vignettes, task-workflow

## Problem solved
The verbose climate output path assumed `climate_info$target_year` was always a single finite numeric value. Direct-`delta_sst` future runs leave that field absent/empty, so `is.finite()` was called on a length-0 value and aborted execution.

## Summary
Inspected `inst/extcode/02-climate-change.R` first, then traced the failing guard to the verbose climate reporting block in `R/hazard_run.R`. Replaced the raw `is.finite(climate_info$target_year)` check with a scalar-safe guard that prints the target year only when it exists as a single finite numeric value.

## Files changed
- `R/hazard_run.R` — 5 lines added, 1 line replaced; defensive guard in verbose climate reporting.
- `dev/codex_reports/2026-03-16__1147__fix-target-year-guard.md` — new task report.

## Commands run
- `rg -n "climate_info|target_year|future_period|delta_sst|scenario" inst\extcode\02-climate-change.R`
- `$i=1; Get-Content inst\extcode\02-climate-change.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ }`
- `rg -n "target_year|climate_info" R inst tests`
- `$i=1; Get-Content R\hazard_run.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 710 -First 40`
- `$i=1; Get-Content tests\testthat\test-climate.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 1 -First 80`
- `Rscript -e "parse(file='inst/extcode/02-climate-change.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "pkgload::load_all('.'); cfg <- make_hazard_cfg(simulation_years = 2L, climate = make_climate_cfg(delta_sst = 1.5, sensitivity_mode = 'fixed')); targets <- tibble::tibble(name = 'Saba', lat = 17.63, lon = -63.23); out <- run_hazard_model(cfg = cfg, targets = targets, seed = 123L, verbose = TRUE); cat('DIRECT_DELTA_OK\n')"`
- `Rscript -e "pkgload::load_all('.'); cfg <- make_hazard_cfg(simulation_years = 2L, climate = make_climate_cfg(scenario = 'ssp245', target_year = 2050, sensitivity_mode = 'fixed')); targets <- tibble::tibble(name = 'Saba', lat = 17.63, lon = -63.23); out <- run_hazard_model(cfg = cfg, targets = targets, seed = 123L, verbose = TRUE); cat('SCENARIO_YEAR_OK\n')"`
- `git diff -- R\hazard_run.R`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`

## Test results
- Parse check `inst/extcode/02-climate-change.R`: PASS.
- Parse check `R/hazard_run.R`: PASS.
- Direct `delta_sst` runtime reproduction with `verbose = TRUE`: PASS; no `target_year` length-0 error, completed and printed `DIRECT_DELTA_OK`.
- Scenario-helper runtime reproduction with `verbose = TRUE`: PASS; still printed `Target year       : 2050.0` and completed with `SCENARIO_YEAR_OK`.
- Residual warnings: `theta.ml(...): iteration limit reached` during model fitting in both runtime checks; pre-existing and unrelated to this guard fix.

## Behavior changes
Only verbose/reporting behavior changed: the `Target year` line is now skipped when `climate_info$target_year` is missing, empty, `NA`, or non-finite. Model calibration, simulation, and climate-response calculations are unchanged.

## Follow-ups/risks
`inst/extcode/02-climate-change.R` still appears to contain unrelated runtime issues later in the script, including references to `future_period` and `out_585`, which were outside the requested fix. Similar scalar guards may be worth auditing in other optional climate metadata output paths if more direct-delta workflows are added.
