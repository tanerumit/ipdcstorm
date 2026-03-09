# Validation Output Cleanup Report

## Goal
Remove validation output-level branching and backward-compatibility code so the suite emits only the retained reduced validation artifact contract.

## Scope
- Refactor `R/hazard_validation.R`
- Regenerate affected Rd files
- No new dependencies

## Summary
- Removed `cfg$output$level` and the `output` field from `make_validation_cfg()`.
- Removed the obsolete `advanced$n_boot` config knob from `make_validation_cfg()` because the retained workflow now always runs hindcast output with `n_boot = 0L`.
- Removed all `minimal`/`full` branching in `run_validation_suite()` and `plot_hindcast_validation()`.
- Removed obsolete internal diagnostics helpers for DIAGNOSTICS A-D and the unused Tier 1B modern blocked-CV helper stack.
- Simplified artifact writing to only:
  - `hindcast_return_levels.png`
  - `bias_decomposition.png`
  - `hindcast_return_levels.csv`
  - `rate_check.csv`
  - `wind_field.csv`
  - `validation_summary.csv`
- Removed markdown table bundle emission and all extra suite-generated plot artifacts tied to the old output-level split.
- Updated roxygen/Rd text to reflect the retained compact summary output and supported `advanced` fields.

## Files Changed
- `R/hazard_validation.R`
- `man/make_validation_cfg.Rd`
- `man/plot_hindcast_validation.Rd`
- `man/validate_hazard_model.Rd`

## Commands Run
1. `Get-Content .agents/skills/r-coding/SKILL.md`
2. `Get-Content DESCRIPTION`
3. `Get-Content NAMESPACE`
4. `rg -n "output\\$level|minimal|full|legacy|DIAGNOSTICS A|DIAGNOSTICS B|DIAGNOSTICS C|DIAGNOSTICS D|BLOCKED CV VALIDATION|RATE SANITY CHECK|run_validation_suite|hazard_validation" R tests`
5. Multiple line-numbered inspections of `R/hazard_validation.R`
6. `Get-Content tools/codex_reports/2026-03-09__1742__hazard-validation-artifact-refactor.md`
7. `Get-Date -Format "yyyy-MM-dd__HHmm"`
8. `Rscript -e "parse(file='R/hazard_validation.R')"` (failed once due BOM introduced during intermediate file rewrite)
9. BOM cleanup rewrite for `R/hazard_validation.R`
10. `Rscript -e "devtools::document()"`
11. `Rscript -e "parse(file='R/hazard_validation.R')"`
12. `Rscript -e "source('R/hazard_validation.R'); cfg <- make_validation_cfg(); stopifnot(is.null(cfg[['output']])); stopifnot(!('n_boot' %in% names(cfg[['advanced']]))); stopifnot(identical(sort(names(cfg[['advanced']])), sort(c('base_size','hindcast_use_raw_rates','xi_bounds'))));"`
13. `rg -n "output_level|cfg\\$output|minimal|full mode|DIAGNOSTICS A|BLOCKED CV VALIDATION|RATE SANITY CHECK|validation_tables\\.md|modern_blocked_cv|hindcast_dist_" R man NAMESPACE`
14. `git diff -- R/hazard_validation.R man/make_validation_cfg.Rd man/validate_hazard_model.Rd man/plot_hindcast_validation.Rd NAMESPACE`

## Test Results
- `parse(file='R/hazard_validation.R')`: PASS after BOM cleanup
- `devtools::document()`: PASS
- Config smoke check for removed `output` and `advanced$n_boot`: PASS
- No dedicated validation-suite integration test was run; this repository does not currently include a targeted `testthat` file for `run_validation_suite()`.

## Behavior Changes
- Removed config fields:
  - `cfg$output`
  - `cfg$output$level`
  - `advanced$n_boot` in `make_validation_cfg()`
- Removed console sections:
  - `DIAGNOSTICS A-D`
  - `BLOCKED CV VALIDATION (annual maxima, TS+)`
  - `RATE SANITY CHECK`
- Removed suite-generated artifacts:
  - `validation_tables.md`
  - per-island `hindcast_dist_*.png`
  - suite wiring for `rate_comparison.png`
  - suite wiring for `wind_field_scatter.png`
  - suite wiring for `qq_annual_max.png`
  - suite wiring for `cdf_comparison.png`
- Removed returned legacy fields from `run_validation_suite()`:
  - `output_level`
  - `bias_diagnostics`

## Follow-ups / Risks
- Exported standalone plotting functions for rate, wind-field, QQ, and CDF views still exist for direct use, but `run_validation_suite()` no longer emits those artifacts.
- End-to-end artifact generation with real hazard-model output was not rerun in this task, so retained file creation was verified by code-path inspection plus parse/config smoke checks rather than a full data-backed suite execution.
