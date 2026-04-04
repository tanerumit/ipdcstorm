## Goal
Fix missing arguments and API mismatches in `inst/extcode/01-baseline-validation.R` so the baseline validation pipeline runs against the current package API.

## Scope
Included:
- Minimal edits in the pipeline script only.
- Runtime verification by executing the pipeline script end to end.

Excluded:
- Changes to package internals in `R/`.
- Refactors to other pipeline or vignette scripts.

## Skills loaded
- always : r-coding
- loaded : task-workflow
- skipped: r-vignettes

## Problem solved
The pipeline had drifted from the current API in two places:
- It defined targets with `location` while downstream daily-generation usage expected the current `name` schema used by the pipeline object itself.
- Its stationary climate call relied on a helper path that now reaches a `run_hazard_model()` branch expecting a finite `target_year`, causing the script to fail before validation.

## Summary
Updated the target table to use `name`, aligned the daily generation call with that column, and changed the baseline climate configuration to an explicit zero-delta stationary configuration with a finite `target_year`. After these changes, `inst/extcode/01-baseline-validation.R` executed successfully through hazard generation, validation, and daily hazard-impact generation.

## Files changed
- `inst/extcode/01-baseline-validation.R` — 7 insertions, 3 deletions; fixed target schema and baseline climate arguments.
- `dev/codex_reports/2026-03-15__1627__fix-baseline-validation-pipeline.md` — new completion report.

## Commands run
- `Get-Content -Path AGENTS.md`
- `git branch --list`
- `git checkout agent`
- `Get-ChildItem -Force -Path .agents`
- `Get-ChildItem -Force -Path .agents\skills`
- `Get-ChildItem -Recurse -Force -Path .agents\skills | Select-Object FullName`
- Printed skill activation summary to console
- `rg --files | rg "01-baseline-validation\.R$|baseline-validation|pipelines"`
- `Get-Content -Path inst\extcode\01-baseline-validation.R`
- `rg -n "validate_hazard_model\s*<-|run_hazard_model\s*<-|make_climate_cfg\s*<-|make_hazard_cfg\s*<-|make_validation_cfg\s*<-|run_validation_suite\s*<-|generate_daily_hazard_impact\s*<-" R`
- `Rscript -e "source('inst/extcode/01-baseline-validation.R')"`
- `git diff --stat -- inst\extcode\01-baseline-validation.R`
- `git diff -- inst\extcode\01-baseline-validation.R`

## Test results
- Pass: `Rscript -e "source('inst/extcode/01-baseline-validation.R')"` completed successfully.
- Coverage: end-to-end runtime validation of the edited pipeline script only; no unit tests were added or modified.

## Behavior changes
- The pipeline now uses the current target naming convention directly (`name`).
- The baseline climate run is constructed via explicit zero warming plus finite target year, avoiding the current stationary-helper runtime failure path.

## Follow-ups/risks
- `run_hazard_model()` still appears to assume `climate_info$target_year` is length-1 numeric even for helper-based stationary runs; that package-level issue remains outside this task’s scope.
- The pipeline runtime is substantial (about 2.5 minutes in this environment), so repeated validation is relatively expensive.
