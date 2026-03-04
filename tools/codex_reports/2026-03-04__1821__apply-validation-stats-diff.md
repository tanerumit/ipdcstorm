## Goal

Apply the user-provided unified diff to `dev/validation_stats.R` exactly, then run the full package test suite with `R.exe -q -e "devtools::test()"`.

## Scope

Only `dev/validation_stats.R` was edited for code changes. No exported R files, generated docs, or tests were modified.

## Summary

Applied the requested changes to add `min_year` and `storm_vmax_min` parameters to `.compute_site_annual_max()` and `run_validation_stats()`, inserted the requested filtering logic before annual-max aggregation, and passed the new arguments through the call site in `run_validation_stats()`.

## Files Changed

- `dev/validation_stats.R`
- `tools/codex_reports/2026-03-04__1821__apply-validation-stats-diff.md`

## Commands Run

- `Get-Content -Path .agents/skills/r-coding/SKILL.md -TotalCount 200`
- `Get-Content -Path dev/validation_stats.R`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `rg -n "run_validation_stats|compute_site_annual_max|validation_stats" tests/testthat R dev`
- `Get-ChildItem -Path tools -Force`
- `Rscript -e "parse(file='dev/validation_stats.R')"`
- `R.exe -q -e "devtools::test()"`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `git status --short`
- `git diff -- dev/validation_stats.R`

## Test Results

- `Rscript -e "parse(file='dev/validation_stats.R')"`: passed
- `R.exe -q -e "devtools::test()"`: passed, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 31 ]`

## Behavior Changes

- `run_validation_stats()` now accepts `min_year` (default `1970`) and `storm_vmax_min` (default `34`).
- Annual maxima are now computed after filtering to years `>= min_year` when `min_year` is finite.
- Annual maxima are now computed after filtering to storms with `Vmax_kt >= storm_vmax_min` when `storm_vmax_min` is finite.

## Follow-ups/Risks

- `dev/validation_stats.R` is currently an untracked dev helper, so this change is not covered by package exports or dedicated tests.
- The repo already had unrelated modified and untracked files; they were left untouched.
