## Goal

Introduce deterministic multiplicative lambda scalers by site and validation class (`TS34plus`, `HUR`), and apply them to count simulation so occurrence rates align with the existing reference climatology plus gate approximation without changing wind physics.

## Scope

- Updated internal rate-scaler helpers and deterministic scaler metadata.
- Applied calibrated lambdas in `run_hazard_model()` simulation and hindcast validation count generation.
- Expanded rate-check reporting to show raw and adjusted lambdas plus scaler status.
- Added focused tests for scaler derivation and simulation usage.
- Added a developer note documenting the numerical behavior.

## Summary

- Added `.lambda_scalers_from_rate_check()`, `.apply_lambda_scalers_to_lambda_table()`, and `.lambda_scaler_id()` as internal helpers.
- Refactored rate-check construction into `.build_rate_check_table()` so the same reference logic is reused for runtime calibration and validation reporting.
- `run_hazard_model()` now computes site/class scalers once from the existing reference table, logs scaler metadata, applies adjusted lambdas before `simulate_twolevel_counts()`, and returns `out$lambda_scalers` plus `out$lambda_scaler_id`.
- Hindcast validation now calibrates the training-period lambda table before simulating annual maxima, so GEV validation uses the adjusted event-rate intensity as well.
- Rate validation now reports `lambda_model_raw`, `lambda_ref`, `lambda_scale`, `lambda_adj`, clamp state, and evaluates flags against adjusted lambdas.

## Files Changed

- `R/hazard_run.R`
- `R/hazard_utils.R`
- `R/hazard_validation.R`
- `tests/testthat/test-lambda_scaling.R`
- `dev/notes/lambda_rate_scaling.md`

## Commands Run

- `Get-Content -Path .agents/skills/r-coding/SKILL.md`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `Get-ChildItem -Path tests\testthat -File | Select-Object -ExpandProperty FullName`
- `Get-Content -Path R\hazard_run.R`
- `Get-Content -Path R\hazard_validation.R`
- `Get-Content -Path R\hazard_core.R`
- `Get-Content -Path R\hazard_utils.R`
- `Get-Content -Path tests\testthat\test-smoke.R`
- `rg -n "rate_check|validate_rates|simulate_twolevel_counts|lambda_model|lambda_ref|Rate sanity" R\hazard_validation.R R\hazard_run.R R\hazard_core.R`
- `rg -n "simulate_twolevel_counts <- function|simulate_twolevel_counts\(" R`
- `Get-Content -Path R\hazard_climate.R` (targeted line ranges)
- `Get-ChildItem -Path dev\notes -File | Select-Object -ExpandProperty Name`
- `Get-ChildItem -Path tools\codex_reports -Force -ErrorAction SilentlyContinue | Select-Object -First 1`
- `git status --short`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `Rscript -e "parse(file='R/hazard_utils.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "devtools::test(filter = 'rate|lambda|validation')"`
- `git diff --stat`

## Test Results

- `Rscript -e "parse(file='R/hazard_utils.R')"`: passed
- `Rscript -e "parse(file='R/hazard_validation.R')"`: passed
- `Rscript -e "parse(file='R/hazard_run.R')"`: passed
- `Rscript -e "devtools::test(filter = 'rate|lambda|validation')"`: passed, `15` tests, `0` failures, `0` warnings, `0` skips

## Behavior Changes

- Simulation event rates are now adjusted multiplicatively by site and validation class before annual counts are generated.
- Scalers target `lambda_ref * expected_ratio`, preserving the existing gate-approximation logic instead of forcing raw equality to the published climatology.
- Scalers are clamped to `[0.25, 4.0]`.
- Missing reference lambdas keep `lambda_scale = 1` with `scale_status = "no_ref"`.
- Severity generation per event is unchanged; only occurrence frequency changes.
- Validation output now reports raw lambda, reference lambda, scale, adjusted lambda, and scaler state, and the rate-check flag uses adjusted lambda.

## Follow-Ups / Risks

- I did not run the full hindcast validation suite, so the acceptance thresholds for hindcast CI coverage and rate-tier improvement still need an end-to-end run on real data.
- If many sites hit the upper clamp or if adjusted hurricane lambda would frequently exceed the adjusted total lambda, the clamp policy may need review before keeping the calibration enabled by default.
