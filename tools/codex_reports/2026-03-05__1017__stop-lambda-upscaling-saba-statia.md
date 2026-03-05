# Task Report

## Goal
Implement a down-only lambda scaling option to prevent frequency upscaling bias (especially Saba/Statia tails), keep backward-compatible default behavior, and propagate mode + run metadata into validation outputs.

## Scope
- Added lambda scaling mode support and guards in internal scaler helpers.
- Threaded mode through hazard run + validation pipelines.
- Extended `rate_check.csv`/summary content with scaling diagnostics and run metadata fields.
- Added/updated lambda scaling tests.

## Summary
- Added internal mode normalization: `target` (existing behavior) and `down_only` (prevents upscaling).
- Updated `.lambda_scalers_from_rate_check()`:
  - In `down_only`, caps raw scale at `<= 1.0` before existing clamps.
  - Added assertions for non-negative finite lambdas.
  - Added outputs: `lambda_scaling_mode`, `lambda_scale_applied`, `was_upscaled`.
  - Added explicit down-only assertion (`was_upscaled` must remain `FALSE`).
- Added one-line warning in `run_hazard_model()` when `target` mode causes any upscaling.
- Added config plumbing via `cfg$advanced$lambda_scaling_mode` with default `"target"`.
- Updated validation flow to use selected mode in hindcast training scaler application and rate-check generation.
- Extended rate-check table to include:
  - `lambda_scaling_mode`
  - `lambda_scale_applied`
  - `was_upscaled`
- Added validation assertion: in `down_only` mode, any upscaled row triggers error.
- Added run metadata propagation to validation outputs:
  - `seed`
  - `ibtracs_data_id` (filename + row count when available)
  - `parameter_id` hash
  - `lambda_scaling_mode`
- Updated wrapper seeding in `validate_hazard_model()` to set seed at entry and attach to `out$run_metadata`.

## Files Changed
- `R/hazard_utils.R`
- `R/hazard_run.R`
- `R/hazard_validation.R`
- `tests/testthat/test-lambda_scaling.R`

## Commands Run
- `Rscript -e "parse(file='R/hazard_utils.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "devtools::test(filter = 'lambda_scaling')"`

## Test Results
- Parse checks: all passed for touched R files.
- `devtools::test(filter = 'lambda_scaling')`: PASS (19), WARN (0), FAIL (0), SKIP (0).

## Behavior Changes
- Default remains `target` (backward-compatible).
- New optional mode `down_only` prevents lambda upscaling and preserves existing downscale/clamp behavior.
- Rate-check outputs now expose explicit scaling diagnostics and scaling mode.
- Validation summary and rate-check outputs now include run metadata fields listed above.

## Follow-ups / Risks
- No full hindcast rerun was executed in this patch step; acceptance thresholds for Saba/Statia/St_Martin still require the requested 2,000-year baseline-vs-patched validation runs.
- Repository had unrelated pre-existing working-tree changes; this patch avoided touching those files.
