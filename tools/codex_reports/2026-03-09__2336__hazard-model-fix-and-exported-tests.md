# Goal
Fix the `run_hazard_model()` lambda-scaling error caused by incomplete scaler tables, and add unit tests covering the exported functions in `hazard_core.R`, `hazard_downscale.R`, and `hazard_validation.R`.

# Scope
- Patch legacy/incomplete lambda-scaler handling in internal hazard utilities.
- Add regression coverage for the scaler failure mode.
- Add public API tests for exported functions in:
  - `R/hazard_core.R`
  - `R/hazard_downscale.R`
  - `R/hazard_validation.R`
- Verify with targeted and full test-suite runs.

# Summary
- Hardened `.apply_lambda_scalers_to_lambda_table()` so missing `lambda_target_total` and `lambda_model_total_raw` columns no longer trigger warnings or `if (...)` errors.
- Hardened `.build_rate_check_table()` so it tolerates lambda-scaler tables that do not include `lambda_scaling_mode`.
- Added three new test files that exercise the exported APIs in the requested hazard modules using deterministic lightweight fixtures.
- Added a regression test for applying lambda scalers when total-scaling columns are absent.

# Files Changed
- `R/hazard_utils.R`
- `R/hazard_validation.R`
- `tests/testthat/test-lambda_scaling.R`
- `tests/testthat/test-hazard-core-exported.R`
- `tests/testthat/test-hazard-downscale-exported.R`
- `tests/testthat/test-hazard-validation-exported.R`

# Commands Run
- `Get-Content .agents/skills/r-coding/SKILL.md -TotalCount 250`
- `rg -n "run_hazard_model|apply_lambda_scalers|lambda_target_total|lambda_model_total_raw|hazard_core|hazard_downscale|hazard_validation" R tests/testthat`
- `Get-ChildItem -Force`
- Multiple targeted `Get-Content` / numbered file inspections for `R/*.R` and `tests/testthat/*.R`
- `git status --short`
- `Rscript -e "testthat::test_dir('tests/testthat', reporter = 'summary', filter = 'hazard-core-exported|hazard-downscale-exported|hazard-validation-exported|lambda_scaling')"`
- `Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter = 'summary', filter = 'hazard-core-exported|hazard-downscale-exported|hazard-validation-exported|lambda_scaling', load_package = 'none')"`
- `Rscript -e "pkgload::load_all('.'); print(names(get_reference_rates()))"`
- `Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter = 'summary', load_package = 'none')"`

# Test Results
- Targeted hazard/scaler subset: passed.
- Full `tests/testthat` suite: passed.
- Residual warnings remain in `test-hindcast-rl-ci.R` about missing `sim_lo_90` / `sim_hi_90` columns. These were warning-only and did not fail the suite; they appear to predate this task.

# Behavior Changes
- `run_hazard_model()` and validation rate-check paths now accept lambda-scaler tables that omit total-scaling columns and/or `lambda_scaling_mode`, falling back to per-class `lambda_scale` where needed.
- Exported hazard APIs now have direct unit-test coverage for core geometry/rates, downscaling/damage, validation/config/plot wrappers, and validation convenience wrappers.

# Follow-ups/Risks
- The hindcast CI tests expose warning-only schema drift around `sim_lo_90` / `sim_hi_90`; that should be cleaned up separately if the intended CI column names changed.
- Plot tests currently verify artifact creation, not image contents. If plotting logic changes materially, snapshot or structural plot assertions may be worth adding later.
