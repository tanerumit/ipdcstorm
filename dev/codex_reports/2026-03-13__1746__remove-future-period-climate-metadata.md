## 1. Goal
Remove `future_period` from the climate workflow, make scenario-helper resolution target-year-only, and trim redundant resolved climate metadata without changing the numerical climate model beyond the required metadata/plumbing cleanup.

## 2. Scope
Included: `R/hazard_climate.R`, `R/hazard_run.R`, `R/hazard_downscale.R`, `R/hazard_validation.R`, affected roxygen outputs in `man/`, and existing in-scope tests.
Excluded: climate science retuning, new dependencies, new tests, non-climate refactors, and any pipeline script smoke checks.

## 3. Problem solved
The climate path carried redundant timing and metadata fields (`future_period`, midpoint inference, duplicate perturb/source/rate-scale storage), which made the DMDU scenario workflow ambiguous and forced downstream code to branch across aliases that no longer carried distinct meaning.

## 4. Summary
Scenario-helper climate resolution now requires `scenario + target_year` and never fabricates a period from `start_year`. `resolve_climate_inputs()` now returns one timing field (`target_year`), one perturb field (`perturb`), one source field (`source`), and a clearer split between raw and applied rate sensitivity/multipliers: `beta_sst_raw`, `beta_sst`, `raw_rate_scale`, and `rate_scale`.
`run_hazard_model()`, daily downscaling, and validation summaries were patched to consume the simplified schema, and the existing climate/validation-api tests were updated to assert the new target-year-only workflow and the renamed metadata.

## 5. Files changed
- `R/hazard_climate.R` (+33/-122): removed `future_period` API/plumbing, deleted midpoint/default-period inference, removed `cc_params` and `scenario_start_year` aliases, and renamed resolved raw/applied rate metadata.
- `R/hazard_run.R` (+15/-21): updated climate metadata assembly, fit attributes, parameter hashing, and run metadata to the simplified schema.
- `R/hazard_downscale.R` (+0/-4): removed fallback to the redundant `cc_params` fit attribute and used canonical `perturb` only.
- `R/hazard_validation.R` (+3/-3): switched validation summary fields to `rate_scale` and `raw_rate_scale`, plus updated the climate example.
- `tests/testthat/test-climate.R` (+32/-45): patched existing tests for target-year-only helpers and the renamed raw/applied climate fields.
- `tests/testthat/test-hazard-validation-api.R` (+1/-1): added the required `target_year` to the scenario-helper config.
- `man/get_scenario_delta.Rd` (+3/-13), `man/make_climate_cfg.Rd` (+2/-7), `man/resolve_climate_inputs.Rd` (+4/-8), `man/run_hazard_model.Rd` (+3/-3), `man/validate_hazard_model.Rd` (+1/-1): regenerated documentation.

## 6. Commands run
- `git branch --show-current`
- `rg -n "future_period|scenario_start_year|sst_source|rate_scale|applied_rate|climate_cfg|resolve_climate_inputs|get_scenario_delta|make_climate_cfg|resolve_climate_target|calibrate_climate_baseline" R tests/testthat`
- `Rscript -e "parse(file='R/hazard_climate.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='R/hazard_downscale.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='tests/testthat/test-climate.R')"`
- `Rscript -e "parse(file='tests/testthat/test-hazard-validation-api.R')"`
- `Rscript -e "devtools::document()"`
- `Rscript -e "devtools::test(filter='climate|hazard-validation-api')"`
- `git diff --numstat -- R/hazard_climate.R R/hazard_run.R R/hazard_downscale.R R/hazard_validation.R tests/testthat/test-climate.R tests/testthat/test-hazard-validation-api.R man/get_scenario_delta.Rd man/make_climate_cfg.Rd man/resolve_climate_inputs.Rd man/run_hazard_model.Rd man/validate_hazard_model.Rd`
- `git status --short`

## 7. Test results
- Parse checks: passed for all touched R files and patched tests.
- `Rscript -e "devtools::document()"`: passed.
- `Rscript -e "devtools::test(filter='climate|hazard-validation-api')"`: passed, `174` tests passed, `0` failed, `0` warnings, `0` skipped.

## 8. Behavior changes
- Removed public/internal argument: `future_period` from `make_climate_cfg()`, `get_scenario_delta()`, and `resolve_climate_inputs()`.
- Scenario-helper runs now require `target_year`; no midpoint inference or default fabricated future period remains.
- Removed resolved/runtime fields: `future_period`, `scenario_start_year`, `cc_params`, `f_rate_climate`, `raw_sst_scale`.
- Renamed resolved/runtime fields:
  `beta_sst_effective` -> `beta_sst`
  `beta_sst` (raw pre-damping meaning in resolver output) -> `beta_sst_raw`
  `sst_scale` -> `rate_scale`
  `raw_sst_scale` -> `raw_rate_scale`
- Downstream compatibility patches:
  `R/hazard_run.R` now stores `rate_scale` on fit attributes and in run metadata.
  `R/hazard_downscale.R` now reads only canonical `perturb`.
  `R/hazard_validation.R` now summarizes `rate_scale` and `raw_rate_scale`.

## 9. Follow-ups/risks
- The scoped suite is green, but there may be out-of-scope callers or notebooks outside `tests/testthat/` still passing `future_period` or reading removed metadata names.
- The deterministic future helper tests now pin explicit `target_year = 2050`; any external workflows that depended on the old implicit timing must now choose a year explicitly.
