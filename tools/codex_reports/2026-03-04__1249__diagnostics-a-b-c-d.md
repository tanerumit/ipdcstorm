# Codex Report

## Goal

Implement validation-time diagnostics for the four Saba/Statia tail-bias hypotheses (A-D), print them in the console, and produce a findings note for Saba and Statia.

## Scope

- Added internal diagnostic helpers only.
- Kept exported function names and signatures unchanged.
- Limited code changes to wind estimation, validation diagnostics, and one targeted test.

## Summary

- Added an internal option gate (`ipdcstorm.disable_r34_calibration`) inside `.estimate_site_wind_holland()` so diagnostics can rerun the same site winds with R34 calibration disabled without changing default behavior.
- Added internal validation helpers in `R/hazard_validation.R` to:
  - rerun site winds with/without calibration,
  - compute A-D diagnostic summaries,
  - print a delimited `DIAGNOSTICS A–D` section with per-site subsections A, B, C, D.
- Hooked the diagnostic printout into `validate_hazard_model()` before `run_validation_suite()`.
- Added a focused regression test that confirms the internal option lowers calibrated winds in a known case.
- Ran the validation workflow once using the packaged IBTrACS file under `inst/extdata/...` and wrote the resulting Saba/Statia summary to `dev/notes/diagnostics_A_B_C_D.md`.

## Files Changed

- `R/hazard_core.R`
- `R/hazard_validation.R`
- `tests/testthat/test-rmw_estimation.R`
- `dev/notes/diagnostics_A_B_C_D.md`

## Commands Run

- `Rscript -e "parse(file='R/hazard_core.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='tests/testthat/test-rmw_estimation.R')"`
- `Rscript -e "testthat::test_file('tests/testthat/test-rmw_estimation.R')"`
- `Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-rmw_estimation.R')"`
- `Rscript -e "pkgload::load_all('.'); targets <- tibble::tribble(...); cfg <- make_hazard_cfg(); out <- run_hazard_model(...)"` (failed because `data/ibtracs/ibtracs.NA.list.v04r01.csv` is absent in this workspace)
- `@' ... '@ | Rscript -` with `cfg <- make_hazard_cfg(data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv")`, `.run_validation_diagnostics_ad(...)`, and `run_validation_suite(...)`
- `@' ... '@ | Rscript -` with `validate_hazard_model(cfg = cfg, targets = targets, validation_cfg = make_validation_cfg(save_plots = FALSE, save_tables = FALSE))` as an end-to-end smoke test of the public validation wrapper

## Test Results

- `parse(file=...)` succeeded for both touched R files and the touched test file.
- `testthat::test_file('tests/testthat/test-rmw_estimation.R')` failed when run without loading the local package namespace first; existing unqualified tests could not find exported functions, and the new option test hit the installed namespace instead of the working tree.
- `pkgload::load_all('.')` followed by `testthat::test_file('tests/testthat/test-rmw_estimation.R')` passed (`7` assertions / checks shown as passing in the reporter).
- Validation run succeeded when pointed at `inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv`.

## Behavior Changes

- Default wind-estimation behavior is unchanged unless the internal option `ipdcstorm.disable_r34_calibration` is set to `TRUE`.
- `validate_hazard_model()` now prints an extra `DIAGNOSTICS A–D` console section before the standard validation suite output.

## Follow-Ups / Risks

- The new diagnostic print hook lives in `validate_hazard_model()`. Calling `run_validation_suite()` directly on a precomputed `out` object will not emit the new diagnostics because target coordinates are not available there.
- Statia's C diagnostic can be inconclusive when top-tail rows lack finite `R34_km`; a follow-up could break out missing-`R34` top-tail rows explicitly if that becomes a recurring analysis need.
