## Goal
Improve the documentation of `make_sst_cfg()` by adding a `Details` section that explains each advanced parameter.

## Scope
- Update roxygen docs for `make_sst_cfg()` in `R/hazard_climate.R`.
- Regenerate Rd documentation affected by that roxygen block.
- Run required validation commands.

## Summary
- Expanded `@details` for `make_sst_cfg()` to describe all expert knobs in `advanced`:
  - `beta_sst`
  - `beta_prior`
  - `gamma_intensity`
  - `gamma_prior`
  - `cc_params` (including `NULL`, `list()`, and named-list behavior)
- Regenerated docs with roxygen; `man/make_sst_cfg.Rd` now contains the expanded details text.

## Files Changed
- `R/hazard_climate.R`
- `man/make_sst_cfg.Rd`

## Commands Run
1. `Rscript -e "roxygen2::roxygenise()"`
2. `Rscript -e "parse(file='R/hazard_climate.R')"`
3. `Rscript -e "testthat::test_dir('tests/testthat')"`
4. `Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat')"`

## Test Results
- Parse check passed:
  - `Rscript -e "parse(file='R/hazard_climate.R')"`
- Initial test run failed because package functions were not loaded in that invocation:
  - Failure: `could not find function "make_hazard_cfg"`
- Re-run in package-loaded context passed:
  - `Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat')"`
  - Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ]`

## Behavior Changes
- No runtime behavior changes.
- Documentation-only change for `make_sst_cfg()` details.

## Follow-ups/Risks
- None identified for runtime behavior.
