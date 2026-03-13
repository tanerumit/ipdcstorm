## 1. Goal
Refactor the hazard climate workflow so `delta_sst` is the canonical climate driver, with equivalent scenario-helper and direct-input pathways, and no downstream dependence on scenario labels or future timing once `delta_sst` is fixed.

## 2. Scope
Included: `R/hazard_climate.R`, `R/hazard_run.R`, regenerated `man/` files for touched exports, and climate tests in `tests/testthat/test-climate.R`.
Excluded: historical recalibration logic, IBTrACS selection, wind/downscaling physics, damage logic, and broader validation framework changes beyond existing climate-facing tests.

## 3. Problem solved
The previous workflow still let scenario period metadata drive applied climate adjustments through response-regime logic and hidden future-period construction in `run_hazard_model()`. That broke DMDU-style equivalence because identical `delta_sst` values could produce different downstream behavior depending on provenance.

## 4. Summary
Added explicit climate input modes to `make_climate_cfg()`: `scenario_helper` and `direct_delta_sst`. Scenario helper mode now resolves `delta_sst` from `scenario` plus `target_year` or `future_period`; direct mode accepts `delta_sst` directly and rejects incompatible mixed specifications.

`resolve_climate_inputs()` now resolves one canonical `delta_sst` and applies fixed `delta_sst`-only guardrails, removing period-dependent response logic. `run_hazard_model()` no longer synthesizes a future period internally, records input mode and provenance as metadata, and hashes runs from resolved climate effects rather than scenario labels.

## 5. Files changed
`R/hazard_climate.R`: +227 / -80. Added target-year resolution helpers, direct `delta_sst` support, mixed-spec validation, `delta_sst`-only response guardrails, and updated roxygen.
`R/hazard_run.R`: +23 / -10. Removed hidden future-period construction, carried input mode/target metadata, added `climate_input_mode` fit attribute, and changed parameter hashing to resolved climate effects.
`tests/testthat/test-climate.R`: +111 / -20. Added direct-mode, helper/direct equivalence, conflict validation, and metadata assertions; updated deterministic reference expectations.
`man/get_scenario_delta.Rd`: +14 / -5. Regenerated for `target_year` support.
`man/make_climate_cfg.Rd`: +33 / -8. Regenerated for dual climate pathways and validation rules.
`man/resolve_climate_inputs.Rd`: +9 / -7. Regenerated for canonical `delta_sst` resolution and metadata fields.

## 6. Commands run
`git branch --show-current`
`git status --short`
`Get-Content -Path C:\Users\taner\WS\shared-tools\01_skills\r-coding\SKILL.md -TotalCount 250`
`Get-Content -Path AGENTS.md -TotalCount 250`
`rg -n "make_climate_cfg|resolve_rate_response_regime|delta_sst|future_period|scenario" R tests/testthat man`
`Rscript -e "parse(file='R/hazard_climate.R')"`
`Rscript -e "parse(file='R/hazard_run.R')"`
`Rscript -e "parse(file='tests/testthat/test-climate.R')"`
`Rscript -e "testthat::test_file('tests/testthat/test-climate.R')"` 
`Rscript -e "testthat::test_file('tests/testthat/test-hazard-validation-api.R')"` 
`Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-climate.R')"`
`Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-hazard-validation-api.R')"`
`Rscript -e "devtools::document()"`
`git status --short`
`git diff --stat`
`git diff --numstat`

## 7. Test results
Pass: `Rscript -e "parse(file='R/hazard_climate.R')"`
Pass: `Rscript -e "parse(file='R/hazard_run.R')"`
Pass: `Rscript -e "parse(file='tests/testthat/test-climate.R')"`
Pass: `Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-climate.R')"` with 157 passing assertions.
Pass: `Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-hazard-validation-api.R')"` with 19 passing assertions.
Pass: `Rscript -e "devtools::document()"`
Note: plain `testthat::test_file(...)` runs failed before `pkgload::load_all('.')` because the package functions were not loaded into the test session.

## 8. Behavior changes
Users can now supply climate by either scenario-helper metadata or direct `delta_sst`.
Equal resolved `delta_sst` values now produce identical downstream climate adjustments and run hashes regardless of scenario/year provenance.
Run metadata now records `input_mode`, `target_year`, `future_period`, and resolved `delta_sst`; scenario provenance remains metadata only.
Mixed direct/helper climate specifications now error when they disagree.

## 9. Follow-ups/risks
The fixed `delta_sst`-only guardrail replaces the old period-based response regime with a constant damping/bounds profile, matching the previous near-term guardrail behavior but intentionally changing late-period sensitivity behavior.
Existing deterministic future-reference tests had to be updated because the refactor now enforces provenance-independent downstream behavior instead of preserving older period-coupled snapshots.
