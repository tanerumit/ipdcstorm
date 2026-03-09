Goal
- Refactor hazard/validation simulation-count inheritance and terminology so validation can inherit `n_sim`, configs prefer `n_sim`, storm-class arguments prefer `storm_class` / `storm_classes`, and live code/output paths no longer use `TS34plus`.

Scope
- Updated hazard config/run, annual-count helpers, validation config/run, lambda-scaler utilities, regenerated affected Rd files, and added targeted tests for inheritance, overrides, alias compatibility, and `TS34plus` removal.
- Left unrelated dirty worktree changes untouched.

Summary
- `make_hazard_cfg()` now prefers `n_sim`, keeps legacy `n_sim_years` in sync for compatibility, and `run_hazard_model()` returns both `cfg` and `config`.
- `make_validation_cfg()` now defaults `n_sim = NULL`; `run_validation_suite()` resolves effective `n_sim` from `validation_cfg$n_sim` first, then `out$config$n_sim`, then `out$config$n_sim_years`, with fallback to legacy `out$cfg` and a hard stop if unresolved or `< 100`.
- Exported/internal counting and wrapper interfaces now prefer `storm_classes`, still accept deprecated `severities`, and normalize `TS34plus` input aliases to `TS`.
- `TS34plus` was previously acting as an aggregate (`TS + HUR`) in the rate-check/lambda-scaling pathway. That pathway was rewritten so exposed rows are `TS` and `HUR` only, while the former total-rate semantics are preserved internally by carrying total-rate context on the `TS` scaler row and deriving exclusive-TS targets from total-minus-hurricane reference targets.

Files changed
- `R/hazard_core.R`
- `R/hazard_run.R`
- `R/hazard_utils.R`
- `R/hazard_validation.R`
- `tests/testthat/test-hazard-validation-api.R`
- `tests/testthat/test-lambda_scaling.R`
- `man/compute_annual_counts.Rd`
- `man/dot-resolve_hazard_n_sim.Rd`
- `man/dot-resolve_validation_n_sim.Rd`
- `man/dot-validate_hindcast.Rd`
- `man/make_hazard_cfg.Rd`
- `man/make_validation_cfg.Rd`
- `man/run_hazard_model.Rd`
- `man/validate_hazard_model.Rd`

Commands run
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='R/hazard_utils.R')"`
- `Rscript -e "source('R/hazard_core.R'); source('R/hazard_utils.R'); source('R/hazard_run.R'); source('R/hazard_validation.R')"`
- `Rscript -e "devtools::document()"`
- `Rscript -e "devtools::test(filter = 'hazard|validation|core')"`

Test results
- Parse checks passed for all touched R files.
- Source/load command passed for the touched hazard files.
- `devtools::document()` completed successfully after fixing one roxygen block boundary.
- `devtools::test(filter = 'hazard|validation|core')` passed: 20 tests, 0 failures, 0 warnings, 0 skips.

Behavior changes
- Validation now inherits simulation length from model output by default instead of using a fixed validation default.
- `n_sim` is now the documented user-facing simulation-length field; legacy `n_sim_years` remains accepted and synchronized.
- Deprecated `severities` aliases still work but now warn and normalize immediately.
- Validation/reference/rate-check outputs now use `TS` and `HUR` only; no active output row uses `TS34plus`.

Follow-ups/risks
- `devtools::document()` also reflected unrelated pre-existing worktree changes outside this task; those files were not modified intentionally here.
- The rate-check `TS` row now represents the exclusive non-hurricane tropical-storm component derived from the former TS-or-higher reference totals, while total-rate calibration semantics are preserved internally for lambda adjustment.
