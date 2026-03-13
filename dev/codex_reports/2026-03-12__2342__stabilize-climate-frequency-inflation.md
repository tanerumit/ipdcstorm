## 1. Goal

Stabilize climate-driven storm frequency scaling so future hazard runs no longer produce implausible order-of-magnitude count inflation, while preserving the stationary baseline path and exposing reproducible climate diagnostics.

## 2. Scope

Included:
- Basin-consistent climate-count calibration changes in `R/hazard_climate.R` and `R/hazard_run.R`.
- Climate diagnostics/metadata additions, plausibility guardrails, and regression tests.
- Roxygen/man refresh for touched exports and a brief technical note under `dev/notes/`.

Excluded:
- Wind-field physics, R34 handling, event intensity downscaling, and damage/impact modules.
- Any redesign of the public hazard API beyond added metadata fields.

## 3. Problem solved

`beta_0` calibration was sensitive to the user target set because climate-rate fitting was fed target-stacked annual counts. That could duplicate basin storms across locations, overstate annual activity, and inflate `exp(beta_sst * delta_sst)` in future runs.

## 4. Summary

`run_hazard_model()` now calibrates climate sensitivity from basin-consistent annual counts derived from de-duplicated storm-year records, not from the location-stacked annual-count table. `estimate_beta_sst()` now rejects target-conditioned count inputs, returns the exact annual total series used for calibration, and applies stronger shrinkage plus a `beta_0` plausibility fallback. `resolve_climate_inputs()` now carries explicit diagnostics including `sst_scale` and caps the SST count multiplier at `4x`.

Observed post-change check on a fixed seed (`2026`) for `scenario = "ssp585"` with the package test target:
- Stationary mean annual total: `1.05`
- SSP5-8.5 mean annual total: `1.73`
- Future/baseline ratio: `1.66`
- `delta_sst = 0.58 degC`, `beta_sst = 1.01 1/degC`, `sst_scale = 1.80`

## 5. Files changed

- `R/hazard_climate.R` — `+144/-67`; added basin-count provenance, `beta_0` guardrails, `sst_scale`, and richer resolved diagnostics.
- `R/hazard_run.R` — `+152/-138`; switched climate calibration to basin de-duplicated counts and exposed climate diagnostics in outputs/run metadata.
- `tests/testthat/test-climate.R` — `+317/-111`; added regression coverage for diagnostics, target-duplication independence, input rejection, and multiplier sanity; updated expected future reference values.
- `man/estimate_beta_sst.Rd` — `+10/-1`; documented count provenance and guardrail outputs.
- `man/make_climate_cfg.Rd` — `+3/-7`; refreshed generated docs to match current source.
- `man/make_hazard_cfg.Rd` — `+5/-0`; refreshed generated docs to match current source.
- `man/resolve_climate_inputs.Rd` — `+10/-2`; documented `sst_scale`, annual-count provenance, and guardrails.
- `man/run_hazard_model.Rd` — `+32/-46`; documented basin-consistent climate calibration and run metadata climate diagnostics.
- `dev/notes/2026-03-12_climate-frequency-guardrails.md` — new 15-line technical note.

Note: the worktree already contained unrelated modifications outside this task scope; they were left in place.

## 6. Commands run

- `Get-Content AGENTS.md -Raw`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-coding\SKILL.md -Raw`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\task-workflow\SKILL.md -Raw`
- `git branch --show-current`
- `rg -n "estimate_beta_sst|beta_sst|delta_sst|sst_scale|gamma|p_hur_base|future_period|climate" R tests/testthat`
- `Get-Content R\hazard_climate.R -Raw`
- `Get-Content R\hazard_run.R -Raw`
- `rg -n "compute_annual_counts <-|compute_annual_counts\(" R`
- `Get-Content tests\testthat\test-climate.R -Raw`
- `Rscript -e "parse(file='R/hazard_climate.R'); parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='tests/testthat/test-climate.R')"`
- `Rscript -e "devtools::test(filter = 'climate|hazard_run|validation')"`
- `Rscript -e "devtools::document()"`
- `Rscript -e "<load_all and fixed-seed stationary/ssp585 diagnostic summary>"`

## 7. Test results

- Parse checks passed for `R/hazard_climate.R`, `R/hazard_run.R`, and `tests/testthat/test-climate.R`.
- Scoped tests passed: `Rscript -e "devtools::test(filter = 'climate|hazard_run|validation')"` → `PASS 271`, `FAIL 0`, `WARN 0`, `SKIP 0`.
- `devtools::document()` completed successfully and regenerated Rd files for touched exports.

## 8. Behavior changes

- Climate-rate calibration is now invariant to duplicated target geometry because it uses basin-consistent annual counts.
- Future runs expose `future_period`, `delta_sst`, `beta_0`, `beta_sst`, `gamma_0`, `gamma`, `p_hur_base`, `sst_scale`, and count-series provenance in `out$cfg$climate` and climate diagnostics in `out$run_metadata$climate`.
- Pathological positive `beta_0` estimates now fall back to the literature prior above `1.2 1/degC`, and resolved future SST count multipliers are capped at `4x`.
- Stationary runs still resolve through the same climate path and keep `delta_sst = 0`, `sst_scale = 1`.

## 9. Follow-ups/risks

- Direct callers of `estimate_beta_sst()` must now supply basin-consistent annual counts; location-stacked inputs fail fast instead of being silently aggregated.
- The `4x` multiplier ceiling is an explicit scientific guardrail, not a re-fit of scenario physics; if future validation shows this ceiling is too strict or too loose, it should be revisited with basin-level evidence.
- The repository has unrelated pre-existing changes in other files; if those are later documented or rebased, rerun `devtools::document()` and the scoped suite to keep generated docs synchronized.
