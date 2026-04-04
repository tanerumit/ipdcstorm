# 1. Goal

Refactor the hazard-model API so climate is always embedded in `hazard_cfg`, baseline runs use `make_climate_cfg(scenario = "stationary")`, and the legacy climate-off path is removed from constructors, runtime resolution, metadata, tests, and docs.

# 2. Scope

Included:

- `R/hazard_climate.R`, `R/hazard_run.R`, and direct validation callers in `R/hazard_validation.R`
- tests covering climate config, `run_hazard_model()`, smoke coverage, and validation wrappers
- regenerated `man/` files for touched exported functions
- a technical note in `dev/notes/unified-climate-workflow.md`

Not included:

- scientific redesign of the climate equations
- unrelated hazard-model refactors
- backwards-compatibility shims, fallback argument translation, or soft deprecations
- README/vignette edits, because the repo currently has no climate-facing examples there beyond placeholders

# 3. Problem solved

The package had two baseline representations:

- `climate = NULL`, which bypassed climate resolution and emitted `"off"` metadata
- `make_climate_cfg(scenario = "stationary")`, which ran the resolver with `delta_sst = 0`

That split baseline behavior, metadata, and configuration state across two incompatible paths. The refactor removes the bypass so one `hazard_cfg` fully defines one hazard run and stationary/future scenarios share the same resolver.

# 4. Summary

`make_hazard_cfg()` now owns an embedded `climate` field and defaults it to `make_climate_cfg(scenario = "stationary")`. `run_hazard_model()` no longer accepts a top-level `climate` argument and fails fast if `cfg$climate` is missing, not a `climate_cfg`, or still contains legacy `enabled` state.

`make_climate_cfg()` and `resolve_climate_inputs()` no longer support climate-off mode or emit `"off"` metadata. Baseline and future runs now always pass through the same climate-resolution pipeline and return the same metadata schema.

# 5. Files changed

- `R/hazard_climate.R`: removed `enabled` support, removed `"off"` resolution branch, updated printing/docs for stationary baseline semantics
- `R/hazard_run.R`: embedded climate in `hazard_cfg`, removed `run_hazard_model(climate = ...)`, required `cfg$climate`, removed off-path metadata/output handling
- `R/hazard_validation.R`: updated direct callers/examples that forwarded the removed `climate` argument
- `tests/testthat/test-climate.R`: rewrote climate API coverage for embedded config, legacy rejection, and stationary/future reference outputs
- `tests/testthat/test-hazard-validation-api.R`: added hazard-config climate assertions
- `tests/testthat/test-hazard-validation-exported.R`: updated mocked `run_hazard_model()` signatures and expectations
- `tests/testthat/test-smoke.R`: aligned smoke coverage with embedded climate config and suppressed unrelated estimation warnings
- `man/make_climate_cfg.Rd`, `man/make_hazard_cfg.Rd`, `man/resolve_climate_inputs.Rd`, `man/run_hazard_model.Rd`, `man/validate_hazard_model.Rd`: regenerated from roxygen
- `dev/notes/unified-climate-workflow.md`: technical note on the API/behavior change

Line-count summary from `git diff --stat` for tracked files: 338 insertions, 352 deletions across 12 tracked files, plus 1 new note file.

# 6. Commands run

- `Get-Content AGENTS.md`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-coding\SKILL.md`
- `git branch --show-current`
- `rg -n "run_hazard_model|make_hazard_cfg|make_climate_cfg|climate\\s*=\\s*NULL|enabled\\s*=\\s*FALSE|climate_mode|off" R tests/testthat README* vignettes dev man`
- `Rscript -e "parse(file='R/hazard_climate.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "devtools::document()"`
- `Rscript -e "devtools::test(filter = 'climate|hazard_cfg|run_hazard_model')"`
- `Rscript -e "devtools::test(filter = 'hazard-validation|smoke')"`
- ad hoc `Rscript`/`pkgload::load_all('.')` comparison scripts to capture pre-refactor reference outputs and compare post-refactor stationary/future runs against them

# 7. Test results

- Parse checks: passed for `R/hazard_climate.R`, `R/hazard_run.R`, and `R/hazard_validation.R`
- Roxygen: `devtools::document()` passed and regenerated the touched Rd files
- Scoped tests: `devtools::test(filter = 'climate|hazard_cfg|run_hazard_model')` passed with `110` passes, `0` warnings, `0` failures
- Additional targeted tests: `devtools::test(filter = 'hazard-validation|smoke')` passed with `178` passes, `0` warnings, `0` failures

# 8. Behavior changes

- `run_hazard_model()` no longer accepts a separate `climate` argument
- `hazard_cfg` now contains a required `climate` field
- `make_hazard_cfg()` defaults baseline runs to embedded `make_climate_cfg(scenario = "stationary")`
- legacy `enabled = FALSE` / climate-off inputs now error immediately
- returned climate metadata is always present and no longer uses `"off"`

Results delta:

- Relative to legacy `climate = NULL`, the checked baseline case kept the same simulated summaries and rate table for the fixed seed comparison used here
- The change is in metadata and config traceability: legacy output reported `climate_mode = "off"` / `climate_scenario = "off"` / `climate_source = "off"`, while the unified baseline now reports `"baseline"` / `"stationary"` / `"builtin"`
- This is expected and acceptable because the removed off-path no longer bypasses climate resolution; baseline and future now share the same resolver and metadata structure

# 9. Follow-ups/risks

- `validate_hazard_model()` and hindcast helpers were updated because they directly called `run_hazard_model()`; this was necessary spillover from the API change
- The smoke/climate workflows still trigger `MASS::glm.nb()` convergence warnings in some unsuppressed interactive runs; tests now suppress those where they are not the target of the assertion
- If downstream scripts outside this repo still call `run_hazard_model(..., climate = ...)` or build `climate_cfg` objects with `enabled`, they will now fail immediately and need to be migrated
