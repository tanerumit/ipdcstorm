## 1. Goal
Stabilize and speed up the hazard-model unit tests around embedded climate config, fixed-seed reference behavior, hazard visualization output, Holland outer-cutoff expectations, and related smoke coverage without changing production behavior.

## 2. Scope
Included:
- Updating stale hazard/climate/smoke/rmw tests under `tests/testthat/`.
- Preserving and tightening fixed-seed reference assertions for stationary and future climate runs.
- Reducing repeated end-to-end hazard runs where equivalent resolver/simulation assertions were sufficient.
- Verifying the requested focused subset and the full suite.

Excluded:
- No production-code edits in `R/`.
- No scientific-method refactor or hazard/climate algorithm redesign.
- No new dependencies or helper files.

## 3. Problem solved
The current branch already contained most climate/embed API updates, but the test suite still had three issues:
- Stale tests targeting removed wind-field interfaces (`R34_source`, `RMW_provenance`) in `test-rmw_estimation.R`.
- Runtime-heavy climate tests doing redundant full `run_hazard_model()` passes for invariants that could be checked at resolver/simulation level.
- Smoke tests using more simulated years than needed for path-level coverage.

## 4. Summary
I updated the stale wind-field tests to assert the current fallback/climatology contract (`R34_missing`, `R34_is_climo`, `R34_eff_km`) instead of removed provenance fields and removed obsolete option-driven assertions tied to the retired `R34_source` argument. I also narrowed the climate tests to one stationary and one future frozen reference run, reused cached outputs more aggressively, replaced redundant end-to-end comparisons with direct `resolve_climate_inputs()` / `simulate_twolevel_counts()` checks where appropriate, and cut smoke fixtures to a single simulated year.

## 5. Files changed
- `tests/testthat/test-climate.R` — 549 lines total; reduced repeated full hazard runs, refreshed fixed-seed frozen expectations for 8-year reference runs, kept strict embedded-climate and negative-validation coverage.
- `tests/testthat/test-rmw_estimation.R` — 188 lines total; replaced stale removed-interface assertions with current fallback diagnostics and deterministic wind-field invariants.
- `tests/testthat/test-smoke.R` — 84 lines total; reduced smoke simulation size from 3 years to 1 year.
- `tests/testthat/test-hazard_viz.R` — already part of the branch’s updated baseline; retained deterministic path/filename assertions.
- `tests/testthat/test-hazard_outer_cutoff.R` — already part of the branch’s updated baseline; retained current `R34_is_fallback` interface assertions.

Git diff summary for touched test files:
- 5 files changed, 300 insertions, 219 deletions.

## 6. Commands run
Exact commands used during validation/debugging:
- `git branch --show-current`
- `git status --short`
- `Rscript -e "devtools::test(filter = 'hazard_viz|outer_cutoff|climate|hazard_cfg|run_hazard_model|smoke')"`
- `Rscript -e "devtools::test()"`
- `Rscript -e "devtools::test(filter = 'climate|rmw_estimation')"`
- `git diff --stat -- tests/testthat/test-climate.R tests/testthat/test-rmw_estimation.R tests/testthat/test-smoke.R tests/testthat/test-hazard_viz.R tests/testthat/test-hazard_outer_cutoff.R`

Exploratory local checks were also run to inspect current fixed-seed summaries and parameter hashes before freezing the updated expectations.

## 7. Test results
Pass status:
- Focused subset: `Rscript -e "devtools::test(filter = 'hazard_viz|outer_cutoff|climate|hazard_cfg|run_hazard_model|smoke')"` passed with 155 tests in 95.3 s.
- Full suite: `Rscript -e "devtools::test()"` passed with 417 tests in 105.6 s.

Runtime notes:
- Earlier focused-subset baseline in this task was 98.0 s.
- Final focused-subset runtime was 95.3 s.
- Climate file runtime dropped from about 70.4 s on the initial focused run to about 61.0 s after removing redundant end-to-end runs.

## 8. Behavior changes
User-facing/package behavior:
- None in production code.

Test behavior:
- Fixed-seed reference tests now freeze only the load-bearing outputs: summary counts, rate table, and run metadata fields.
- Wind-field tests now assert current fallback diagnostics rather than removed provenance columns.
- Smoke coverage still exercises the end-to-end packaged-data path, but with `simulation_years = 1L` to keep routine runs practical.

## 9. Follow-ups/risks
- The focused subset is still dominated by `test-climate.R`; additional reductions would likely require either more fixture-level reuse or splitting truly end-to-end reference coverage from resolver-only checks.
- The worktree already had substantial uncommitted climate/viz/test updates before this task; this report treats those as the baseline and documents only the incremental stabilization performed here.
- No “Results delta” note is needed because no production behavior changed.
