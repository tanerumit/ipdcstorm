## 1. Goal
Revise the default climate-response logic so near-term future runs keep basinwide tropical-cyclone counts close to baseline, shift the climate signal toward hurricane fraction and local redistribution, and expose validation diagnostics that separate basin rate change from redistribution.

## 2. Scope
Included:
- `R/hazard_climate.R`
- `R/hazard_run.R`
- `R/hazard_validation.R`

Not included:
- test file updates
- unrelated wind-field, hindcast, or late-century retuning
- dependency changes

## 3. Problem solved
Default `ssp585` runs centered on 2035 were inflating local annual counts by roughly `1.7x` to `1.9x`, which is stronger than literature support for near-term North Atlantic / Caribbean TC projections. The patch constrains near-term basinwide count response and pushes the climate signal toward hurricane mix and modest redistribution instead.

## 4. Summary
Added a climate-response regime resolver that classifies future periods as `near_term`, `transition`, or `late_century`, damps the raw SST-driven count multiplier, and hard-bounds the applied basinwide rate multiplier. The run path now applies a separate mean-preserving location redistribution term, records regime and parameter hash metadata, and the validation suite now reports basin count change, hurricane-fraction change, and redistribution diagnostics.

## 5. Files changed
- `R/hazard_climate.R` — added regime-aware damping/bounds for future count scaling; added raw vs applied rate metadata; minimal localized edit.
- `R/hazard_run.R` — switched simulation to applied basin multiplier, added explicit location redistribution scaling, and expanded run metadata with regime/hash fields.
- `R/hazard_validation.R` — added climate-response diagnostics summarizing basin count ratio, hurricane-fraction shift, and redistribution metrics.

## 6. Commands run
- `Rscript -e "parse(file='R/hazard_climate.R')"`
- `Rscript -e "parse(file='R/hazard_run.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "devtools::test(filter = 'climate|hazard-validation-api|smoke')"`
- `Rscript -e "devtools::test(filter = 'climate')"`
- `Rscript -e "devtools::test(filter = 'hazard-validation-api|smoke')"`
- `Rscript -e "<load package and compare stationary vs ssp585 start_year=2035 outputs across Saba/St_Eustatius/St_Maarten>"`
- `git diff --unified=3 -- R/hazard_climate.R R/hazard_run.R R/hazard_validation.R`

## 7. Test results
- Parse checks: passed for all three touched files.
- `devtools::test(filter = 'hazard-validation-api|smoke')`: passed.
- `devtools::test(filter = 'climate')`: one expected failure remains because an existing snapshot-style assertion encodes the old inflated future-count behavior.
- Deterministic comparison run with fixed seed:
  - baseline basin count ratio: `0.995`
  - future near-term basin count ratio: `1.057`
  - future raw multiplier: `1.800`
  - future applied basin multiplier: `1.064`
  - baseline simulated hurricane fraction: `0.170`
  - future simulated hurricane fraction: `0.184`
  - hurricane-fraction change: `+0.958 pp`

## 8. Behavior changes
- Near-term future count response is now constrained by default instead of following the full raw SST scalar.
- Local count changes now come from a separate redistribution term rather than a uniform count inflation.
- Validation output now includes explicit climate-response diagnostics for basin count change, hurricane share, and redistribution.
- Run metadata now records climate regime, raw/applied rate multipliers, redistribution settings, and a parameter hash.

## 9. Follow-ups/risks
- Existing climate tests that assert old future-count totals should be updated if the new behavior is accepted.
- Redistribution is intentionally modest and tied to baseline hurricane fraction; if literature-backed spatial diagnostics become available, that term can be replaced with a more physical driver.
- I did not add new user-facing `run_hazard_model()` documentation because the public signature did not change, but interpretation text may merit a follow-up doc update.
