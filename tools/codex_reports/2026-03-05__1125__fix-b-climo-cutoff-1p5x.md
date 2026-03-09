# Goal
Apply Fix B by tightening only the climatological R34 outer cutoff multiplier from `1.8x` to `1.5x`, keep observed pathway behavior unchanged, and expose effective multipliers in diagnostics.

# Scope
- Updated Holland outer-cutoff multiplier handling in climo pathway only.
- Added explicit diagnostic reporting for observed vs climo cutoff multipliers.
- Updated cutoff-focused tests to match the frozen spec.

# Summary
- Introduced `.holland_outer_cutoff_multipliers()` and set:
  - observed multiplier = `1.5`
  - climo multiplier = `1.5`
- `.resolve_holland_outer_cutoff_km()` now uses this config and enforces finite, non-negative cutoff radii.
- Validation diagnostics now carry and print both effective multipliers (`mult_obs`, `mult_climo`) in section B output, including skipped-path messages.
- Adjusted tests in `test-hazard_outer_cutoff.R` so fallback/climo expectations are `1.5x`.

# Files Changed
- `R/hazard_core.R`
- `R/hazard_validation.R`
- `tests/testthat/test-hazard_outer_cutoff.R`

# Commands Run
- `rg -n "1\\.8|R34_climo|r_cut|cutoff|outer" R tests`
- `rg -n "climo|observed radii|radii" R tests`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "testthat::test_file('tests/testthat/test-hazard_outer_cutoff.R')"`
- `Rscript -e "devtools::test(filter = 'validation|tier1|hindcast|cutoff')"`
- `Rscript -e "devtools::test()"`

# Test Results
- Parse checks:
  - `R/hazard_core.R`: pass
  - `R/hazard_validation.R`: pass
- Filtered tests (`validation|tier1|hindcast|cutoff`): pass (`PASS 22`, `FAIL 0`, `WARN 0`, `SKIP 0`)
- Cutoff test file (`test-hazard_outer_cutoff.R`): pass (`PASS 6`, `FAIL 0`, `WARN 0`, `SKIP 0`)
- Full test suite: pass (`PASS 51`, `FAIL 0`, `WARN 0`, `SKIP 0`)

# Behavior Changes
- Climo/fallback outer cutoff multiplier changed from `1.8x` to `1.5x`.
- Observed pathway multiplier remains `1.5x`.
- Validation diagnostics now explicitly report both effective cutoff multipliers.
- No changes made to lambda/frequency scaling, RMW logic, gamma steepening, blend weights, KDE pooling, or GEV fitting.

# Follow-ups / Risks
- Acceptance criteria tied to Tier 2/Tier 3 site bias outcomes require running the full validation workflow with target datasets; test suite alone does not quantify Statia/St_Martin RP bias deltas.
