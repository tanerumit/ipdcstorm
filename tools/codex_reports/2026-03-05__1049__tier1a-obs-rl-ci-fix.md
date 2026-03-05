## Goal
Fix Tier 1A hindcast RL-in-CI logic so observed return-level 90% CIs are computed when feasible, and CI-unavailable rows are not mislabeled as outside CI.

## Scope
- `R/hazard_validation.R`
- `tests/testthat/test-hindcast-rl-ci.R`

## Summary
- Added robust observed RL CI helper with method order:
  - Delta-method CI from a numerically estimated GEV parameter covariance matrix.
  - Deterministic parametric bootstrap fallback (hurdle-GEV simulation + refit).
  - CI unavailable when sample support is insufficient (`n_pos < 10` by default).
- Wired Tier 1A hindcast to compute observed CIs via the new helper even when full bootstrap output mode is disabled.
- Updated Tier 1A console reporting:
  - If CI exists: keep inside/outside classification.
  - If CI unavailable: print `CI unavailable (n_pos=...)` and do not print outside mark.
- Added focused testthat coverage for:
  - Moderate sample finite CI.
  - Small-to-moderate sample finite CI.
  - Very small sample unavailable CI + NA outside-flag behavior.
  - Deterministic bootstrap fallback path.

## Files Changed
- `R/hazard_validation.R`
- `tests/testthat/test-hindcast-rl-ci.R`

## Commands Run
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "devtools::test(filter = 'hindcast|tier1a|ci|gev')"`

## Test Results
- `devtools::test(filter = 'hindcast|tier1a|ci|gev')`: PASS (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`).

## Behavior Changes
- Tier 1A observed RL CI is now available in many cases where it was previously always `NA` in minimal mode.
- If observed CI is unavailable, `obs_in_90ci` remains `NA` and console output no longer reports `✗ outside 90% CI` for that row.
- Comparison table now includes `obs_ci_method` (`delta`, `bootstrap`, or `unavailable`) and hindcast per-island result now stores `obs_gev` for diagnostics.

## Follow-ups / Risks
- Delta-method covariance uses a numerical Hessian and can fail in ill-conditioned cases; bootstrap fallback covers these but adds runtime when invoked.
- Full package check (`devtools::check()`) was not run here due expected heavier runtime; run it if you want full CRAN-style validation.
