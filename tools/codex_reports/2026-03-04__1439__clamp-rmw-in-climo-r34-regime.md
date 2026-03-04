# Codex Report

## Goal

Apply a minimal internal patch to clamp the working RMW in the missing-radii / climatological-R34 regime inside `.estimate_site_wind_holland()`.

## Scope

- Touched only `R/hazard_core.R`.
- Added a climo-only `RMW_km0` guard immediately after climatological `R34` infill.
- Did not change any exported function signatures or return structures.

## Summary

- Added `rmw_over_r34_cap <- 3.0`.
- For rows with `ok`, `R34_is_climo`, and finite positive `R34_eff`, the working solver RMW is now clamped as:
  `RMW_km0 <- pmax(5, pmin(RMW_km0, R34_eff / rmw_over_r34_cap))`.
- This only affects the missing-radii regime where `R34` was infilled climatologically.

## Files Changed

- `R/hazard_core.R`

## Commands Run

- `Rscript -e "parse(file='R/hazard_core.R')"`
- `R.exe -q -e "devtools::test(filter = 'rmw|wind|holland|hazard')"`
- `R.exe -q -e "devtools::test()"`

## Test Results

- Parse check passed.
- `devtools::test(filter = 'rmw|wind|holland|hazard')` failed in `test-hazard_outer_cutoff.R`:
  - `fallback R34 decays immediately beyond the 1.8x cutoff`
  - actual `19.5` vs expected `32.5`
- `devtools::test()` failed on the same test and otherwise passed the remaining suites.

## Behavior Changes

- Only the internal climo-`R34` fallback path uses the new RMW clamp.

## Follow-Ups / Risks

- The new clamp changes fallback-path winds before the outer-cutoff decay is applied, so tests that assume the old fallback-path baseline wind now need to be updated or re-evaluated.
- I did not run the full manual validation suite after seeing the regression test failure, so no validated return-level delta is recorded in this patch report.
