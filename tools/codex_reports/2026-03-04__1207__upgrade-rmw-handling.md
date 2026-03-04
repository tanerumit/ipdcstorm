# Goal
Upgrade radius-of-maximum-wind handling used by wind downscaling so observed `USA_RMW` is used first, radii-based inference uses fixed calibrated coefficients, and fallback estimates are bounded.

# Scope
- Updated [R/hazard_ibtracs.R](C:/Users/taner/WS/ipdcstorm/R/hazard_ibtracs.R) to convert `USA_RMW` from nm to km and mark implausible values as missing.
- Updated [R/hazard_core.R](C:/Users/taner/WS/ipdcstorm/R/hazard_core.R) to add deterministic RMW resolution helpers, calibrated radii mapping, and guarded fallback logic inside `compute_site_winds_full()`.
- Added [tests/testthat/test-rmw_estimation.R](C:/Users/taner/WS/ipdcstorm/tests/testthat/test-rmw_estimation.R) for precedence, invalid-observation handling, coefficient stability, caps, and determinism.

# Summary
- Observed `rmw_km` now comes from `USA_RMW` converted to km, with valid range `(5, 150)` km.
- `compute_site_winds_full()` now resolves `RMW_km` in this order: valid observed `rmw_km`, calibrated storm-wide mean-radii mapping, then guarded Knaff fallback.
- Fixed coefficients were fit once from the packaged North Atlantic IBTrACS sample and stored in code comments:
  - `0.6517550 * mean(R64_km)`
  - `0.6676334 * mean(R50_km)` when `R64` is unavailable
  - `0.4106665 * mean(R34_km)` when `R64` and `R50` are unavailable
- All inferred RMW values now share the same bounds: minimum `8` km and an intensity-dependent upper cap `max(50, 140 - 0.75 * (Vmax_kt - 34))`.

# Files Changed
- [R/hazard_core.R](C:/Users/taner/WS/ipdcstorm/R/hazard_core.R)
- [R/hazard_ibtracs.R](C:/Users/taner/WS/ipdcstorm/R/hazard_ibtracs.R)
- [tests/testthat/test-rmw_estimation.R](C:/Users/taner/WS/ipdcstorm/tests/testthat/test-rmw_estimation.R)

# Commands Run
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `Rscript -e "parse(file='R/hazard_ibtracs.R')"`
- `Rscript -e "devtools::test(filter = 'rmw|RMW|downscale')"`

# Test Results
- Parse checks passed for both touched R files.
- Targeted test filter passed: `PASS 5`, `FAIL 0`, `WARN 0`, `SKIP 0`.

# Behavior Changes
- Valid observed `USA_RMW` now directly controls `RMW_km` in wind downscaling.
- Invalid observed `USA_RMW` values no longer leak through as zero/oversized radii; they fall back to calibrated inference.
- Old fixed directional multipliers (`0.35`, `0.40`, `0.50`) are no longer used for RMW inference.
- Knaff fallback can no longer return unbounded large RMW values for weak/low-latitude edge cases.

# Follow-ups/Risks
- The stored coefficients were fit from the packaged North Atlantic IBTrACS sample in this repo; if calibration should be restricted further (for example, Caribbean-only), the constants should be refit from that narrower subset.
- Event summarization still aggregates the original `rmw_km` track column; this patch changes the wind-model `RMW_km` used during site-wind computation, not downstream event summaries.
