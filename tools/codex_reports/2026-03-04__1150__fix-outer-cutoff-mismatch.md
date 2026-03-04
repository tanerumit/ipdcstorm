## Goal
Make the Holland-profile outer cutoff (`R_outer`) match the intended multipliers: `1.5 * R34` for observed radii and `1.8 * R34` for fallback/climatological radii, with regression tests to lock that behavior in.

## Scope
Limited to the internal Holland wind implementation in `R/hazard_core.R` and a new focused test file in `tests/testthat/test-hazard_outer_cutoff.R`.

## Summary
Added an internal helper, `.resolve_holland_outer_cutoff_km()`, to centralize the cutoff calculation and make the observed-vs-fallback rule explicit. Updated `.estimate_site_wind_holland()` to reuse the existing `R34_is_climo` flag so rows infilled by `estimate_R34_climo()` use the fallback multiplier while caller-supplied positive `R34_km` values use the observed multiplier. Added deterministic tests for both ratios and for the existing exponential-decay behavior immediately beyond each cutoff.

## Files Changed
- `R/hazard_core.R`
- `tests/testthat/test-hazard_outer_cutoff.R`

## Commands Run
- `Get-Content -Path .agents/skills/r-coding/SKILL.md`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `Get-ChildItem -Path R -Filter *.R | Select-Object -ExpandProperty FullName`
- `Get-ChildItem -Path tests\\testthat -Filter *.R | Select-Object -ExpandProperty FullName`
- `rg -n "estimate_site_wind_holland|R_outer|outer cutoff|R34" R tests\\testthat`
- `git status --short`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `Rscript -e "parse(file='tests/testthat/test-hazard_outer_cutoff.R')"`
- `R -q -e "devtools::test(filter = 'outer_cutoff|holland|wind')"` (failed in PowerShell because `R` is an alias there)
- `Rscript -e "devtools::test(filter = 'outer_cutoff|holland|wind')"`

## Test Results
- `Rscript -e "parse(file='R/hazard_core.R')"`: passed
- `Rscript -e "parse(file='tests/testthat/test-hazard_outer_cutoff.R')"`: passed
- `Rscript -e "devtools::test(filter = 'outer_cutoff|holland|wind')"`: passed, `FAIL 0 | WARN 0 | SKIP 0 | PASS 7`

## Behavior Changes
- Observed, caller-supplied positive `R34_km` now maps to `R_outer = 1.5 * R34_km`.
- Fallback/climatological `R34_km` generated inside `.estimate_site_wind_holland()` now maps to `R_outer = 1.8 * R34_km`.
- The post-cutoff behavior is unchanged: winds still use the existing exponential decay once `r_km > R_outer`.

## Follow-ups/Risks
- The observed-vs-fallback distinction is currently local to `.estimate_site_wind_holland()`: it treats only internally infilled radii as fallback. If a future caller starts pre-imputing `R34_km` upstream, that caller would need to pass equivalent provenance into this function to preserve the same branching.
