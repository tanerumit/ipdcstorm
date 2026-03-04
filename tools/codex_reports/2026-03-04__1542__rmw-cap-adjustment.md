Goal

Apply the provided unified diff exactly to the repository and validate the change with the requested test runs.

Scope

Minimal change in `R/hazard_core.R` only. No API, documentation, dependency, or test-file changes.

Summary

Updated the internal `rmw_over_r34_cap` constant from `3.0` to `4.0` in two places:
- within `.estimate_site_wind_holland()`
- within `compute_site_winds_full()`

Files Changed

- `R/hazard_core.R`
- `tools/codex_reports/2026-03-04__1542__rmw-cap-adjustment.md`

Commands Run

- `Get-Content -Path .agents/skills/r-coding/SKILL.md`
- `Get-Content -Path R\hazard_core.R`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `rg -n "hazard_outer_cutoff|rmw|wind|holland" tests\testthat`
- `Get-Content -Path tests\testthat\test-hazard_outer_cutoff.R`
- `Get-Content -Path tests\testthat\test-rmw_estimation.R`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `R -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"` attempted first, but in PowerShell `R` is an alias to `Invoke-History`
- `Get-Command R`
- `Get-Command Rscript`
- `Get-Command R.exe`
- `R.exe -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"`
- `R.exe -q -e "devtools::test()"`

Test Results

- Parse check passed for `R/hazard_core.R`.
- Filtered tests passed: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 12 ]`
- Full test suite passed: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`

Behavior Changes

Yes. For rows using climatological `R34`, the internal clamp now allows a larger `RMW` by capping at `R34 / 4.0` instead of `R34 / 3.0`.

Follow-ups/Risks

- No additional risks found from the requested change.
- The direct `R -q -e ...` form is not usable in this PowerShell environment because `R` is aliased; `R.exe` is required.
