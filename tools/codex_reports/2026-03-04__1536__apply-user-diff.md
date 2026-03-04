Goal

Apply the provided unified diff exactly to `R/hazard_core.R`, then run the requested targeted and full `devtools::test()` commands.

Scope

Only `R/hazard_core.R` was edited. No additional code or test files were modified by this turn.

Summary

Applied the requested hunk to `compute_site_winds_full()`:
- added diagnostic/traceability columns for effective climatological `R34` and effective clamped `RMW`
- routed both the Holland symmetric wind solver and the forward-motion asymmetry through `RMW_used_km`

Files changed

- `R/hazard_core.R`
- `tools/codex_reports/2026-03-04__1536__apply-user-diff.md`

Commands run

- `Get-Content -Path .agents/skills/r-coding/SKILL.md`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `Get-Content -Path R\\hazard_core.R`
- `rg -n "hazard_outer_cutoff|rmw|wind|holland" tests/testthat`
- `Get-Content -Path tests\\testthat\\test-hazard_outer_cutoff.R`
- `Get-Content -Path tests\\testthat\\test-rmw_estimation.R`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `R -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"` (failed in PowerShell because bare `R` resolved to `Invoke-History`)
- `R.exe -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"`
- `R.exe -q -e "devtools::test()"`
- `if (-not (Test-Path tools\\codex_reports)) { New-Item -ItemType Directory -Path tools\\codex_reports | Out-Null } ; Get-Date -Format "yyyy-MM-dd__HHmm"`
- `git status --short`

Test results

- Parse validation passed for `R/hazard_core.R`
- Filtered tests passed: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 12 ]`
- Full tests passed: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`

Behavior changes

- `compute_site_winds_full()` now returns additional intermediate columns: `lat0`, `R34_is_climo`, `R34_eff_km`, and `RMW_used_km`
- The effective clamped `RMW_used_km` is now explicitly used in the downstream wind solver and asymmetry adjustment

Follow-ups/risks

- The repository already had unrelated modified and untracked files before this change; they were left untouched.
- The added intermediate columns are user-visible output changes from `compute_site_winds_full()`, even though the public function signature is unchanged.
