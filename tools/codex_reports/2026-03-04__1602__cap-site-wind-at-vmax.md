# Goal

Apply the provided unified diff exactly to cap `V_site_kt` at `Vmax_kt` after forward-motion asymmetry in `compute_site_winds_full()`, then run the requested tests.

# Scope

- Read relevant package files and the target implementation file.
- Apply the minimal requested change in `R/hazard_core.R`.
- Validate syntax for the touched R file.
- Run the two requested `devtools::test()` commands.

# Summary

Added a post-asymmetry physical guard in `compute_site_winds_full()` so final site wind remains bounded by storm `Vmax_kt` when both values are finite. This preserves the existing symmetric Holland cap and extends the same constraint to the final `V_site_kt`.

# Files Changed

- `R/hazard_core.R`
- `tools/codex_reports/2026-03-04__1602__cap-site-wind-at-vmax.md`

# Commands Run

- `Get-Content -Path .agents/skills/r-coding/SKILL.md -TotalCount 220`
- `Get-Content -Path DESCRIPTION -TotalCount 220`
- `Get-Content -Path NAMESPACE -TotalCount 220`
- `Get-Content -Path R\hazard_core.R -TotalCount 260`
- `rg -n "compute_site_winds_full|V_site_symmetric_kt|V_site_kt = \.add_forward_motion_asymmetry" R\hazard_core.R`
- `rg -n "hazard_outer_cutoff|rmw|wind|holland" tests\testthat`
- `$i=1; Get-Content -Path R\hazard_core.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 790 -First 40`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `R.exe -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"`
- `R.exe -q -e "devtools::test()"`

# Test Results

- `Rscript -e "parse(file='R/hazard_core.R')"`: passed
- `R.exe -q -e "devtools::test(filter = 'hazard_outer_cutoff|rmw|wind|holland')"`: passed, `PASS 12`, `FAIL 0`, `WARN 0`, `SKIP 0`
- `R.exe -q -e "devtools::test()"`: passed, `PASS 16`, `FAIL 0`, `WARN 0`, `SKIP 0`

# Behavior Changes

- `compute_site_winds_full()` now caps final `V_site_kt` at `Vmax_kt` after forward-motion asymmetry is added, when both values are finite.
- No API, arguments, defaults, return types, or column names changed.

# Follow-Ups/Risks

- No dedicated test was added for the new post-asymmetry cap because the request was to apply the supplied diff exactly. If this guard becomes a stable requirement, add a focused regression test covering forward-motion enhancement exceeding `Vmax_kt`.
