Goal
- Apply the user-provided unified diff exactly to `R/hazard_core.R`, then run `R.exe -q -e "devtools::test()"`.

Scope
- Read the relevant package metadata and target implementation.
- Apply only the requested hunk in `compute_site_winds_full()`.
- Validate syntax for the touched file.
- Run the full package test suite.

Summary
- Updated `compute_site_winds_full()` so `R34_nm` now uses a quality gate based on the count of finite quadrant radii:
- Use directional `R34` when at least 3 quadrants are present.
- Fall back to the mean `R34` when 2 quadrants are present.
- Leave `R34` missing when fewer than 2 quadrants are present.
- Reused the computed `R34_nm_mean` for `R34_mean_km`.

Files changed
- `R/hazard_core.R`
- `tools/codex_reports/2026-03-04__1841__apply-provided-diff.md`

Commands run
- `Get-Content -Path .agents/skills/r-coding/SKILL.md -TotalCount 220`
- `Get-Content -Path DESCRIPTION -TotalCount 220`
- `Get-Content -Path NAMESPACE -TotalCount 220`
- `Get-Content -Path R\hazard_core.R -TotalCount 220`
- `rg -n "compute_site_winds_full|R34_mean_km|R34_nm" tests/testthat R -g "*.R"`
- `Get-ChildItem -Path tools\codex_reports -Force | Select-Object -First 5 | ForEach-Object { $_.Name }`
- `$i=1; Get-Content -Path R\hazard_core.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 710 -First 70`
- `Get-Content -Path tests\testthat\test-rmw_estimation.R -TotalCount 220`
- `Rscript -e "parse(file='R/hazard_core.R')"`
- `R.exe -q -e "devtools::test()"`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `$i=1; Get-Content -Path R\hazard_core.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 728 -First 22`

Test results
- `Rscript -e "parse(file='R/hazard_core.R')"`: passed.
- `R.exe -q -e "devtools::test()"`: passed, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 31 ]`.

Behavior changes
- For sparse `R34` quadrant input, `compute_site_winds_full()` no longer always uses the directional quadrant radius.
- With exactly 2 finite `R34` quadrants, it now uses the mean radius instead of a single directional value.
- With 0 or 1 finite `R34` quadrants, it now leaves `R34_nm` missing so downstream climatological infill can handle it.

Follow-ups/risks
- No targeted test was added for the new 0/1/2-quadrant `R34` gating behavior because the request was to apply the provided diff exactly and the existing full test suite passed.
- The new intermediate columns `nq34`, `R34_nm_dir`, and `R34_nm_mean` remain in the returned tibble because the requested diff introduces them within `mutate()` and does not remove them.
