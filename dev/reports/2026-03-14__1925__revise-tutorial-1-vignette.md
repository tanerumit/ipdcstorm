# Goal

Revise `vignettes/tutorial_1_introduction_setup.qmd` using the vignette workflow so it is aligned with the current `ipdcstorm` API and renders successfully to HTML.

# Scope

Included:
- Rewrite of the existing tutorial vignette narrative and runnable code chunks.
- API alignment for configuration, model run, and daily downscaling examples.
- Quarto render validation to standalone HTML.

Not included:
- Changes to R source, tests, package metadata, or generated man files.
- New vignettes or renaming of the existing vignette.

# Skills loaded

- always : r-coding
- loaded : r-vignettes
- skipped: task-workflow

# Problem solved

The vignette content had drifted away from the exported package interface. It referenced outdated argument names and calling patterns, which made the tutorial misleading and at risk of failing during render.

# Summary

Reworked the tutorial into a compact quickstart centered on one stationary baseline run. Updated the text to explain the current workflow, restored executable examples that match the package API, clarified that climate is embedded in `hazard_cfg`, and documented the named-list return shape of `generate_daily_hazard_impact()`. Verified the revised vignette by rendering `vignettes/tutorial_1_introduction_setup.qmd` to HTML without errors.

# Files changed

- `vignettes/tutorial_1_introduction_setup.qmd` — 180 lines; substantive rewrite of prose, chunk labels, and code examples to match the current API and vignette conventions.
- `dev/codex_reports/2026-03-14__1925__revise-tutorial-1-vignette.md` — new task report.

# Commands run

- `git branch --list agent master`
- `Get-ChildItem -Name .agents\skills`
- `Get-Content AGENTS.md`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-coding\SKILL.md`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-vignette\SKILL.md`
- `Get-ChildItem C:\Users\taner\WS\shared-tools\01_skills\r-vignette\references`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-vignette\references\quarto-defaults.yml`
- `Get-Content C:\Users\taner\WS\shared-tools\01_skills\r-vignette\references\chunk-conventions.md`
- `Write-Output '── Skills ─────────────────────────────────────────────────'; Write-Output '  always  : r-coding'; Write-Output '  loaded  : r-vignettes'; Write-Output '  skipped : task-workflow'; Write-Output '───────────────────────────────────────────────────────────'`
- `Get-Content vignettes\tutorial_1_introduction_setup.qmd`
- `rg "make_hazard_cfg|run_hazard_model|generate_daily_hazard_impact|read_ibtracs_clean|compute_site_winds_full|make_storm_events|compute_annual_counts|compute_lambda_table|simulate_twolevel_counts" R tests man NAMESPACE`
- `rg --files vignettes`
- `git status --short`
- `Get-Content R\hazard_run.R | Select-String -Pattern "make_hazard_cfg <- function|run_hazard_model <- function" -Context 0,35`
- `Get-Content R\hazard_downscale.R | Select-String -Pattern "generate_daily_hazard_impact <- function" -Context 0,45`
- `git diff -- vignettes\tutorial_1_introduction_setup.qmd`
- `Get-Content man\make_hazard_cfg.Rd`
- `Get-Content man\make_climate_cfg.Rd`
- `Get-Content man\generate_daily_hazard_impact.Rd`
- `quarto render vignettes\tutorial_1_introduction_setup.qmd --to html`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `git diff --stat -- vignettes\tutorial_1_introduction_setup.qmd`
- `Get-Content vignettes\tutorial_1_introduction_setup.qmd | Measure-Object -Line`
- `git status --short`

# Test results

- Pass: `quarto render vignettes\tutorial_1_introduction_setup.qmd --to html`
- Coverage: vignette render only; no unit tests were run because no package source or test files were changed.

# Behavior changes

- The tutorial now uses the current exported interface: `make_climate_cfg()`, `make_hazard_cfg(historical_start_year = ..., simulation_years = ..., climate = ...)`, `run_hazard_model(..., seed = ..., verbose = FALSE)`, and `generate_daily_hazard_impact(damage = list(...))`.
- The daily downscaling example now correctly shows that `generate_daily_hazard_impact()` returns a named list, even for one location.
- The narrative is now a quickstart baseline workflow with a short explanation of how future climate configurations plug into the same pipeline.

# Follow-ups/risks

- The tutorial still uses the packaged demo IBTrACS file for speed; results are illustrative rather than analysis-grade.
- If the package climate API changes again, the future-climate section should be rechecked because it intentionally summarizes concepts without running full future scenarios.
