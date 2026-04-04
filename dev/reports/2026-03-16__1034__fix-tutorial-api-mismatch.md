# Goal

Fix runtime bugs and API mismatches in `inst/tutorials/tutorial_1_introduction_setup.qmd`, then render the tutorial successfully to HTML.

# Scope

Included:
- Checked the tutorial against the current `ipdcstorm` API.
- Updated stale tutorial calls and parameter names.
- Applied the smallest tutorial-side workaround needed for a package verbosity bug encountered during render.
- Rendered the tutorial to `inst/tutorials/tutorial_1_introduction_setup.html`.

Excluded:
- Changes to package source files under `R/`.
- Broader tutorial rewrites beyond the broken call sites and adjacent explanatory text.
- Fixes to unrelated files already modified in the working tree.

# Skills loaded

- always : r-coding
- loaded : r-vignettes
- skipped: defuddle, flowchart-creator, json-canvas, obsidian-bases, obsidian-cli, obsidian-markdown, pptx, task-workflow, skill-creator, skill-installer

# Problem solved

The tutorial referenced outdated function arguments and output assumptions, causing Quarto execution to fail during configuration and daily downscaling steps. A second failure came from a package-side verbose logging path that assumes `target_year` exists even for stationary climate runs.

# Summary

Updated the tutorial to use `historical_start_year` and `simulation_years` in `make_hazard_cfg()`, removed the obsolete `sst_cfg = NULL` argument from `run_hazard_model()`, passed `verbose = FALSE` to avoid the stationary-climate logging bug, and changed the daily hazard example to use `damage = list(method = "powerlaw")` plus extraction of the single-location tibble from the returned named list. The tutorial now renders successfully to standalone HTML.

# Files changed

- `inst/tutorials/tutorial_1_introduction_setup.qmd` — 44 insertions, 34 deletions; fixed stale API calls, corrected parameter names in nearby prose, and adjusted the daily example to the current return shape.
- `inst/tutorials/tutorial_1_introduction_setup.html` — rendered output regenerated successfully.
- `dev/codex_reports/2026-03-16__1034__fix-tutorial-api-mismatch.md` — new task report required by repo policy.

# Commands run

```powershell
Get-Content -Path AGENTS.md -Raw
git branch --list agent
@'
── Skills ─────────────────────────────────────────────────
  always  : r-coding
  loaded  : r-vignettes
  skipped : defuddle, flowchart-creator, json-canvas, obsidian-bases, obsidian-cli, obsidian-markdown, pptx, task-workflow, skill-creator, skill-installer
───────────────────────────────────────────────────────────
'@ | Write-Output
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-coding\SKILL.md' -Raw
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-vignettes\SKILL.md' -Raw
rg -n "make_hazard_cfg|run_hazard_model|run_daily_hazard|build_event_library|generate_daily|plot_|read_ibtracs|compute_" inst/tutorials/tutorial_1_introduction_setup.qmd
rg -n "make_hazard_cfg\s*<-|run_hazard_model\s*<-|run_daily_hazard\s*<-|build_event_library_from_out\s*<-|generate_daily_hazard_impact\s*<-|plot_" R
quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html
```

# Test results

- Quarto render: PASS
- Output file created: `inst/tutorials/tutorial_1_introduction_setup.html`
- No R source files were modified, so parse checks and unit tests for package code were not applicable to this task.

# Behavior changes

- The tutorial now uses the current configuration API and current daily hazard API.
- The daily example now correctly treats `generate_daily_hazard_impact()` output as a named list and extracts the `"Saba"` tibble before summarising.
- The rendered HTML tutorial is available alongside the `.qmd`.

# Follow-ups/risks

- The `verbose = FALSE` addition is a tutorial-side workaround for a package bug in `run_hazard_model()` logging for stationary climate runs. If desired, that should be fixed in package code separately.
- The target tutorial file already contained prior user edits in the working tree; this task preserved and built on those changes.
