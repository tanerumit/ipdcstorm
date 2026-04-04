# Goal

Revise Tutorial 2 into a focused multi-site baseline and climate-change walkthrough that uses the current climate API and renders end to end.
Task pack: inline user task pack in chat (no `dev/codex/*.task.md` path was provided).

# Scope

Included:
- Full rewrite of `inst/tutorials/tutorial_2_multisite_validation_climate.qmd`.
- Multi-site baseline example for St_Martin, Saba, and Statia.
- Climate comparison using `make_climate_cfg()` embedded in `make_hazard_cfg()`.
- Daily Saba comparison for baseline vs future scenario.
- End-to-end Quarto render verification.

Excluded:
- Validation content and validation reference tables.
- Any R source, tests, or pipeline script edits.
- Level 3 perturbation demonstration.

# Skills loaded

- always : r-coding
- loaded : r-vignettes, task-workflow
- skipped: design-doc-mermaid, flowchart-creator, obsidian-skills, pptx

# Problem solved

The previous Tutorial 2 content was outdated and referenced obsolete climate and daily-impact APIs (`make_sst_cfg()`, `sst_cfg`, and string-based `damage_method`). It also mixed in validation material that no longer belongs in this tutorial.

# Summary

Replaced Tutorial 2 with a concise Quarto tutorial aligned to Tutorial 1's tone and formatting. The new content covers a brief three-island baseline run, a side-by-side stationary vs `delta_sst = 1.5` climate comparison, one `ggplot2` activity chart, a daily Saba disruption comparison, compact output-column tables, and a short parameter guide for `delta_sst`, `scenario`, and `sensitivity_mode`.

# Files changed

- `inst/tutorials/tutorial_2_multisite_validation_climate.qmd` — new file, 254 lines, complete tutorial rewrite.
- `inst/tutorials/tutorial_2_multisite_validation_climate.html` — generated render artifact from Quarto.

# Commands run

- `Get-Content -Path AGENTS.md -Raw`
- `git branch --list agent`
- `Get-ChildItem -Path .agents\skills -Force | Select-Object -ExpandProperty Name`
- `Get-ChildItem -Path .agents\skills\r-coding -Recurse | Select-Object FullName`
- `Get-ChildItem -Path .agents\skills\r-vignettes -Recurse | Select-Object FullName`
- `Get-ChildItem -Path .agents\skills\task-workflow -Recurse | Select-Object FullName`
- `Get-Content -Path .agents\skills\r-coding\SKILL.md -Raw`
- `Get-Content -Path .agents\skills\r-coding\references\r-workflow.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\SKILL.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\references\chunk-conventions.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\references\quarto-defaults.yml -Raw`
- `Get-Content -Path .agents\skills\task-workflow\SKILL.md -Raw`
- printed skill activation summary to console
- `rg --files . | rg "02-climate-change\.R|tutorial_2_multisite_validation_climate\.qmd$|sst_scenario_info|make_climate_cfg|generate_daily_hazard_impact"`
- `Get-ChildItem -Path inst\tutorials -Force`
- `Get-Content -Path inst\tutorials\tutorial_1_introduction_setup.qmd -Raw`
- `Get-Content -Path inst\tutorials\repo2\tutorial_2_multisite_validation_climate.qmd -Raw`
- `Get-Content -Path inst\extcode\02-climate-change.R -Raw`
- `Get-Content -Path man\generate_daily_hazard_impact.Rd -Raw`
- `Get-Content -Path man\make_climate_cfg.Rd -Raw`
- `Get-Content -Path man\sst_scenario_info.Rd -Raw`
- `rg -n "n_tc|n_total|p_hur|p_hurricane|simulation_years|historical_start_year|make_hazard_cfg\(|run_hazard_model\(" R man inst\extcode`
- `Get-Content -Path man\make_hazard_cfg.Rd -Raw`
- `Get-Content -Path man\run_hazard_model.Rd -Raw`
- `Rscript -e "library(ipdcstorm); targets <- data.frame(name=c('St_Martin','Saba','Statia'), lat=c(18.0708,17.6350,17.4890), lon=c(-63.0501,-63.2300,-62.9740)); cfg <- make_hazard_cfg(data_path = system.file('extdata','ibtracs_demo.csv', package='ipdcstorm'), historical_start_year = 1970L, simulation_years = 5L, climate = make_climate_cfg(scenario='stationary')); out <- run_hazard_model(cfg = cfg, targets = targets, seed = 123, verbose = FALSE); print(names(out$sim)); print(head(out$sim));"`
- `R -q -e "quarto::quarto_render('inst/tutorials/tutorial_2_multisite_validation_climate.qmd')"` (failed in PowerShell argument parsing)
- `Rscript -e "quarto::quarto_render('inst/tutorials/tutorial_2_multisite_validation_climate.qmd')"` (sandbox failure)
- `Rscript -e "quarto::quarto_render('inst/tutorials/tutorial_2_multisite_validation_climate.qmd')"` (rerun with escalation; passed)
- `Rscript -e "quarto::quarto_render('inst/tutorials/tutorial_2_multisite_validation_climate.qmd')"` (final rerender after style pass; passed)

# Test results

- `quarto::quarto_render('inst/tutorials/tutorial_2_multisite_validation_climate.qmd')` — PASS.
- Tutorial rendered to `inst/tutorials/tutorial_2_multisite_validation_climate.html` with `execute: eval: true`.
- Coverage not applicable for vignette-only changes.

# Behavior changes

- Tutorial 2 now targets the current climate workflow centered on `make_climate_cfg()` and `make_hazard_cfg(climate = ...)`.
- Validation sections and reference tables were removed from Tutorial 2.
- Daily impact examples now use `damage = list(method = "intensity")` and extract from the named list returned by `generate_daily_hazard_impact()`.
- The rendered tutorial is now available at the root tutorial path instead of only under `inst/tutorials/repo2/`.

# Follow-ups/risks

- The current package runtime exposes annual simulation columns such as `n_total` and `p_hurricane`; the tutorial summarizes those into `mean_total` and `mean_p_hur` tables to keep the narrative clear.
- The old `inst/tutorials/repo2/` copy remains in the repository and may still confuse future readers if not retired separately.
