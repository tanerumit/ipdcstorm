# Goal

Create a new standalone Quarto tutorial at `inst/tutorials/tutorial_3_validation.qmd` that walks users through the stationary-baseline validation workflow using `make_validation_cfg()` and `run_validation_suite(out, cfg)`.
Task pack: inline user task pack in chat (no `dev/codex/*.task.md` path was provided).

# Scope

Included:
- New Tutorial 3 Quarto file.
- Standalone stationary hazard-model setup using the bundled demo dataset.
- Validation configuration, full-suite run, tier-by-tier walkthrough, inline plots, and saved-output summary.
- End-to-end Quarto render verification.

Excluded:
- Any edits to R source, tests, package APIs, or existing tutorials.
- Climate change workflow content.
- Separate demonstrations of individual tier validation functions.

# Skills loaded

- always : r-coding
- loaded : r-vignettes, task-workflow
- skipped: design-doc-mermaid, flowchart-creator, obsidian-skills, pptx

# Problem solved

The repo did not have a dedicated end-user tutorial for the model-validation workflow. Users needed a standalone walkthrough that explains what the three validation tiers test, how to run them with the current API, and how to interpret the outputs.

# Summary

Added a new 324-line tutorial aligned to the tone and structure of Tutorials 1 and 2. The tutorial recreates a baseline `out` object from the demo dataset, builds `val_cfg`, runs `run_validation_suite()`, shows the key tables for hindcast, rate, and wind-field diagnostics, and displays all six validation plots inline. It also documents the saved CSV/PNG outputs and includes brief guidance on what to check first when a tier fails.

# Files changed

- `inst/tutorials/tutorial_3_validation.qmd` — new file, 324 lines, full tutorial content.
- `inst/tutorials/tutorial_3_validation.html` — generated render artifact.

# Commands run

- `Get-Content -Path AGENTS.md -Raw`
- `git branch --list agent`
- `Get-ChildItem -Path .agents\skills -Force | Select-Object -ExpandProperty Name`
- `Get-Content -Path .agents\skills\r-coding\SKILL.md -Raw`
- `Get-Content -Path .agents\skills\r-coding\references\r-workflow.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\SKILL.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\references\chunk-conventions.md -Raw`
- `Get-Content -Path .agents\skills\r-vignettes\references\quarto-defaults.yml -Raw`
- `Get-Content -Path .agents\skills\task-workflow\SKILL.md -Raw`
- printed skill activation summary to console
- `Get-Content -Path inst\tutorials\tutorial_1_introduction_setup.qmd -Raw`
- `Get-Content -Path inst\tutorials\tutorial_2_multisite_validation_climate.qmd -Raw`
- `Get-Content -Path man\make_validation_cfg.Rd -Raw`
- `Get-Content -Path man\run_validation_suite.Rd -Raw`
- `Get-Content -Path man\validate_hazard_model.Rd -Raw`
- `Rscript -e "library(ipdcstorm); library(dplyr); ibtracs_demo <- system.file('extdata','ibtracs_demo.csv', package='ipdcstorm'); targets <- data.frame(name=c('St_Martin','Saba','Statia'), lat=c(18.0708,17.6350,17.4890), lon=c(-63.0501,-63.2300,-62.9740)); cfg <- make_hazard_cfg(data_path=ibtracs_demo, historical_start_year=1970L, search_radius_km=800, simulation_years=500L, climate=make_climate_cfg(scenario='stationary')); out <- run_hazard_model(cfg=cfg, targets=targets, seed=123L, verbose=FALSE); val_cfg <- make_validation_cfg(holdout_years=10L, n_sim=2000L, return_periods=c(5,10,25,50), conf_level=0.90, seed=123L, out_dir='output/validation_smoke', save_plots=TRUE, save_tables=TRUE); val <- run_validation_suite(out=out, cfg=val_cfg); cat('summary names\\n'); print(names(val)); cat('hindcast cols\\n'); print(names(val$hindcast$comparison)); cat('rate cols\\n'); print(names(val$rate_check)); cat('wind cols\\n'); print(names(val$wind_field)); cat('artifacts plots\\n'); print(val$artifacts$plots); cat('artifacts tables\\n'); print(val$artifacts$tables); cat('summary head\\n'); print(head(val$summary));"`
- `Get-ChildItem -Path output\validation_smoke | Select-Object Name,Length`
- `Get-Content -Path R\hazard_validation.R -Raw`
- `rg -n "rate_comparison|hindcast_return_levels|wind_field_scatter|bias_decomposition|qq_annual_max|exceedance_comparison|validation_summary.csv|run_validation_suite\(|make_validation_cfg\(" R\hazard_validation.R`
- `Get-Content -Path R\hazard_validation.R | Select-Object -Skip 3828 -First 45` (with line numbering)
- `Rscript -e "quarto::quarto_render('inst/tutorials/tutorial_3_validation.qmd')"` (rerun with escalation; passed)
- `Get-ChildItem inst\tutorials\tutorial_3_validation* | Select-Object Name,Length,LastWriteTime`
- `(Get-Content inst\tutorials\tutorial_3_validation.qmd | Measure-Object -Line).Lines`
- `Get-ChildItem output\validation | Select-Object Name,Length`

# Test results

- `quarto::quarto_render('inst/tutorials/tutorial_3_validation.qmd')` — PASS.
- Tutorial rendered successfully with `execute: eval: true`.
- The rendered run produced non-empty outputs for all three validation tiers on the demo dataset, so rollback criteria were not triggered.
- The rendered output directory contained the six expected plot filenames and four expected CSV files.

# Behavior changes

- Added a new user-facing validation tutorial at the root tutorial path.
- The tutorial uses the current typed validation API and displays all six validation plots inline.
- The tutorial includes a small runtime fallback to populate `val$artifacts$plots$rate_comparison` if the current suite does not register that path automatically, while leaving package code unchanged.

# Follow-ups/risks

- `run_validation_suite()` currently saves `rate_comparison.png` but does not register it in `val$artifacts$plots`; the tutorial compensates for this locally. If the package behavior is later fixed, that fallback becomes unnecessary but remains harmless.
- The validation run writes into `output/validation/`, which may accumulate artifacts from repeated tutorial renders unless cleaned separately.
