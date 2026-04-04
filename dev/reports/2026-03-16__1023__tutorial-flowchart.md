# Goal

Add one clear process flowchart to `inst/tutorials/tutorial_1_introduction_setup.qmd` near the `### What the model does` section so the three-step `ipdcstorm` workflow is easier to understand for climate-water experts new to hurricane modeling.

# Scope

Included:
- Read the tutorial section and relevant vignette/skill guidance.
- Cross-checked the workflow against the package tutorial narrative and the temporal downscaling source file.
- Inserted one Quarto-compatible workflow diagram in the target tutorial.

Excluded:
- Rewriting or reorganizing tutorial prose beyond the figure insertion.
- Changes to package code, other tutorials, or pipeline scripts.
- Fixing unrelated tutorial render issues discovered during validation.

# Skills loaded

- always : r-coding
- loaded : flowchart-creator, r-vignettes
- skipped: defuddle, json-canvas, obsidian-bases, obsidian-cli, obsidian-markdown, pptx, task-workflow, skill-creator, skill-installer

# Problem solved

The tutorial described the three workflows in prose but did not provide a compact visual map of how historical calibration, stochastic simulation, and temporal downscaling connect end to end.

# Summary

Inserted a single Mermaid flowchart directly under "The modeling scheme consists of three workflows:". The diagram summarizes the three main steps and the main user-facing substeps under each step without adding implementation-level detail or changing the surrounding tutorial explanation.

# Files changed

- `inst/tutorials/tutorial_1_introduction_setup.qmd` — 21 insertions, 13 deletions already present in the working tree relative to `HEAD`; this task added one Mermaid flowchart block near the workflow overview.
- `dev/codex_reports/2026-03-16__1023__tutorial-flowchart.md` — new task report required by repo policy.

# Commands run

```powershell
Get-Content -Path AGENTS.md -Raw
git branch --list agent
Get-ChildItem -Path .agents\skills -Recurse -Filter SKILL.md | Select-Object -ExpandProperty FullName
Get-ChildItem -Path 'C:\Users\taner\WS\_shared\skills' -Recurse -Filter SKILL.md | Where-Object { $_.FullName -match 'r-coding|r-vignettes|flowchart-creator' } | Select-Object -ExpandProperty FullName
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\flowchart-creator\SKILL.md' -Raw
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-coding\SKILL.md' -Raw
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-vignettes\SKILL.md' -Raw
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-vignettes\references\quarto-defaults.yml' -Raw
Get-Content -Path 'C:\Users\taner\WS\_shared\skills\r-vignettes\references\chunk-conventions.md' -Raw
@'
── Skills ─────────────────────────────────────────────────
  always  : r-coding
  loaded  : flowchart-creator, r-vignettes
  skipped : defuddle, json-canvas, obsidian-bases, obsidian-cli, obsidian-markdown, pptx, task-workflow, skill-creator, skill-installer
───────────────────────────────────────────────────────────
'@ | Write-Output
rg -n "Historical calibration|Stochastic simulation|Temporal downscaling|hazard_core|hazard_climate|hazard_downscale" -S R inst\extcode\pipelines inst\tutorials\tutorial_1_introduction_setup.qmd
rg -n "mermaid|flowchart|graph TD|::: \{=html\}" inst/tutorials vignettes R
git diff -- inst/tutorials/tutorial_1_introduction_setup.qmd
quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html --output-dir C:\Users\taner\AppData\Local\Temp\ipdcstorm-quarto-check --no-project
```

# Test results

- Quarto render: failed, but not because of the inserted flowchart.
- Failure encountered at tutorial chunk `cfg`: `make_hazard_cfg()` reported unused arguments `start_year = 1970, n_sim_years = 1000`.
- Quarto also reported a cleanup issue for `inst/tutorials/.quarto` after the failed render.

# Behavior changes

- The tutorial now includes a compact workflow figure that visually connects the three package workflows and their main substeps near the conceptual overview.

# Follow-ups/risks

- The tutorial currently has an unrelated render-time API mismatch in the `cfg` chunk that prevents full HTML validation.
- The working tree was already dirty before this task; this edit was layered onto the current version of the tutorial without reverting user changes.
