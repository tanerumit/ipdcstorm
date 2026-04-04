## Goal

Revise the flowchart in `inst/tutorials/tutorial_1_introduction_setup.qmd` near the `### What the model does` section so the three-step `ipdcstorm` workflow is easier to understand for climate-water experts new to hurricane modeling, with explicit links between substeps.

## Scope

Included:
- Kept the three main stages arranged horizontally.
- Connected substeps vertically within each stage.
- Added explicit handoff links from the last substep of one stage to the first substep of the next.
- Rendered the tutorial HTML to validate the updated Mermaid diagram.

Not included:
- No changes to the surrounding explanatory text.
- No changes to R code, tests, or pipeline scripts.

## Skills loaded

- always : r-coding
- loaded : flowchart-creator
- skipped: r-vignettes, task-workflow, defuddle, json-canvas, obsidian-bases, obsidian-cli, obsidian-markdown, pptx, skill-creator, skill-installer

## Problem solved

The prior revision grouped substeps visually inside each stage, but the workflow path itself was only shown at the stage level. That left the order of operations between substeps implicit rather than explicit.

## Summary

The Mermaid flowchart now shows a continuous process path at the substep level: `1.1 -> 1.2 -> 1.3 -> 2.1 -> 2.2 -> 2.3 -> 3.1 -> 3.2 -> 3.3`. The main stages remain as horizontal containers, while the visible arrows now convey both within-stage progression and cross-stage handoff.

## Files changed

- `inst/tutorials/tutorial_1_introduction_setup.qmd` — 87 insertions, 34 deletions. Added vertical substep links and inter-stage substep connections in the Mermaid block.
- `inst/tutorials/tutorial_1_introduction_setup.html` — 6402 insertions, 3827 deletions. Regenerated tutorial output after rendering.
- `dev/codex_reports/2026-03-16__1103__tutorial-flowchart-substep-links.md` — new file. Added task audit report.

## Commands run

- `git branch --list agent`
- `Get-Content -Raw "C:\Users\taner\WS\_shared\skills\r-coding\SKILL.md"`
- `Get-Content -Raw "C:\Users\taner\WS\_shared\skills\flowchart-creator\SKILL.md"`
- `$i=1; Get-Content inst\tutorials\tutorial_1_introduction_setup.qmd | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 44 -First 44`
- `@' ... '@ | Write-Output` (skill activation summary)
- `quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `git diff --numstat -- inst/tutorials/tutorial_1_introduction_setup.qmd inst/tutorials/tutorial_1_introduction_setup.html`

## Test results

- Pass: `quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html`
- No R source files were changed, so R parse checks and unit tests were not applicable for this task.

## Behavior changes

- The tutorial diagram now displays explicit arrows between substeps, not only between the three main stages.
- Readers can follow the end-to-end workflow order directly from the visual path.

## Follow-ups/risks

- Mermaid may route cross-subgraph arrows slightly differently across versions, though the current render completed successfully in this environment.
- If the diagram later gains more substeps, spacing may need to be revisited to keep the layout readable.
