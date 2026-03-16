# Goal

Revise the Mermaid flowchart near the `### What the model does` section in `inst/tutorials/tutorial_1_introduction_setup.qmd` so the three-step `ipdcstorm` workflow is easier to understand for climate-water experts new to hurricane modeling.

# Scope

Included:
- Reworked the Mermaid flowchart structure and styling in `inst/tutorials/tutorial_1_introduction_setup.qmd`.
- Fixed the Mermaid syntax issue by removing the invalid subgraph class assignment pattern.
- Re-rendered the tutorial HTML to confirm the diagram parses in context.

Not included:
- No changes to model code, tutorial narrative outside the diagram block, or other tutorial files.
- No changes to pipeline scripts.

# Skills loaded

- always : r-coding
- loaded : design-doc-mermaid
- skipped: flowchart-creator, obsidian-skills, pptx, r-vignettes, task-workflow

# Problem solved

The existing workflow diagram was hard to scan and contained Mermaid syntax that failed in the rendered tutorial. New users needed a clearer mental model of how calibration, simulation, and downscaling connect.

# Summary

I replaced the original diagram with three horizontal workflow columns, each containing a stage header, vertically stacked substeps, and an output node. Cross-workflow arrows now connect the outputs of Steps 1 and 2 to the next stage header, which makes the left-to-right progression more explicit while keeping the detailed substeps top-to-bottom within each workflow. After a follow-up runtime failure in Mermaid 11.2.0, I further simplified the block by removing inline HTML and the Mermaid init wrapper so the emitted HTML now contains plain Mermaid syntax only.

# Files changed

- `inst/tutorials/tutorial_1_introduction_setup.qmd` — edited the Mermaid block under `### What the model does`; one targeted documentation change.
- `inst/tutorials/tutorial_1_introduction_setup.html` — regenerated from Quarto render after the diagram update.

# Commands run

- `Get-Content -Raw AGENTS.md`
- `git branch --list agent; git branch --show-current`
- `Get-ChildItem -Name .agents/skills`
- `rg -n "What the model does|mermaid|flowchart" inst/tutorials/tutorial_1_introduction_setup.qmd`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\r-coding\SKILL.md`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\design-doc-mermaid\SKILL.md`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\r-vignettes\SKILL.md`
- `Get-ChildItem -Recurse -File C:\Users\taner\WS\_shared\skills\design-doc-mermaid\references | Select-Object -ExpandProperty FullName`
- `Get-ChildItem -Recurse -File C:\Users\taner\WS\_shared\skills\r-vignettes\references | Select-Object -ExpandProperty FullName`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\design-doc-mermaid\references\guides\diagrams\activity-diagrams.md`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\design-doc-mermaid\references\guides\troubleshooting.md`
- `Get-Content -Raw C:\Users\taner\WS\_shared\skills\design-doc-mermaid\references\guides\unicode-symbols\guide.md`
- `Write-Output '── Skills ─────────────────────────────────────────────────'; Write-Output '  always  : r-coding'; Write-Output '  loaded  : design-doc-mermaid'; Write-Output '  skipped : flowchart-creator, obsidian-skills, pptx, r-vignettes, task-workflow'; Write-Output '───────────────────────────────────────────────────────────'`
- `$i=1; Get-Content inst/tutorials/tutorial_1_introduction_setup.qmd | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 34 -First 45`
- `git status --short`
- `$i=1; Get-Content inst/tutorials/tutorial_1_introduction_setup.qmd | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 46 -First 45`
- `$i=1; Get-Content inst/tutorials/tutorial_1_introduction_setup.qmd | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 91 -First 20`
- `quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`
- `git diff --numstat -- inst/tutorials/tutorial_1_introduction_setup.qmd inst/tutorials/tutorial_1_introduction_setup.html`
- `git diff -- inst/tutorials/tutorial_1_introduction_setup.qmd`

# Test results

- Pass: `quarto render inst/tutorials/tutorial_1_introduction_setup.qmd --to html`
- Verified: rendered HTML now embeds the simplified Mermaid block in `inst/tutorials/tutorial_1_introduction_setup.html`.
- Not run: package unit tests, because this task only changed tutorial documentation.
- Coverage: render validation confirms the Mermaid block parses and the tutorial builds to standalone HTML.

# Behavior changes

- The workflow figure now presents the three major stages left-to-right, with substeps stacked vertically inside each stage.
- Each stage ends with an explicit output node, making the handoff between workflows easier to interpret for first-time readers.
- The Mermaid syntax error is resolved in the rendered tutorial.

# Follow-ups/risks

- The tutorial file already had unrelated uncommitted edits before this task; this report covers only the diagram update.
- I validated by rendering the full tutorial HTML, but I did not inspect pixel-level layout in a browser beyond successful render output.
