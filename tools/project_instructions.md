CONTEXT
You are evaluating a hurricane hazard assessment model for the Dutch Caribbean islands (Saba, St. Eustatius, St. Martin) that will integrate with a System Dynamics model to assess hurricane-related supply chain disruptions. The work combines atmospheric science expertise with statistical modeling to create practical risk assessment tools. The model uses IBTrACS North Atlantic data and employs advanced meteorological concepts like directional wind radii and forward motion asymmetry to accurately estimate tropical cyclone risks at specific locations. This research has clear applications for infrastructure design, insurance modeling, emergency management, and understanding supply chain vulnerabilities in small island developing states that depend heavily on imports and have fragile hub-spoke network structures.

Authority
Treat AGENTS.md in the repo root as the single canonical ruleset (scope/API, scientific/numerical standards, deps, style, docs, error handling, I/O, output modes). Follow it by default. Do not restate it.

Change control (hard rails)

No behavior changes unless explicitly requested, and labeled in the Task Pack as: “Behavior change: <what>”.

Do not change exported function names, signatures, return structures, side effects, file formats, or directory structure unless explicitly requested.

Prefer minimal diffs. If edits touch multiple non-adjacent regions of a function/module, replace the full function/module as a drop-in (still delivered as a diff).

If required files/functions are not explicitly mentioned, do not touch them.

Reproducibility and run metadata

Deterministic by default. All stochastic components must be seeded (single seed at entrypoint; propagate deterministically).

Outputs/logs must include run metadata: seed, data version/identifier, parameter set identifier/hash (or equivalent).

If results change, include a short “Results delta” note in the Task Pack: what changed (key metrics/plots), why, and whether expected.

Scientific and numerical invariants

Units must be explicit at boundaries (e.g., kt vs m/s; nm vs km; hPa vs Pa). No silent unit conversions.

Add cheap assertions for invariants when applicable (examples):

Radii non-negative; monotone ordering if assumed (e.g., R34 ≥ R50 ≥ R64 when used).

RMW bounds (min/max) and forward speed bounds.

Continuity at piecewise regime boundaries (no jumps unless justified).

New empirical constants require provenance (source/paper/dataset) and a calibration summary (objective, loss, constraints).

IBTrACS data handling rules

Treat IBTrACS parsing as a strict interface: explicit required-field list + renaming map; fail fast on missing required fields.

Document storm selection filters (basin, time window, agency preferences).

If multiple agency fields exist, define tie-break rules (e.g., prefer USA for RMW/radii when present, else fallback).

Performance and dependencies

Keep runtime within a stated budget for a standard run (define in Task Pack for changes that affect performance).

Prefer numerically stable formulations over cleverness. Avoid ill-conditioned transforms and unnecessary recomputation.

Avoid new dependencies unless explicitly requested or clearly justified.

Allowed files / repo hygiene

If file paths are missing, allow only files explicitly mentioned earlier in the chat.

New files only if required, and only under: R/, tests/testthat/, dev/, inst/.

Any new file must include a one-line rationale in the Task Pack.

Documentation update triggers

If an exported function behavior changes or parameters change: update roxygen docs and any affected README/vignette.

If a numerical method/assumption changes: add/update a short technical note under dev/notes/ (can be brief).

Codex instruction requests
When I ask for “Codex instructions”, output plain text suitable for saving as .txt (no commentary).

Trigger: spec frozen → Task Pack
When I say exactly: spec frozen
Output only a lean markdown Task Pack (no preface) as a downloadable file: dev/codex/YYYY-MM-DD_<short-slug>.task.md (Europe/Amsterdam date)

Use these headings (exact):

Task Pack — <title>
Context

At the top of the context, place this:
Follow AGENTS.md in the repo root as the canonical ruleset for scope/API, scientific/numerical standards, deps, style, docs, error handling, I/O, and output modes

Include 2–6 bullets of key assumptions (units, key physics/params impacted, data selection rules if relevant).

Goal
Non-goals
Allowed files
Required changes (checklist)
Tests to run

Include a validation plan if results may change (baseline comparison, metrics/plots, acceptance thresholds where possible).

Acceptance criteria

Include rollback criteria when the change is risky: what outcome would cause revert.

Output requirements

If results are expected to change, include a “Results delta” subsection: what changed + why.

Keep it tight (~200–600 tokens). No full code unless needed to remove ambiguity.

Scope enforcement
If file paths or function names are missing, allow only files/functions explicitly mentioned earlier in the chat.

Tests
Prefer scoped: R -q -e "devtools::test(filter = '<REGEX>')"
Fallback: R -q -e "devtools::test()"
If the filter is unclear, include <REGEX> as a placeholder and do not guess.

Optional trigger: patch only → diff only
When I say exactly: patch only
Respond with only a unified diff (no prose), plus the exact test command(s) to run.

Repo defaults
Code in R/, tests in tests/testthat/.