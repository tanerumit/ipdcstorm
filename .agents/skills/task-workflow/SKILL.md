---
name: task-workflow
description: Conversation-to-task packaging conventions. Use when the user says "task pack", "patch only", or asks for a prompt file, or diff-only output. Governs how chat discussions are converted into actionable agent instructions.
---

# Task Workflow Skill

Rules for converting chat discussions into structured agent-ready deliverables.

---

## Trigger: "task pack" → Task Pack

When the user says exactly: **task pack**

Output only a lean markdown Task Pack (no preface) as a downloadable file:
`dev/codex/YYYY-MM-DD_<short-slug>.task.md` (Europe/Amsterdam date)

Use these headings (exact):

```
Task Pack — <title>
```

### Context
At the top of the context, place this line:

> Follow AGENTS.md in the repo root as the canonical ruleset for scope/API, scientific/numerical standards, deps, style, docs, error handling, I/O, and output modes.

Then include 2–6 bullets of key assumptions (units, key physics/params impacted, data selection rules if relevant).

### Goal
### Non-goals
### Allowed files

### Required changes (checklist)

### Tests to run
- Prefer scoped: `R -q -e "devtools::test(filter = '<REGEX>')"`
- Fallback: `R -q -e "devtools::test()"`
- If the filter is unclear, include `<REGEX>` as a placeholder — do not guess.
- Include a validation plan if results may change (baseline comparison, metrics/plots, acceptance thresholds where possible).

### Acceptance criteria
- Include rollback criteria when the change is risky: what outcome would cause revert.

### Output requirements
- If results are expected to change, include a **Results delta** subsection: what changed + why.

### Constraints
- Keep it tight (~200–600 tokens). No full code unless needed to remove ambiguity.
- **Scope enforcement:** if file paths or function names are missing, allow only files/functions explicitly mentioned earlier in the chat.

---

## Trigger: "patch only" → Diff only

When the user says exactly: **patch only**

Respond with only a unified diff (no prose), plus the exact test command(s) to run.

---

## Repo defaults

- Code in `R/`, tests in `tests/testthat/`.
- Agent work on the `draft` branch, never `master`.
