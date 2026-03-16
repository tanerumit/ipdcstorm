# ipdc storm — Agent Instructions (canonical)

Generic R coding constraints (API stability, minimal diff, deps, style, validation) live in `.agents/skills/r-coding/SKILL.md` and always apply.

## Repo map

| Path | Contents |
|------|----------|
| `R/` | Package source |
| `tests/testthat/` | Unit tests |
| `man/` (generated) | Roxygen2 documentation |
| `vignettes/` | Long-form docs (Quarto) |
| `inst/extcode/pipelines/` | Entry-point pipeline scripts |
| `.agents/skills/` | Agent skill files |
| `dev/codex_reports/` | Completion reports |

## Pre-flight (mandatory, before any code)

Before writing or editing any file, complete these steps in order:

1. Read this file (`AGENTS.md`) in full.
2. Check out the `agent` branch (create from `master` if missing).
3. Scan `.agents/skills/` and identify which skills apply to the task.
4. Read every applicable `SKILL.md` and its `references/` folder.
5. Print a skill activation summary to the console (see below).

### Skill activation protocol

At the start of every task, print the following to the console:

```
── Skills ─────────────────────────────────────────────────
  always  : r-coding
  loaded  : <skill-1>, <skill-2>        ← task-relevant
  skipped : <skill-3>                   ← present but not relevant
───────────────────────────────────────────────────────────
```

Rules:
- `r-coding` is **always** loaded — it is not optional.
- A skill is **loaded** if its `SKILL.md` "When to Use" criteria match the current task. Read the skill's `references/` files as well.
- A skill is **skipped** if it exists in `.agents/skills/` but its activation criteria do not match.
- If no additional skills match beyond `r-coding`, print `loaded: (none)`.
- If a task pack references a specific skill by name, that skill is loaded regardless of its trigger criteria.

## Scientific correctness and transparency

- Disclose scientific/numerical assumptions and simplifications before code
  when explanation is requested.
- No silent simplifications of equations, physics, or algorithms.
- State key numerical/physical/model limitations when explanation is requested.

## Allowed dependencies

Core (already used throughout): `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `data.table`, `sf`.
Do not add anything outside this set without explicit approval.

## Branch and commit conventions

- **Never commit directly to `main`.** All agent work happens on the `agent` branch.
- If `agent` does not exist, create it from `master` before starting.
- Commit message: `<type>: <what changed>` — one line, max 72 chars.
  - Types: `fix`, `feat`, `refactor`, `docs`, `test`, `chore`
  - Example: `fix: correct Holland B parameter bounds check`
- One logical change per commit. Do not bundle unrelated fixes.
- Merge to `master` is manual (by the user), never by the agent.

## Pipeline compatibility after refactors 
- Any refactor that may affect downstream usage must be checked against inst/extcode/pipelines/. 
- If pipeline scripts are broken by changes to functions, configs, return objects, column names, or workflow resolution, patch the affected scripts within scope using minimal-diff edits. 
- Do not broaden changes to unrelated scripts.

## Completion gates and reporting (mandatory)

**Purpose:**
- Create an audit trail for reproducibility, accountability, and future context. Reports enable rapid onboarding to existing work, retrospective debugging, and clear communication of scope boundaries.

**Report structure:**
- Ensure `dev/codex_reports/` exists.
- Create exactly one report file per task: `dev/codex_reports/YYYY-MM-DD__HHMM__<task-slug>.md`

Required sections:

1. **Goal**: What outcome was intended.
   - If triggered by a task pack, include its path: `Task pack: dev/codex/YYYY-MM-DD_slug.task.md`
2. **Scope**: What was/wasn't included.
3. **Skills loaded**: List every skill that was activated for this task (must match the console output from pre-flight). Format:
   ```
   - always : r-coding
   - loaded : r-vignettes, task-workflow
   - skipped: (none)
   ```
4. **Problem solved**: Specific pain point or blocker this addresses.
5. **Summary**: Concise overview of approach and results.
6. **Files changed**: Paths, line counts, nature of edits.
7. **Commands run**: Exact reproduction steps.
8. **Test results**: Pass/fail status and coverage.
9. **Behavior changes**: User-facing or integration impacts.
10. **Follow-ups/risks**: Debt, edge cases, dependencies not yet addressed.