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
| `tools/codex_reports/` | Completion reports |

## Scientific correctness and transparency

- Disclose scientific/numerical assumptions and simplifications before code when explanation is requested.
- No silent simplifications of equations/physics/algorithms.
- State key numerical/physical/model limitations when explanation is requested.
- Climate modification levels follow a three-tier hierarchy (L1 → SST-rate, L2 → intensity shift, L3 → storm perturbation) defined in `R/hazard_climate.R`.
- Internal units: wind in m/s, distances in km, pressure in hPa. Convert at output boundaries only.

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


## Completion gates and reporting (mandatory)

- Ensure `tools/codex_reports/` exists.
- Create exactly one report file per task: `tools/codex_reports/YYYY-MM-DD__HHMM__<task-slug>.md`
- Report must include: goal, scope, summary, files changed, commands run, test results, behavior changes, follow-ups/risks.
