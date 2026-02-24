# AGENTS.md — weathergenr repo guidance

## Repo map
- Code: R/
- Tests: tests/testthat/
- Docs: man/ (generated), vignettes/
- Skills: .agents/skills/
- Reports: tools/codex_reports/

## Non-negotiable constraints
- Encoding: UTF-8 (no BOM).
- Minimal diff: only necessary changes; no style-only refactors; no renames/moves unless asked.
- Public API stability: do not change exported function names, arguments, defaults, return types, column names, or classes unless explicitly instructed.
- Dependencies: do not add new external dependencies; do not add runtime reliance on Suggests. If unavoidable, guard with requireNamespace() and justify.
- Style: avoid purrr; prefer explicit loops and base R idioms (unless the codebase already uses something else locally).
- Determinism: preserve deterministic behavior; control randomness/seeds explicitly.
- User-visible behavior: no silent changes to outputs/warnings/errors/messages; document any required change.
- Data integrity: do not silently drop columns/fields/rows.

## Required workflow
- Read relevant files first: R/, tests/testthat/, DESCRIPTION, NAMESPACE.
- Implement the smallest safe fix.
- If behavior/assumptions change: update tests and (if exported) roxygen.
- Do not conclude until validation is done (or explicitly skipped with reason).

## Required validation (minimum)
- Rscript -e "parse(file='R/<touched-file>.R')" for each touched R file
- If roxygen changed: regenerate docs and record it.

## Completion gates and reporting (mandatory)
- Ensure tools/codex_reports/ exists.
- Create/update exactly one report file:
  tools/codex_reports/YYYY-MM-DD__HHMM__<task-slug>.md
- Report must include: goal, scope, summary, files changed, commands run, test results, behavior changes, follow-ups/risks.