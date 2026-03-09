---
name: r-coding
description: Implement, debug, test, and refactor R code and R packages. Use for R functions, failing tests, performance issues, roxygen2 documentation, dependency hygiene, and package checks.
---

# R Coding Skill

R-specific implementation tactics. Repo-wide rules (API stability, minimal-diff policy, reporting, required validation, encoding, etc.) live in AGENTS.md and always apply.

## When to use

- Implement or refactor R functions (base R / tidyverse).
- Diagnose and fix failing testthat tests.
- Improve performance (vectorization, preallocation, algorithmic complexity).
- Update roxygen2 docs and regenerate Rd/NAMESPACE.
- Resolve package-check issues (imports, NOTES/WARNINGS, Rd problems).

## HARD CONSTRAINTS (NON-NEGOTIABLE)

- UTF-8 only (no BOM).
- No new dependencies. Do not require Suggests packages at runtime. If unavoidable, guard with requireNamespace() and justify.
- No purrr / functional-programming rewrites. Use explicit for loops and base R idioms.
- Deterministic behavior. Control randomness with explicit seeds.
- Minimal diff. Modify only what is necessary — no sweeping refactors, style-only reformatting, or file renames unless instructed.
- No silent changes to user-visible outputs, warnings, errors, columns, fields, or rows. Document any necessary changes in the final summary.
- No backwards compatibility unless explicitly requested.

## Working approach

1. Start from a minimal reproduction (single function call or single test file).
2. Read the smallest set of relevant files first (implementation + tests + package metadata).
3. Prefer the narrowest safe fix; avoid broad rewrites unless requested.
4. Stay locally consistent with conventions in the touched files.

## Implementation heuristics

- **Argument validation:** fail fast with specific error messages; validate type, length, and range.
- **Side effects:** isolate I/O from pure computation; avoid hidden global state.
- **NSE:** avoid non-standard evaluation unless the codebase already uses it for that API.

## Error handling & conditions

- `stop()` for unrecoverable problems, `warning()` for recoverable issues, `message()` for diagnostics only — never for results.
- Use `tryCatch()` with specific condition classes; never `try()`. Scope handlers to the narrowest operation possible.
- Match on condition *class*, not message text — messages change across R versions and locales.
- Error messages must state what was expected, what was received, and which argument failed.

## NAMESPACE hygiene

- Never call `library()` or `require()` inside package code (`R/`).
- Use `@importFrom pkg fun` for frequent calls; `pkg::fun()` for one-offs. Avoid blanket `@import pkg`.
- After any NAMESPACE change, run `devtools::document()` and verify — do not hand-edit `NAMESPACE`.
- When removing a dependency, grep `R/` and `tests/` for stale `::` calls and `@importFrom` directives.

## S3 method discipline

- Name methods `generic.class`; avoid dots in class names (dispatch ambiguity).
- Always register methods via `@export` or explicit `S3method()` — unregistered methods break under namespaced calls.
- Method formals must match the generic exactly (including `...`).
- Assign class in the constructor (`structure()` or `class<-`); provide at minimum a `print` method for user-facing classes.

## Numerical and data-edge safety

- Handle `NA`/`NaN`/`Inf` intentionally (don't rely on accidental propagation).
- Make assumptions explicit when thresholds, scaling, or aggregation are involved.
- Division by zero / undefined baselines: choose conservative, documented behavior.

## Performance checklist

- Avoid quadratic loops on time-series unless unavoidable; justify if O(n²).
- Preallocate vectors/matrices in loops; avoid repeated `rbind/cbind` in loops.
- Prefer `vapply()` over `sapply()` when shapes matter; use `for` loops when clarity wins.
- Subset large objects once; reuse — avoid unnecessary copies.

## Roxygen2 guidance (exported functions)

- `@description`: one precise paragraph.
- `@param`: type, length/dimensions, units, NA handling.
- `@return`: exact structure, classes, names, units, and invariants.
- `@examples`: minimal and fast; wrap expensive examples in `\dontrun{}`.

## Testing guidance

- Add/adjust tests when behavior or assumptions change.
- Use explicit seeds for any randomness; prefer deterministic tests.
- Cover: input validation, NA/empty/constant cases, return structure (class/names/dims), and key edge conditions.
- Group tests by function or function family; keep fixtures small and local.

## Reference

If present in the repo, consult: `references/r-workflow.md`.
