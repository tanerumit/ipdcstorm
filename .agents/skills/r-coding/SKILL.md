---
name: r-coding
description: Implement, debug, test, and refactor R code and R packages. Use for R functions, failing tests, performance issues, roxygen2 documentation, dependency hygiene, and package checks.
---

# R Coding Skill

- R-specific implementation tactics for package development.
- Repo-specific rules (project layout, scientific domain, reporting format) live in AGENTS.md and always apply on top of this skill.
- No automated linting enforced. Stay locally consistent with the touched file.
- New files: follow tidyverse style guide (snake_case, 2-space indent, `<-` for assignment).
- Do not reformat untouched code.

## When to use

- Implement or refactor R functions (base R / tidyverse).
- Diagnose and fix failing testthat tests.
- Improve performance (vectorization, preallocation, algorithmic complexity).
- Update roxygen2 docs and regenerate Rd/NAMESPACE.
- Resolve package-check issues (imports, NOTES/WARNINGS, Rd problems).

---

## HARD CONSTRAINTS (NON-NEGOTIABLE)

These apply to every change. No exceptions without explicit user instruction.

### Change discipline
- **Minimal diff.** Touch only what is necessary. No style-only reformats, renames, or file moves.
- **Public API stability.** Do not change exported function names, arguments, defaults, return types, column names, or classes.
- **No silent behavior changes.** Any change to user-visible outputs, warnings, errors, or messages must be documented in the final summary.

### Code style
- **No purrr / functional-programming rewrites.** Use explicit `for` loops and base R idioms unless the touched file already uses a different convention locally.
- **UTF-8 only** (no BOM).

### Dependencies
- **No new dependencies.** Do not require Suggests packages at runtime.
- If unavoidable: guard with `requireNamespace()`, justify in summary.

### Determinism
- Preserve deterministic behavior. Control randomness with explicit seeds.
- Fixed seeds for all stochastic code unless the user explicitly requests random behavior.

### Data integrity
- Do not silently drop columns, fields, or rows.
- Handle `NA`/`NaN`/`Inf` intentionally — never rely on accidental propagation.

---

## Required workflow

1. **Read first.** Open the narrowest set of relevant files: implementation, tests, DESCRIPTION, NAMESPACE.
2. **Smallest safe fix.** Prefer the narrowest change that solves the problem; avoid broad rewrites unless requested.
3. **Stay locally consistent** with conventions in the touched files.
4. **Update co-artifacts.** If behavior or assumptions change → update tests. If exported → update roxygen.
5. **Validate before concluding** (see Required Validation below). Do not mark done until validation passes or is explicitly skipped with reason.

---

## Required validation (minimum)

Run these for every change. Report results.

```bash
# 1. Parse check — every touched R file
Rscript -e "parse(file='R/<touched-file>.R')"

# 2. Unit tests — if a matching test file exists
Rscript -e "testthat::test_file('tests/testthat/test_<touched>.R')"

# 3. Roxygen — if any documentation changed
Rscript -e "devtools::document()"
# Then verify NAMESPACE and Rd files are consistent.
```

If a validation step fails, fix it before concluding.

---

## Implementation heuristics

- **Argument validation:** fail fast with specific error messages. Validate type, length, range. State what was expected, what was received, and which argument failed.
- **Side effects:** isolate I/O from pure computation. Avoid hidden global state.
- **NSE:** avoid non-standard evaluation unless the codebase already uses it for that API.

## Error handling & conditions

| Severity | Function | Use case |
|----------|----------|----------|
| Fatal | `stop()` | Unrecoverable; must halt |
| Recoverable | `warning()` | Something is wrong but execution can continue |
| Diagnostic | `message()` | Progress/info only — never for results |

- Use `tryCatch()` with specific condition classes; never bare `try()`.
- Scope handlers to the narrowest operation possible.
- Match on condition *class*, not message text (messages change across R versions and locales).

## NAMESPACE hygiene

- Never call `library()` or `require()` inside package code (`R/`).
- `@importFrom pkg fun` for frequent calls; `pkg::fun()` for one-offs. No blanket `@import pkg`.
- After any NAMESPACE change: run `devtools::document()` and verify. Never hand-edit `NAMESPACE`.
- When removing a dependency: grep `R/` and `tests/` for stale `::` calls and `@importFrom` directives.

## S3 method discipline

- Name methods `generic.class`; avoid dots in class names (dispatch ambiguity).
- Always register via `@export` or explicit `S3method()`.
- Method formals must match the generic exactly (including `...`).
- Assign class in the constructor (`structure()` or `class<-`); provide at minimum a `print` method for user-facing classes.

## Numerical safety

- Make assumptions explicit when thresholds, scaling, or aggregation are involved.
- Division by zero / undefined baselines: choose conservative, documented behavior.
- Track units/dimensions when applicable; enforce dimensional consistency.

## Performance checklist

- Avoid quadratic loops on time-series unless unavoidable; justify if O(n²).
- Preallocate vectors/matrices in loops; no repeated `rbind`/`cbind` inside loops.
- Prefer `vapply()` over `sapply()` when output shape matters; use `for` loops when clarity wins.
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

---

## Reference

If present in the repo, consult: `references/r-workflow.md`.

