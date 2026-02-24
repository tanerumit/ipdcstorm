---
name: r-coding
description: Implement, debug, test, and refactor R code and R packages. Use for R functions, failing tests, performance issues, roxygen2 documentation, dependency hygiene, and package checks.
---

# R Coding Skill

This skill provides **R-specific implementation tactics**. Repo-wide rules (API stability, minimal-diff policy, reporting, required validation, encoding, etc.) are defined in `AGENTS.md` and must be followed.

## When to use

- Implement or refactor R functions (base R / tidyverse).
- Diagnose and fix failing `testthat` tests.
- Improve performance (vectorization, preallocation, algorithmic complexity).
- Update roxygen2 docs to match behavior and regenerate Rd/NAMESPACE as needed.
- Resolve package-check issues (imports, NOTES/WARNINGS, Rd problems).

## Working approach

- Start from a minimal reproduction (single function call or single test file).
- Read the smallest set of relevant files first (implementation + tests + package metadata).
- Prefer the narrowest safe fix; avoid broad rewrites unless requested.
- Keep changes locally consistent with existing conventions in the touched files.

## Implementation heuristics

- **Argument validation:** fail fast; validate type/length/range; keep error messages specific.
- **Side effects:** isolate I/O; keep pure computation testable; avoid hidden global state.
- **State & randomness:** make seeds explicit where used; avoid implicit randomness.
- **NSE:** avoid non-standard evaluation unless the codebase already uses it for that API.

## Numerical and data-edge safety

- Handle `NA`/`NaN`/`Inf` intentionally (don’t rely on accidental propagation).
- Make assumptions explicit when thresholds, scaling, or aggregation are involved.
- For division by zero / undefined baselines: choose conservative, documented behavior.

## Performance checklist

- Avoid quadratic loops on time-series unless unavoidable; justify if O(n²).
- Preallocate vectors/matrices in loops; avoid repeated `rbind/cbind` in loops.
- Prefer `vapply()` over `sapply()` when shapes matter; use `for` loops when clarity wins.
- Avoid unnecessary copies of large objects; subset once, reuse.

## Roxygen2 guidance (exported functions)

When updating documentation:
- `@description`: one precise paragraph.
- `@param`: type, length/dimensions, units, NA handling.
- `@return`: exact structure, classes, names, units, and invariants.
- `@examples`: minimal and fast; wrap expensive examples in `\dontrun{}`.

## Testing guidance

- Add/adjust tests when behavior or assumptions change.
- Prefer deterministic tests; set seed explicitly if randomness is involved.
- Cover: input validation, NA/empty/constant cases, return structure (class/names/dims), and key edge conditions.
- Group tests by function or function family; keep fixtures small and local.

## Reference
If present in the repo, consult: `references/r-workflow.md`.