---
name: r-coding
description: Implement, debug, test, and refactor R code and R packages. Use for R functions, failing tests, performance issues, roxygen2 documentation, dependency hygiene, and package checks.
---

# R Coding Skill

This skill is for implementing, debugging, testing, and refactoring R code in research, analytical, and package-development contexts.

Repo-specific rules such as project layout, scientific domain constraints, reporting format, and task-specific conventions live in `AGENTS.md` and always apply on top of this skill.

No automated linting is assumed unless the repo says otherwise. Stay locally consistent with the touched file. For new files, use tidyverse style guide basics: snake_case, 2-space indent, and `<-` for assignment. Do not reformat untouched code.

## Operating mode

Default to **prototype-first** design unless the user explicitly requests an **API-focused** approach.

Prototype-first means:
- optimize for clarity, directness, inspectability, and speed of iteration
- prefer flat inputs, simple function signatures, and transparent control flow
- prefer direct parameter passing over nested config objects
- bias toward deleting, flattening, and merging before adding abstractions
- keep only the minimum structure needed for the current workflow
- avoid generalized APIs, wrapper layers, compatibility shims, factories, registries, and abstraction for hypothetical future reuse

Switch to an API-focused or production-oriented approach only when the user explicitly asks for it, or when the task clearly concerns stable public interfaces, packaging, deployment, release, multi-user reuse, or interface stability.

Do not infer CRAN-grade package standards unless the task explicitly concerns release, checks, documentation, exports, or interface stability.

When choosing between a prototype-oriented solution and a production-oriented solution, choose the prototype-oriented solution by default unless the user explicitly asks otherwise.

Prototype-first affects architecture choice, not scientific rigor, numerical care, namespace discipline, or error hygiene.

## When to use

- Implement or refactor R functions.
- Diagnose and fix failing `testthat` tests.
- Improve performance through vectorization, preallocation, or algorithmic simplification.
- Update `roxygen2` docs and regenerate `Rd`/`NAMESPACE` where relevant.
- Resolve package-check issues such as imports, `NOTES`/`WARNINGS`, or `Rd` problems.
- Simplify over-engineered interfaces or unnecessary abstraction in current workflows.

---

## HARD CONSTRAINTS (NON-NEGOTIABLE)

These apply to every change unless the user explicitly overrides them.

### Change discipline
- **Minimal diff.** Touch only what is necessary. No style-only reformats, renames, or file moves.
- **Public API stability by default.** Do not change exported function names, arguments, defaults, return types, column names, or classes unless explicitly instructed.
- **Prototype-first exception.** If the user explicitly says backward compatibility is not needed, remove compatibility layers, legacy wrappers, and unnecessary API-preservation scaffolding rather than preserving them.
- **No silent behavior changes.** Any change to user-visible outputs, warnings, errors, or messages must be documented in the final summary.

### Code style
- **No purrr / functional-programming rewrites.** Use explicit `for` loops and base R idioms unless the touched file already uses a different convention locally.
- **UTF-8 only** (no BOM).

### Dependencies
- **No new dependencies** unless necessary and justified.
- Do not require `Suggests` packages at runtime.
- If a new dependency is unavoidable, guard with `requireNamespace()` and justify it in the summary.

### Determinism
- Preserve deterministic behavior.
- Control randomness with explicit seeds.
- Use fixed seeds for stochastic code unless the user explicitly requests random behavior.

### Data integrity
- Do not silently drop columns, fields, or rows.
- Handle `NA`, `NaN`, and `Inf` intentionally. Never rely on accidental propagation.

### Scientific and numerical correctness
- Do not silently simplify equations, logic, or scientific assumptions.
- Make assumptions explicit when thresholds, scaling, normalization, or aggregation are involved.
- Track units and dimensions where relevant.
- Prefer numerically stable formulations over cleverness.
- If something is unknown from the available code or context, say so.

---

## Required workflow

1. **Read first.** Open the narrowest relevant set of files: implementation, tests, `DESCRIPTION`, `NAMESPACE`, and docs if relevant.
2. **Determine the mode.** Decide whether the task is prototype-first or explicitly API-focused.
3. **Inspect before rewriting.** Identify which parts contain substantive logic and which parts are wrappers, restructuring layers, validation plumbing, or dead weight.
4. **Smallest safe fix.** Prefer the narrowest change that solves the problem. Avoid broad rewrites unless requested.
5. **Simplify first.** Before adding a new layer, ask:
   - can something be deleted?
   - can nested inputs be flattened?
   - can thin wrappers be merged or removed?
   - can direct parameter passing replace a config object or transformation layer?
6. **Stay locally consistent** with conventions in the touched files.
7. **Update co-artifacts.** If behavior or assumptions change, update tests. If exported behavior or docs change, update `roxygen`.
8. **Validate before concluding** using the minimum required validation below. Do not mark the task done until validation passes or was explicitly skipped with reason.

---

## Required validation (minimum)

Run these for every change and report results.

1. Parse check — every touched R file
   Rscript -e "parse(file='R/<touched-file>.R')"

2. Unit tests — if a matching test file exists
   Rscript -e "testthat::test_file('tests/testthat/test_<touched>.R')"

3. Roxygen — if any documentation changed
   Rscript -e "devtools::document()"
   Then verify NAMESPACE and Rd files are consistent.

If a validation step fails, fix it before concluding, unless the user explicitly asked for analysis only.

---

## Implementation heuristics

- **Argument validation:** fail fast with specific error messages. Validate type, length, range, and structure where relevant. State what was expected, what was received, and which argument failed.
- **Side effects:** isolate I/O from pure computation. Avoid hidden global state.
- **NSE:** avoid non-standard evaluation unless the codebase already uses it for that API.
- **Code structure:** prefer fewer, clearer functions over many thin wrappers.
- **Abstractions:** every new helper, abstraction layer, or object structure must solve a current problem, not a hypothetical future one.
- **Data flow:** keep control flow and argument plumbing easy to trace.

---

## Error handling & conditions

| Severity | Function | Use case |
|----------|----------|----------|
| Fatal | `stop()` | Unrecoverable; execution must halt |
| Recoverable | `warning()` | Something is wrong but execution can continue |
| Diagnostic | `message()` | Progress/info only; never for results |

- Use `tryCatch()` with specific condition classes. Never use bare `try()` as the default error-handling pattern.
- Scope handlers to the narrowest operation possible.
- Match on condition **class**, not message text.
- Do not add elaborate condition frameworks unless the current task requires them.

---

## NAMESPACE hygiene

- Never call `library()` or `require()` inside package code under `R/`.
- Use `@importFrom pkg fun` for repeated calls; use `pkg::fun()` for one-offs.
- Do not use blanket `@import pkg`.
- After any `NAMESPACE`-relevant change, run `devtools::document()` and verify the result.
- Never hand-edit `NAMESPACE`.
- When removing a dependency, search `R/` and `tests/` for stale `::` calls and `@importFrom` directives.

---

## S3 method discipline

- Name methods `generic.class`.
- Avoid dots in class names where they create dispatch ambiguity.
- Register methods via `@export` or explicit `S3method()`.
- Method formals must match the generic exactly, including `...` where required.
- Assign class in the constructor using `structure()` or `class<-`.
- For user-facing classes, provide at minimum a `print` method when appropriate.

---

## Numerical safety

- Make assumptions explicit when thresholds, scaling, aggregation, interpolation, or normalization are involved.
- For division by zero, undefined baselines, or empty denominators, choose conservative, documented behavior.
- Track units and dimensions where applicable; enforce dimensional consistency.
- Handle `NA`, `NaN`, and `Inf` intentionally.
- When numerical behavior changes, state the reason and expected effect.

---

## Performance checklist

- Avoid quadratic loops on time-series or large vectors unless unavoidable; justify if complexity is inherently `O(n²)`.
- Preallocate vectors and matrices in loops.
- Do not repeatedly `rbind()` or `cbind()` inside loops when avoidable.
- Prefer `vapply()` over `sapply()` when output shape matters.
- Use `for` loops when they improve clarity and control.
- Subset large objects once and reuse them to avoid unnecessary copies.

---

## Roxygen2 guidance (exported functions)

- `@description`: one precise paragraph.
- `@param`: include type, length or dimensions, units, and `NA` handling where relevant.
- `@return`: state exact structure, classes, names, units, and invariants.
- `@examples`: keep minimal and fast. Wrap expensive examples in `\dontrun{}`.
- If exported behavior changes, update `roxygen` and ensure generated `Rd` and `NAMESPACE` stay consistent.

---

## Testing guidance

- Add or adjust tests when behavior or assumptions change.
- Use explicit seeds for stochastic code.
- Prefer deterministic tests.
- Cover:
  - input validation
  - `NA` / empty / constant cases
  - return structure (`class`, names, dimensions)
  - key edge conditions
- Group tests by function or function family.
- Keep fixtures small and local.

---

## Review criteria

When reviewing or refactoring code, evaluate against:
- correctness
- numerical stability
- scientific transparency
- unnecessary abstraction
- hidden assumptions
- API complexity relative to current task needs
- maintainability for the current project stage

Flag:
- accidental complexity
- dead layers
- weak naming
- mismatched abstractions
- fragile argument plumbing
- unnecessary backward-compatibility burden
- unverified assumptions

For refactor proposals, default to the simplest scientifically adequate design, not the most software-engineered one.

---

## Reference

If present in the repo, consult: `references/r-workflow.md`.