# R Workflow Reference

## Fast Command Map

- Run one test file:
  - `Rscript -e "testthat::test_file('tests/testthat/test-foo.R')"`
- Run tests by filter:
  - `Rscript -e "testthat::test_dir('tests/testthat', filter='foo')"`
- Run all tests:
  - `Rscript -e "testthat::test_local('.')"`
- Load package in dev session:
  - `Rscript -e "pkgload::load_all('.')"`
- Regenerate docs:
  - `Rscript -e "devtools::document()"`
- Run package check:
  - `Rscript -e "devtools::check()"`

## Debugging Pattern

1. Reproduce with one failing test.
2. Add temporary prints or `browser()` at the narrowest branch.
3. Inspect classes, lengths, names, and missing values (`str()`, `class()`, `is.na()`).
4. Remove debug scaffolding before finalizing.

## Design Heuristics

- Separate pure data transforms from file/network side effects.
- Validate inputs at function boundaries.
- Keep vectorized paths clear; avoid hidden recycling.
- Use informative errors with `stop()`/`rlang::abort()`.

## Test Heuristics

- Prefer behavior assertions over implementation details.
- Add regression tests for bug fixes.
- Use fixed seeds for stochastic code.
- Keep fixtures minimal and local to the test context.