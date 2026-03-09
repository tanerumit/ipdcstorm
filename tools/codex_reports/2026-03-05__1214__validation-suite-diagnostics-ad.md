# Validation Suite Diagnostics A-D Patch

## Goal
Ensure the `DIAGNOSTICS A–D` console block (including the `B)` line with `mult_obs`/`mult_climo`) is emitted during the validation suite workflow (`run_hazard_model()` + `run_validation_suite()`), not only via `validate_hazard_model()`.

## Scope
- Modified: `R/hazard_validation.R`
- No exported API/signature/return-shape changes.
- No dependency changes.
- No test file changes.

## Summary
- Added a guarded diagnostics block near the start of `run_validation_suite()`:
  - Reads `output_level <- .validation_output_level(cfg)` (existing variable).
  - Runs diagnostics only when `output_level == "full"`.
  - Resolves diagnostics targets from:
    1) `out$targets` (preferred)
    2) `cfg$targets` (fallback)
  - If targets are unavailable, prints the standard diagnostics header and a conservative skip reason.
  - Wraps diagnostics execution in `tryCatch` and prints the standard header + `skipped: <error>` on failure.
- Removed the pre-suite diagnostics invocation from `validate_hazard_model()` to avoid duplicate `DIAGNOSTICS A–D` blocks (suite now owns diagnostics printing).

## Files Changed
- `R/hazard_validation.R`

## Commands Run
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "devtools::test()"`

## Test Results
- Parse: success.
- Test suite: success.
- `devtools::test()` summary: `FAIL 0 | WARN 0 | SKIP 0 | PASS 51`.

## Behavior Changes
- `run_validation_suite(out, cfg)` now prints `DIAGNOSTICS A–D` in `full` output mode before tier output when targets are available from `out`/`cfg`.
- If targets are missing, it prints a clear skip reason:
  - `targets not available in run_validation_suite(); diagnostics A–D require targets (lat/lon)...`
- Non-`full` modes remain unchanged (diagnostics not printed).
- `validate_hazard_model()` no longer prints diagnostics separately; diagnostics are printed via the suite call, avoiding duplicate output.

## Follow-ups / Risks
- No automated console-output assertion added; existing tests do not currently include an explicit snapshot/capture pattern for this exact diagnostics banner in `run_validation_suite()`.
- Residual risk is limited to message formatting regressions rather than computational outputs.
