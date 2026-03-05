# Goal
Implement minimal/full output-level gating in `R/hazard_validation.R` so default validation output is compact and deterministic, while preserving full-mode behavior.

# Scope
- Modified only: `R/hazard_validation.R`
- No API signature changes to exported functions.
- Added internal `cfg$output$level` contract and runtime gating.

# Summary
- Added `cfg$output$level` to `make_validation_cfg()` with default `"minimal"`.
- Added internal validator `.validation_output_level()` with allowed values `minimal|full`.
- Enforced minimal-mode bootstrap disable (`n_boot = 0`) in `run_validation_suite()` regardless of `advanced$n_boot`.
- Gated Tier 1B blocked CV so `.run_tier1b_modern_blocked_cv()` is only called in full mode.
- Added minimal diagnostics table `.minimal_diagnostics_from_hindcast()` with exactly:
  - `delta_top1_p50`
  - `delta_overall_p99`
  per `location x storm_class` (`storm_class = "TS34plus"`).
- Wired minimal diagnostics into:
  - `val$summary`
  - `val$bias_diagnostics`
  - `validation_summary.csv`
- Minimal-mode artifacts now limited to:
  - CSV: `validation_summary.csv`, `hindcast_return_levels.csv`, `rate_check.csv`, `wind_field.csv`
  - PNG: `hindcast_return_levels.png`, `plot_bias_diagnostics.png`
- Minimal mode no longer writes markdown tables.
- `plot_hindcast_validation()` now avoids per-location extra plots in minimal mode.
- `plot_bias_diagnostics()` now has a minimal-mode path producing only `plot_bias_diagnostics.png` from the two required diagnostics.
- `validate_hazard_model()` now skips A-D diagnostics execution unless `output$level == "full"`.

# Files Changed
- `R/hazard_validation.R`

# Commands Run
- `rg -n "make_validation_cfg|output\$|validation_summary|hindcast_return_levels|rate_check|wind_field|plot_bias_diagnostics|validation_tables|tier1b|blocked_cv|bootstrap|n_boot|plot_cfg" R/hazard_validation.R`
- `Get-Content DESCRIPTION`
- `Get-Content NAMESPACE`
- `Get-Content .agents/skills/r-coding/SKILL.md`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `git status --short`
- `git diff -- R/hazard_validation.R`

# Test Results
- Parse check passed:
  - `Rscript -e "parse(file='R/hazard_validation.R')"`
- Full validation suite runtime tests were not executed in this pass.

# Behavior Changes
- Default (`make_validation_cfg()`): now behaves as minimal-output mode via `cfg$output$level = "minimal"`.
- Minimal mode:
  - Forces no bootstrap CI in hindcast comparisons.
  - Skips Tier 1B modern blocked CV.
  - Saves only canonical 4 CSV + 2 PNG outputs.
  - No markdown tables.
- Full mode (`cfg$output$level = "full"`): retains existing broader diagnostics/plots/tables paths.

# Follow-ups / Risks
- Existing pre-existing edits in `R/hazard_validation.R` (outside this task) are present in the working tree; this patch was layered on top.
- Recommend one end-to-end smoke run in both modes to confirm exact output file lists against acceptance criteria.
