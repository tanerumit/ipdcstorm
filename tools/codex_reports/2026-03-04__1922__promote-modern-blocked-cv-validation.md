Goal
- Promote modern-era blocked CV diagnostics to a Tier 1 companion in validation output, fix duplicate blocked-CV folds in `dev/validation_stats.R`, and add a compact per-site manuscript/report table.

Scope
- Limited to the validation and diagnostics layer.
- No hazard physics, simulation, exported function signatures, or top-level `run_validation_suite()` return fields were changed.
- No new dependencies were added.

Summary
- Removed the extra post-loop append in `dev/validation_stats.R` so blocked CV no longer duplicates the last fold.
- Added internal Tier 1B blocked-CV helpers in `R/hazard_validation.R` to:
  - derive annual maxima from existing per-site trackpoints
  - run blocked CV with defaults (`1970+`, TS34+, 10-year blocks)
  - filter scored folds to `test_n_pos_years >= 2`
  - print and return per-site summary metrics
  - build a compact manuscript table and optionally save it as CSV when `save_tables = TRUE`
- Integrated a new console section in `run_validation_suite()`:
  - `TIER 1A: HINDCAST VALIDATION (single-split RL-in-CI; secondary)`
  - `TIER 1B: MODERN BLOCKED CV (annual maxima, TS+)`
- Added Tier 1B tables to saved validation table artifacts when `save_tables = TRUE`.

Files Changed
- `R/hazard_validation.R`
- `dev/validation_stats.R`

Commands Run
- `Get-Content -Path .agents/skills/r-coding/SKILL.md`
- `Get-Content -Path DESCRIPTION`
- `Get-Content -Path NAMESPACE`
- `Get-Content -Path dev\validation_stats.R`
- `Get-Content -Path R\hazard_validation.R`
- `rg -n "run_validation_suite|TIER|hindcast RL|single-split|save_tables|validation" R\hazard_validation.R`
- `rg -n "validation_stats|run_validation_suite|hazard_validation|blocked_cv|fold_id" tests\testthat`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='dev/validation_stats.R')"`
- `Rscript -e "devtools::test()"`
- Synthetic validation for duplicate blocked-CV folds in `dev/validation_stats.R`
- Synthetic validation for `.run_tier1b_modern_blocked_cv()` and `.modern_blocked_cv_paper_table()`

Test Results
- `Rscript -e "parse(file='R/hazard_validation.R')"`: passed
- `Rscript -e "parse(file='dev/validation_stats.R')"`: passed
- `Rscript -e "devtools::test()"`: passed (`PASS 31`, `FAIL 0`, `WARN 0`, `SKIP 0`)
- Synthetic blocked-CV check in `dev/validation_stats.R`: passed (`nrow(res$summary) == 3`, unique `site + fold_id`)
- Synthetic Tier 1B helper check: passed; printed:
  - per-site summary columns: `site`, `folds_total`, `folds_used`, `mean_test_p0`, `mean_test_q90`, `max_test_q95`, `era_min_year`, `storm_vmax_min`, `block_size`, `threshold_kt`
  - compact table columns: `site`, `era`, `folds_used_total`, `p0_mean`, `q90_mean`, `q95_min`, `q95_max`, `q90_min`, `q90_max`

Behavior Changes
- `run_validation_suite()` now prints Tier 1 hindcast explicitly as a single-split secondary check and adds a new Tier 1B blocked-CV section.
- When `save_tables = TRUE`, the validation workflow now also saves:
  - `modern_blocked_cv_summary.csv`
  - `modern_blocked_cv_compact.csv`
  - markdown table sections for the two Tier 1B summaries
- `dev/validation_stats.R` blocked-CV output now yields exactly one summary row per site-fold instead of duplicating the last fold.

Follow-ups / Risks
- Tier 1B uses existing trackpoint columns and fills a few optional diagnostics columns with `NA` when absent; this keeps the summary deterministic but assumes `iso_time`, site wind, and storm wind columns are present.
- The compact table uses `era = "1970+"` for stability; if a closed year range is preferred later, it can be derived from the actual annual-max record without changing the current scoring logic.
