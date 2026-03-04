# Codex Report

## Goal

Upgrade `.run_validation_diagnostics_tail_geometry()` to use a more stable tail definition and add `r/RMW` plus event-max diagnostics.

## Scope

- Modified only the internal helper in `R/validation_diag_tail_geometry.R`.
- No exported function or documented public API changed.

## Summary

- Added optional arguments: `tail_prob`, `top_k`, `tail_mode`, and `storm_filter`.
- Replaced the old q99-only single-row tail logic with configurable quantile / top-k / union selection.
- Added `r/RMW` tail metrics to the returned tibble.
- Added per-site event-max printing (`SID` grouped maxima, top 10 by `V_site_kt`).
- Kept the helper deterministic and dependency-free.

## Files Changed

- `R/validation_diag_tail_geometry.R`

## Commands Run

- `Rscript -e "parse(file='R/validation_diag_tail_geometry.R')"`
- `Rscript -e "pkgload::load_all('.'); testthat::test_file('tests/testthat/test-rmw_estimation.R')"`
- `@' ... '@ | Rscript -` with:
  - `pkgload::load_all('.')`
  - `cfg <- make_hazard_cfg(data_path = 'inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv')`
  - `targets <- tibble::tribble(...)`
  - `out <- run_hazard_model(cfg = cfg, targets = targets, verbose = FALSE)`
  - `.run_validation_diagnostics_tail_geometry(out$trackpoints, targets, sites = c('Saba', 'Statia'))`

## Test Results

- Parse check passed for the touched R file.
- The existing `tests/testthat/test-rmw_estimation.R` subset passed after `pkgload::load_all('.')`.
- A direct smoke run of the updated helper on `Saba` and `Statia` completed successfully and printed the expected summary plus event-max sections.

## Behavior Changes

- Internal helper now defaults to `tail_mode = "both"` and `storm_filter = "storm_vmax"`.
- Returned tibble now includes `frac_rnorm_lt1`, `rnorm_p10`, `rnorm_p50`, and `rnorm_p90`.
- Console title now reads `TAIL GEOMETRY DIAGNOSTICS (TS+ tail)`.

## Follow-Ups / Risks

- The event-max print uses `message()` for headers and `print()` for tables, so stdout/stderr interleaving can make the table appear before the header line in some terminals.
- If exact console ordering becomes important, the next step would be to standardize all diagnostic output onto a single stream.
