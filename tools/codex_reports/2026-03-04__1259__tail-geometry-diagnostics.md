# Codex Report

## Goal

Run a focused TS+ tail-geometry diagnostic for Saba and Statia (with St_Martin as a control), checking radius quantization and directional-vs-mean geometry consistency.

## Scope

- Added one new internal helper under `R/`.
- Did not change any exported APIs, return structures, or existing workflow side effects.
- Wrote one analysis note under `dev/notes/`.

## Summary

- Added `.run_validation_diagnostics_tail_geometry()` as a non-exported helper in `R/validation_diag_tail_geometry.R`.
- The helper:
  - accepts `out$trackpoints` and `targets`,
  - reruns `compute_site_winds_full()` for requested sites,
  - restricts to TS+ rows with finite positive `R34_km`,
  - computes the requested q99-tail geometry metrics,
  - prints a console section and returns the tibble invisibly.
- Ran the helper once using the packaged IBTrACS file under `inst/extdata/...`.
- Results show quantized `R34_km` values in the TS+ subset, but not a dominant 90 nm (`166.68 km`) median for Saba/Statia in this restricted regime. The directional-vs-mean geometry gap is materially larger for Saba/Statia than for St_Martin.

## Files Changed

- `R/validation_diag_tail_geometry.R`
- `dev/notes/diagnostics_tail_geometry.md`

## Commands Run

- `Rscript -e "parse(file='R/validation_diag_tail_geometry.R')"`
- `@' ... '@ | Rscript -` with:
  - `pkgload::load_all('.')`
  - `cfg <- make_hazard_cfg(data_path = 'inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv')`
  - `targets <- tibble::tribble(...)`
  - `out <- run_hazard_model(cfg = cfg, targets = targets, verbose = FALSE)`
  - `res <- .run_validation_diagnostics_tail_geometry(out$trackpoints, targets)`
- `@' ... '@ | Rscript -` to inspect the top rounded `R34_km` frequencies inside the TS+ subset for `Saba`, `Statia`, and `St_Martin`

## Test Results

- `parse(file='R/validation_diag_tail_geometry.R')` succeeded.
- The helper executed successfully after `pkgload::load_all('.')`.
- No dedicated `testthat` file was added because this is an internal analysis helper and the task asked to prefer minimal/no code changes.

## Behavior Changes

- No existing exported behavior changed.
- New functionality is opt-in only via direct calls to `.run_validation_diagnostics_tail_geometry(...)`.

## Follow-Ups / Risks

- The helper assumes `out_trackpoints[[site]]` contains the per-site tables produced by `run_hazard_model()`. If callers pass a different object shape, it will stop early.
- The quantization itself likely reflects discrete source radii in nautical miles. If this becomes a modeling issue rather than just a diagnostic observation, the next places to inspect are `.get_directional_radius()`, the quadrant radius ingestion path, `mean_radius_nm()`, and the bearing/quadrant logic in `R/hazard_core.R`.
