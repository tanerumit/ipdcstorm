# Improve hazard_climate documentation
- Timestamp: 2026-02-24 17:32 (local)
- Task slug: improve-hazard-climate-docs

## Goal
- Improve the documentation quality for `R/hazard_climate.R` and keep generated Rd docs consistent.

## Scope
- Touched: `R/hazard_climate.R`, `man/make_sst_cfg.Rd`, `man/prepare_sst_data.Rd`, `man/simulate_twolevel_counts.Rd`
- Not touched: Runtime logic/algorithms in climate simulation functions.

## Implementation summary
- Clarified scenario documentation in `make_sst_cfg()` to match actual accepted values from `sst_scenario_info("all")`.
- Added `@details` for `advanced$cc_params` behavior (`NULL`, `list()`, named overrides).
- Documented missing `cc_params` return field in `prepare_sst_data()`.
- Improved formula formatting/readability in `simulate_twolevel_counts()` documentation.
- Regenerated affected Rd files with `devtools::document()`.

## Files changed
- R/hazard_climate.R
- man/make_sst_cfg.Rd
- man/prepare_sst_data.Rd
- man/simulate_twolevel_counts.Rd

## Commands run
```text
Get-ChildItem -Force
rg --files R man | rg hazard_climate
Get-Content/rg inspections for R/hazard_climate.R and related man files
git diff -- R/hazard_climate.R
Rscript -e "devtools::document()"
Rscript -e "files <- list.files('man', pattern='\\.Rd$', full.names = TRUE); for (f in files) tools::checkRd(f); cat('checked', length(files), 'Rd files\\n')"
Get-Date -Format "yyyy-MM-dd HH:mm"
```