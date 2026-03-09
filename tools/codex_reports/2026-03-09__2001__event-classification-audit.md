Goal
- Audit the TS/HUR event-classification pipeline and determine whether label/peak mismatches are cosmetic only or affect upstream/downstream calculations.

Scope
- Analysis only.
- No source files edited.
- Audited `R/hazard_core.R`, `R/hazard_downscale.R`, `R/hazard_viz.R`, `R/hazard_validation.R`, and `R/hazard_run.R` as the helper file directly constructing canonical event tables.

Summary
- Canonical `storm_class` for model event tables is assigned from realized peak site wind in `run_hazard_model()` via `classify_severity(peak_wind_kt, 34, 64)` after `make_storm_events()` computes `peak_wind_kt = max(V_site_kt)`.
- Downscaled daily `event_class` is not recomputed from realized daily/event peak. It is propagated from sampled severity in `sample_events_for_year_extended()` and copied into daily rows in `generate_daily_year_extended()`.
- Therefore, core continuous-wind calculations and canonical event-table rate/intensity summaries are not contaminated by the daily-label mismatch path.
- The mismatch affects downstream class-specific summaries/plots built from daily `event_class` plus realized daily/event peak, especially `prep_events()` in `R/hazard_viz.R`.
- Minimal safe fix location is the downscaled daily/event path: recompute or validate `event_class` from the final realized event peak in `sampled_events` or immediately after `generate_daily_year_extended()`, not in visualization code.

Files Changed
- None.

Commands Run
- `rg -n "classify_severity|event_class|storm_class|sampled_events\\$event_class|V_peak|peak_wind_kt|max_wind_kt|V_site_kt|prep_events|make_storm_events|simulate|downscale" R/hazard_core.R R/hazard_downscale.R R/hazard_viz.R R/hazard_validation.R R`
- `Rscript -e "library(dplyr); library(tidyr); library(ggplot2); library(magrittr); source('R/hazard_core.R'); source('R/hazard_downscale.R'); source('R/hazard_viz.R'); source('R/hazard_validation.R')"`
- `Rscript -e "devtools::test(filter = 'hazard|validation|viz|downscale')"`

Test Results
- Source/load command: passed.
- `devtools::test(filter = 'hazard|validation|viz|downscale')`: passed existing matched tests (`hazard_outer_cutoff`, 6 passes).
- Warning: `Objects listed as exports, but not present in namespace: • plot_summary_panel` was emitted during test startup; unrelated to this audit and not modified.

Behavior Changes
- None.

Follow-ups/Risks
- If the package relies on downscaled daily `event_class` for analytical summaries outside `R/hazard_viz.R`, those summaries may be biased by carried labels when realized pulse peaks cross thresholds.
- A later fix should be applied before any daily/event summaries are derived, ideally at sampled-event construction or daily generation time, to keep all downstream consumers consistent.
