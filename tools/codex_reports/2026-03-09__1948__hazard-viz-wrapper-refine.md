Goal
- Refine selected hazard visualization behavior in `R/hazard_viz.R` and add a wrapper that standardizes labeling/theme sizing and saves exactly five PNG plots, while excluding seasonality from wrapper outputs.

Scope
- Edited plotting logic only in `R/hazard_viz.R`.
- Kept `plot_seasonality_doy()` present and unchanged in wrapper coverage.
- Added no new dependencies.
- Did not regenerate roxygen outputs because the task explicitly constrained edits to `R/hazard_viz.R`.

Summary
- Added internal helpers to standardize plot labels/theme sizing and save PNG files with `ggsave()`.
- Changed `plot_monthly_events()` from grouped bars to zero-filled year-month-class boxplots with explicit mean markers.
- Changed `plot_annual_counts()` to empirical probability mass by annual count, with Poisson probabilities overlaid for event counts.
- Added deterministic vertical-only jitter to `plot_intensity_duration()`.
- Added exported wrapper `save_hazard_viz_plots()` to build and save five standardized plots.

Files Changed
- `R/hazard_viz.R`

Commands Run
- `Rscript -e "parse(file='R/hazard_viz.R')"`
- `Rscript -e "source('R/hazard_viz.R')"`
- `Rscript -e "devtools::test(filter = 'viz|hazard_viz')"`
- Smoke test via `Rscript -` with explicit `dplyr`, `tidyr`, `ggplot2`, and `magrittr` loads, synthetic daily data, and `save_hazard_viz_plots(...)`.
- `git diff -- R/hazard_viz.R`

Test Results
- Parse: passed.
- Source: passed.
- `devtools::test(filter = 'viz|hazard_viz')`: failed because no matching test files exist (`No test files found.`).
- Smoke test: passed; created exactly five files:
- `monthly_events_boxplot.png`
- `annual_event_count_probability.png`
- `wind_timeseries.png`
- `wind_distribution.png`
- `intensity_duration.png`

Behavior Changes
- `plot_monthly_events()` now shows the across-year distribution of monthly event starts by class rather than aggregated grouped bars.
- `plot_annual_counts()` now uses probability in a given year on the y-axis instead of number of years.
- `plot_intensity_duration()` now applies small deterministic vertical jitter to reduce overlap.
- New wrapper `save_hazard_viz_plots()` saves standardized PNG outputs for five plots and excludes `plot_seasonality_doy()`.

Follow-ups/Risks
- The new wrapper is documented with roxygen in `R/hazard_viz.R`, but `NAMESPACE`/`man/` were not regenerated due the user’s single-file edit constraint. Package export/documentation sync should be done in a follow-up if this function is intended for installed-package use immediately.
- `plot_monthly_events(normalize = TRUE)` is retained for API compatibility, but the new boxplot interpretation is already per-year, so the flag no longer changes the plotted values.
