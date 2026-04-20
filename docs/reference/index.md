# Package index

## Core hazard workflow

- [`make_hazard_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md)
  : Create a hazard model configuration
- [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
  : Run the stochastic tropical-cyclone hazard model for one or more
  target locations

## Climate scenario tools

- [`make_climate_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_climate_cfg.md)
  : Build a climate configuration object for the hazard model
- [`resolve_climate_inputs()`](https://tanerumit.github.io/ipdcstorm/reference/resolve_climate_inputs.md)
  : Resolve climate scalars for the hazard model

## Validation

- [`make_validation_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_validation_cfg.md)
  : Create a validation configuration
- [`validate_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/validate_hazard_model.md)
  : Run hazard model and validate in one step
- [`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md)
  : Run the three-tier hazard model validation suite

## IBTrACS and storm inputs

- [`ibtracs_demo`](https://tanerumit.github.io/ipdcstorm/reference/ibtracs_demo.md)
  : North Atlantic IBTrACS Demo Data
- [`read_ibtracs_clean()`](https://tanerumit.github.io/ipdcstorm/reference/read_ibtracs_clean.md)
  : Read and clean IBTrACS CSV (USA fields + pressure/structure when
  present)
- [`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)
  : Find synthetic years where a specific storm's track was sampled
- [`lookup_storm_id()`](https://tanerumit.github.io/ipdcstorm/reference/lookup_storm_id.md)
  : Look up IBTrACS storm IDs from the hazard model event record
- [`make_storm_events()`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md)
  : Aggregate track points into storm events.
- [`compute_storm_heading()`](https://tanerumit.github.io/ipdcstorm/reference/compute_storm_heading.md)
  : Compute storm heading by track point.

## Downscaling and Impact

- [`build_event_library()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md)
  : Build an empirical event library for resampling
- [`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
  : Background wind configuration for correlated Weibull generation
- [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)
  : Generate spatially coherent daily synthetic hazard and impact time
  series
- [`damage_rate_from_wind()`](https://tanerumit.github.io/ipdcstorm/reference/damage_rate_from_wind.md)
  : Bounded power-law damage rate from wind speed
- [`save_daily_hazard_csvs()`](https://tanerumit.github.io/ipdcstorm/reference/save_daily_hazard_csvs.md)
  : Save daily hazard output to per-location CSV files

## Query and Stress Testing

- [`query_impact_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md)
  : Find synthetic years with at least a reference storm's impact level
- [`query_aftermath_impact()`](https://tanerumit.github.io/ipdcstorm/reference/query_aftermath_impact.md)
  : Measure post-event impact in the days following a primary storm
- [`compute_stress_year_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md)
  : Compute per-year multi-metric summary for stress-test
  characterisation
- [`aggregate_stress_metrics()`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md)
  : Compute a weighted composite stress score from per-year metrics
- [`select_stress_years()`](https://tanerumit.github.io/ipdcstorm/reference/select_stress_years.md)
  : Select k representative stress-test years with good sample coverage

## Visualization

- [`plot_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/plot_annual_counts.md)
  : Plot distribution of annual event totals
- [`plot_doy_wind()`](https://tanerumit.github.io/ipdcstorm/reference/plot_doy_wind.md)
  : Plot mean and maximum wind by day of year
- [`plot_intensity_duration()`](https://tanerumit.github.io/ipdcstorm/reference/plot_intensity_duration.md)
  : Plot event intensity versus duration
- [`plot_monthly_events()`](https://tanerumit.github.io/ipdcstorm/reference/plot_monthly_events.md)
  : Plot monthly counts of event starts
- [`plot_monthly_quantiles()`](https://tanerumit.github.io/ipdcstorm/reference/plot_monthly_quantiles.md)
  : Plot monthly wind quantile summary
- [`plot_return_levels()`](https://tanerumit.github.io/ipdcstorm/reference/plot_return_levels.md)
  : Plot empirical wind return levels
- [`plot_seasonality_doy()`](https://tanerumit.github.io/ipdcstorm/reference/plot_seasonality_doy.md)
  : Plot day-of-year distribution of event timing
- [`plot_wind_distribution()`](https://tanerumit.github.io/ipdcstorm/reference/plot_wind_distribution.md)
  : Plot distribution of daily wind speed
- [`plot_wind_timeseries()`](https://tanerumit.github.io/ipdcstorm/reference/plot_wind_timeseries.md)
  : Plot daily wind speed with event-duration overlays
- [`save_hazard_viz_plots()`](https://tanerumit.github.io/ipdcstorm/reference/save_hazard_viz_plots.md)
  : Save standard hazard visualization plots as PNG files
