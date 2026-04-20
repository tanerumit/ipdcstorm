# Compute per-year multi-metric summary for stress-test characterisation

Aggregates the daily hazard-impact series to one row per
`(location, sim_year)`, capturing seven complementary properties of each
simulated year:

- `peak_wind_kt`:

  Maximum daily sustained wind (kt). Primary intensity indicator.

- `n_ts_days`:

  Days with sustained wind \\\ge\\ 34 kt. Captures tropical-storm-level
  duration exposure.

- `n_hur_days`:

  Days with sustained wind \\\ge\\ 64 kt. Captures hurricane-level
  duration exposure.

- `n_events`:

  Count of distinct storm events in the year. Proxy for compound-event
  risk; years with multiple events have compounding recovery and damage
  burdens.

- `max_event_dur_days`:

  Duration in days of the longest single event (number of days
  attributed to one `event_id`). Identifies slow-moving or stalling
  storms.

- `cum_damage`:

  Total cumulative damage fraction over the year (`max(cum_damage)`
  across all days). Aggregate impact measure.

- `max_damage_rate`:

  Maximum single-day damage rate. Captures acute shock intensity
  independently of duration.

## Usage

``` r
compute_stress_year_metrics(daily, sim_years = NULL, location = NULL)
```

## Arguments

- daily:

  Named list of tibbles returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md),
  or a single tibble.

- sim_years:

  Optional filter. Either:

  - An integer vector of `sim_year` values applied to all locations.

  - A tibble with `location` and `sim_year` columns (the direct output
    of
    [`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)
    or
    [`query_impact_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md))
    for per-location filtering.

  - `NULL` (default) to include all years.

- location:

  Character vector of location names to include, or `NULL` for all.

## Value

Tibble with columns `location`, `sim_year`, and the seven metrics listed
above. Metric values are zero for calm years (no TC activity).

## See also

[`aggregate_stress_metrics`](https://tanerumit.github.io/ipdcstorm/reference/aggregate_stress_metrics.md),
[`select_stress_years`](https://tanerumit.github.io/ipdcstorm/reference/select_stress_years.md),
[`query_impact_years`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md),
[`query_storm_track_years`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute metrics for all impact-based candidate years
filtered   <- query_impact_years(daily_out, storm_id = "AL112017", out = hazard_out)
metrics    <- compute_stress_year_metrics(daily_out, sim_years = filtered)
} # }
```
