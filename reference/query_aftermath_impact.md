# Measure post-event impact in the days following a primary storm

For each candidate `(location, sim_year)`, locates the end of the
primary event, then summarises the `window_days` days that follow. This
captures compound risk: recovery-period storms, re-intensification, or
successive hits that arrive while the site is already damaged.

**Primary-event anchor:**

- With `storm_id`:

  The last calendar day on which an event whose `event_id` starts with
  `"<storm_id>_y"` is active. Simulation years where that storm was not
  sampled are silently dropped.

- Without `storm_id` (`NULL`):

  The last day of the event that produced the year's peak wind
  (`max(wind_kt)`). Calm years (no TC activity) anchor on the day of
  peak wind even when `event_id` is `NA`.

The aftermath window runs from `anchor_date + 1` through
`anchor_date + window_days`, clipped to the end of the simulated year.

**Aftermath damage:** `aftermath_cum_damage` is the additional
cumulative damage fraction accrued during the window, computed as
\\\max(\texttt{cum\\damage}) -
\texttt{cum\\damage}\[\texttt{anchor\\date}\]\\. Because `cum_damage` is
monotonically non-decreasing within a year, this delta is always \\\ge
0\\.

## Usage

``` r
query_aftermath_impact(
  daily,
  sim_years = NULL,
  storm_id = NULL,
  window_days = 60L,
  location = NULL
)
```

## Arguments

- daily:

  Named list of tibbles returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md),
  or a single tibble.

- sim_years:

  Optional filter applied before computing aftermath. Either:

  - An integer vector of `sim_year` values applied to all locations.

  - A tibble with `location` and `sim_year` columns (the direct output
    of
    [`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)
    or
    [`query_impact_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md)).

  - `NULL` (default) — all years are included.

- storm_id:

  Character scalar IBTrACS SID used to anchor the aftermath window (e.g.
  `"2017242N16333"` for Irma). When `NULL` the window is anchored on the
  year's peak-wind event.

- window_days:

  Positive integer: length of the aftermath window in calendar days
  (default `60L`).

- location:

  Character scalar or vector of location names, or `NULL` for all
  locations.

## Value

Tibble with one row per `(location, sim_year)` and columns:

- `location`:

  Location name.

- `sim_year`:

  Simulation-year index.

- `primary_event_id`:

  The `event_id` of the anchoring event (`NA` for calm years when
  `storm_id = NULL`).

- `event_end_date`:

  Last calendar date of the primary event; the aftermath window starts
  the following day.

- `aftermath_peak_wind_kt`:

  Maximum sustained wind (kt) in the window.

- `aftermath_n_events`:

  Count of distinct storm events active during the window (excluding the
  primary event).

- `aftermath_n_hur_days`:

  Days with sustained wind \\\ge\\ 64 kt during the window.

- `aftermath_cum_damage`:

  Cumulative damage fraction accrued during the window (delta from
  primary event end).

- `aftermath_max_damage_rate`:

  Maximum single-day damage rate during the window.

- `aftermath_rank`:

  Integer rank within location; 1 = worst aftermath by
  `aftermath_cum_damage`.

Rows are sorted by `location` then `aftermath_rank`.

## See also

[`query_storm_track_years`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md),
[`query_impact_years`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md),
[`compute_stress_year_metrics`](https://tanerumit.github.io/ipdcstorm/reference/compute_stress_year_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Aftermath in the 60 days following an Irma-track event
track_years  <- query_storm_track_years(daily_out, storm_id = "2017242N16333")
aftermath_60 <- query_aftermath_impact(
  daily       = daily_out,
  sim_years   = track_years,
  storm_id    = "2017242N16333",
  window_days = 60L
)

# 90-day aftermath anchored on each year's own peak-wind event
impact_years  <- query_impact_years(daily_out, storm_id = "2017242N16333",
                                     out = hazard_out)
aftermath_90  <- query_aftermath_impact(
  daily       = daily_out,
  sim_years   = impact_years,
  storm_id    = NULL,
  window_days = 90L
)
} # }
```
