# Find synthetic years with at least a reference storm's impact level

Returns every simulation year in which the site's annual impact metric
equalled or exceeded the level caused by a named reference storm,
regardless of which storm(s) were responsible.

The reference impact threshold is resolved in the following order of
priority:

1.  An explicit `threshold` value supplied by the caller.

2.  For `metric = "peak_wind_kt"`: the reference storm's site-level peak
    wind taken directly from `out$events`.

3.  For all metrics (and as a fallback): the median annual metric across
    the simulation years in which that storm appeared in `daily`, as
    identified by
    [`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md).

## Usage

``` r
query_impact_years(
  daily,
  storm_id,
  out = NULL,
  location = NULL,
  metric = c("peak_wind_kt", "cum_damage", "max_damage_rate"),
  threshold = NULL,
  percentile = NULL,
  min_threshold = NULL
)
```

## Arguments

- daily:

  Named list of tibbles returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md),
  or a single tibble.

- storm_id:

  Character scalar IBTrACS SID of the reference storm (e.g.
  `"2017242N16333"` for Hurricane Irma — IBTrACS native SID, not ATCF
  "AL112017").

- out:

  Optional list returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).
  Used to look up the reference storm's historical peak wind when
  `metric = "peak_wind_kt"` and `threshold = NULL`.

- location:

  Character scalar or vector of location names to query, or `NULL` for
  all locations in `daily`.

- metric:

  Character scalar defining how annual impact is measured.

  `"peak_wind_kt"`

  :   Maximum daily sustained wind in the year (kt). Default.

  `"cum_damage"`

  :   Cumulative annual damage fraction (final value of the `cum_damage`
      column).

  `"max_damage_rate"`

  :   Maximum daily damage rate in the year.

- threshold:

  Optional numeric scalar or named numeric vector. When a scalar is
  supplied it is applied to all locations. When a named vector is
  supplied the names must match the queried locations. Overrides
  automatic storm-based threshold derivation; used as the lower bound
  when `percentile` is also specified.

- percentile:

  Optional numeric scalar in `(0, 1)`. When supplied, the empirical
  `percentile` of the annual metric distribution is computed per
  location and used as an additional selection gate. The effective
  threshold becomes `max(ref_threshold, percentile_threshold)`, so only
  years that clear *both* the reference-storm level and the
  distributional percentile are returned. For example
  `percentile = 0.95` retains at most the top 5% of years. Default
  `NULL` disables the percentile gate.

- min_threshold:

  Optional numeric scalar applied as an absolute floor on the effective
  threshold after all other computations. Useful for ensuring a minimum
  physical severity level regardless of the reference storm or
  percentile gate (e.g. `min_threshold = 64` for
  `metric = "peak_wind_kt"` to require at least Cat-1 strength). Default
  `NULL`.

## Value

Tibble with columns:

- `location`:

  Target location name.

- `sim_year`:

  Simulation year index.

- `annual_metric`:

  Annual metric value for that year (column name matches the chosen
  `metric`).

- `threshold`:

  The effective threshold applied for that location, reflecting the
  combined result of storm-based, percentile, and `min_threshold` logic.

Rows are sorted by `location` then `sim_year`.

## See also

[`query_storm_track_years`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md),
[`generate_daily_hazard_impact_spatial`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Years where Saba experienced at least Irma-level peak wind
impact_years <- query_impact_years(
  daily    = daily_out,
  storm_id = "2017242N16333",
  out      = out,
  location = "Saba",
  metric   = "peak_wind_kt"
)
impact_years

# Hybrid: Irma wind as lower bound, 95th-percentile as selection gate
query_impact_years(
  daily      = daily_out,
  storm_id   = "2017242N16333",
  out        = out,
  metric     = "peak_wind_kt",
  percentile = 0.95
)

# Hybrid with explicit Cat-2 floor (83 kt) regardless of Irma's site wind
query_impact_years(
  daily         = daily_out,
  storm_id      = "2017242N16333",
  out           = out,
  metric        = "peak_wind_kt",
  percentile    = 0.95,
  min_threshold = 83
)

# Cumulative damage at least as large as Irma's, with explicit threshold
query_impact_years(
  daily     = daily_out,
  storm_id  = "2017242N16333",
  location  = "Saba",
  metric    = "cum_damage",
  threshold = 0.15
)
} # }
```
