# Find synthetic years where a specific storm's track was sampled

Scans the daily hazard-impact output and returns the simulation years in
which a named historical storm was drawn from the event library,
preserving its original approach geometry, temporal wind profile, and
intensity.

The event-library sampler encodes the source storm in the `event_id`
column as `"<SID>_y<year>_<counter>"`. This function matches rows where
`event_id` begins with `"<storm_id>_y"`, so any simulation year that
received that track will be returned.

## Usage

``` r
query_storm_track_years(daily, storm_id, location = NULL)
```

## Arguments

- daily:

  Named list of tibbles returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md),
  or a single tibble for one location.

- storm_id:

  Character scalar IBTrACS SID of the target storm (e.g.
  `"2017242N16333"` for Hurricane Irma — IBTrACS native SID, not ATCF
  "AL112017").

- location:

  Character scalar or vector of location names to query, or `NULL` to
  query all locations present in `daily`.

## Value

Tibble with columns:

- `location`:

  Target location name.

- `sim_year`:

  Synthetic simulation-year index in which the storm was sampled.

Rows are sorted by `location` then `sim_year`. An empty tibble is
returned (with the same schema) when the storm never appears.

## See also

[`query_impact_years`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md),
[`generate_daily_hazard_impact_spatial`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Find all synthetic years where Hurricane Irma's track was resampled
track_years <- query_storm_track_years(
  daily    = daily_out,
  storm_id = "2017242N16333"
)
track_years

# Restrict to one site
query_storm_track_years(daily_out, storm_id = "2017242N16333", location = "Saba")
} # }
```
