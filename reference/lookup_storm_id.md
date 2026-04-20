# Look up IBTrACS storm IDs from the hazard model event record

The hazard model stores storms under their **IBTrACS native SID**
(format: `YYYYDDDLLLBBB`, e.g. `"2017242N16333"` for Hurricane Irma),
*not* the ATCF identifier familiar from NHC advisories
(`"2017242N16333"`). This function filters `out$events` so you can find
the correct `storm_id` to pass to
[`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)
and
[`query_impact_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md)
using human-readable criteria.

All filter arguments are optional and combined with AND logic. Omitting
all arguments returns the full event record.

## Usage

``` r
lookup_storm_id(
  out,
  year = NULL,
  location = NULL,
  date_range = NULL,
  min_wind_kt = NULL
)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- year:

  Integer scalar or vector: calendar season year(s) to include.

- location:

  Character scalar or vector: location name(s) to filter on. `NULL`
  aggregates across locations (returns distinct storm IDs).

- date_range:

  Length-2 vector of `Date` or character scalars
  (`c("2017-09-01", "2017-09-15")`) constraining `start_time`.

- min_wind_kt:

  Numeric scalar: minimum `peak_wind_kt` at the site.

## Value

Tibble of matching rows from `out$events` with columns `storm_id`,
`location`, `start_time`, `peak_wind_kt`, and `storm_intensity_kt`,
sorted by `start_time`. Use `storm_id` values directly in subsequent
query functions.

## See also

[`query_storm_track_years`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md),
[`query_impact_years`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Find Irma by approximate date and intensity at Saba
lookup_storm_id(
  out          = hazard_out,
  year         = 2017,
  location     = "Saba",
  date_range   = c("2017-09-01", "2017-09-15"),
  min_wind_kt  = 50
)
# Returns storm_id "2017242N16333" — pass this to query_storm_track_years()

# All major hurricanes (>= 96 kt) at any location
lookup_storm_id(out = hazard_out, min_wind_kt = 96)
} # }
```
