# Aggregate track points into storm events.

Summarises track-point data to one row per storm with peak site wind,
storm intensity, pressure, timing, and radius-of-maximum-wind
diagnostics.

## Usage

``` r
make_storm_events(track_df)
```

## Arguments

- track_df:

  Track-point tibble or data frame with at least `SID` and `iso_time`,
  plus optional wind and pressure fields.

## Value

Tibble with one row per storm and event-level summary attributes.

## See also

[`classify_severity`](https://tanerumit.github.io/ipdcstorm/reference/classify_severity.md),
[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md)

## Examples

``` r
track_df <- data.frame(
  SID = c("AL012000", "AL012000"),
  iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
  V_site_kt = c(40, 55),
  wind_kt = c(60, 65)
)
make_storm_events(track_df)
#> # A tibble: 1 × 11
#>   storm_id start_time          end_time            n_points peak_wind_kt
#>   <chr>    <dttm>              <dttm>                 <int>        <dbl>
#> 1 AL012000 2000-08-01 00:00:00 2000-08-01 06:00:00        2           55
#> # ℹ 6 more variables: storm_intensity_kt <dbl>, min_pressure_hpa <dbl>,
#> #   pressure_deficit_hpa <dbl>, rmw_min_km <dbl>, rmw_mean_km <dbl>, year <dbl>
```
