# Compute storm heading by track point.

Adds a `heading_deg` column based on great-circle motion between
successive track points within each storm.

## Usage

``` r
compute_storm_heading(df)
```

## Arguments

- df:

  Data frame with columns `SID`, `iso_time`, `lat`, and `lon`.

## Value

Data frame with the original columns plus numeric `heading_deg` in
degrees clockwise from north.

## See also

[`compute_site_winds_full`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md)

## Examples

``` r
df <- data.frame(
  SID = c("AL012000", "AL012000"),
  iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
  lat = c(18, 18.2),
  lon = c(-63, -63.2)
)
compute_storm_heading(df)
#> # A tibble: 2 × 5
#>   SID      iso_time              lat   lon heading_deg
#>   <chr>    <dttm>              <dbl> <dbl>       <dbl>
#> 1 AL012000 2000-08-01 00:00:00  18   -63          316.
#> 2 AL012000 2000-08-01 06:00:00  18.2 -63.2        316.
```
