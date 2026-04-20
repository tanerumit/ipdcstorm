# Compute daily exposure hours above a wind threshold

Converts a daily wind threshold exceedance test into daily exposure
hours at package daily resolution.

## Usage

``` r
exposure_hours(daily, threshold_kt, use_gust = FALSE)
```

## Arguments

- daily:

  Tibble/data.frame returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md).

- threshold_kt:

  Numeric scalar wind threshold (kt).

- use_gust:

  Logical scalar; if `TRUE`, use `wind_gust_kt`.

## Value

Numeric vector of daily exposure hours at package daily resolution.

## Examples

``` r
daily <- tibble::tibble(wind_kt = c(20, 40), wind_gust_kt = c(24, 48))
exposure_hours(daily, threshold_kt = 34)
#> Error in exposure_hours(daily, threshold_kt = 34): could not find function "exposure_hours"
```
