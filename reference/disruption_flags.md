# Compute disruption flags from daily hazard output

Adds logical disruption indicators to a daily hazard table by comparing
daily wind or surge values with user-supplied thresholds.

## Usage

``` r
disruption_flags(
  daily,
  thr_port = NA_real_,
  thr_infra = NA_real_,
  thr_surge = NA_real_,
  use_gust = FALSE
)
```

## Arguments

- daily:

  Tibble/data.frame returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md).

- thr_port:

  Numeric scalar wind threshold (kt) for port disruption.

- thr_infra:

  Numeric scalar wind threshold (kt) for infrastructure disruption.

- thr_surge:

  Numeric scalar surge threshold (m) for surge disruption.

- use_gust:

  Logical scalar; if `TRUE`, use `wind_gust_kt` instead of `wind_kt`.

## Value

Input tibble with added logical disruption columns.

## See also

[`generate_daily_hazard_impact_spatial`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)

## Examples

``` r
daily <- tibble::tibble(
  wind_kt = c(20, 45),
  wind_gust_kt = c(24, 54),
  surge_m = c(0.1, 0.8)
)
disruption_flags(daily, thr_port = 34, thr_infra = 50, thr_surge = 0.5)
#> Error in disruption_flags(daily, thr_port = 34, thr_infra = 50, thr_surge = 0.5): could not find function "disruption_flags"
```
