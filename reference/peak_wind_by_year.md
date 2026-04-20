# Summarise peak wind per simulation year

Aggregates a daily hazard table to annual peak sustained wind by
location, simulation year, and scenario.

## Usage

``` r
peak_wind_by_year(daily)
```

## Arguments

- daily:

  Tibble/data.frame returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md).

## Value

Tibble with one row per location and simulation year and a
`peak_wind_kt` column.

## Examples

``` r
daily <- tibble::tibble(
  location = c("Saba", "Saba"),
  sim_year = c(1, 1),
  scenario = c("baseline", "baseline"),
  wind_kt = c(20, 55)
)
peak_wind_by_year(daily)
#> Error in peak_wind_by_year(daily): could not find function "peak_wind_by_year"
```
