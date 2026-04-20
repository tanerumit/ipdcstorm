# Flag tropical cyclone or hurricane days

Flags days in a daily hazard table that are associated with any tropical
storm or hurricane event.

## Usage

``` r
is_tc_day(daily)

is_hur_day(daily)
```

## Arguments

- daily:

  Tibble/data.frame returned by
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md).

## Value

Logical vector identifying tropical-storm or hurricane days.

## See also

[`generate_daily_hazard_impact_spatial`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)

## Examples

``` r
daily <- tibble::tibble(event_class = c(NA, "TS", "HUR"))
is_tc_day(daily)
#> Error in is_tc_day(daily): could not find function "is_tc_day"
is_hur_day(daily)
#> Error in is_hur_day(daily): could not find function "is_hur_day"
```
