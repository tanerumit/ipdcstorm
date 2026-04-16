# Add day/month/year fields to daily hazard data

Add day/month/year fields to daily hazard data

## Usage

``` r
prep_daily(daily_impact)
```

## Arguments

- daily_impact:

  A data frame or tibble with at least a `date` column coercible by
  [`format()`](https://rdrr.io/r/base/format.html) (typically `Date`),
  one row per day.

## Value

`daily_impact` with three added integer columns: `doy` (1-366), `month`
(1-12), and `year` (4-digit calendar year).
