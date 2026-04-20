# Add day/month/year fields to daily hazard data

Add day/month/year fields to daily hazard data

## Usage

``` r
prep_daily(daily_impact)
```

## Arguments

- daily_impact:

  A data frame or tibble with either a `date` column or both `sim_year`
  and `doy` columns (from
  [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)).

## Value

`daily_impact` with three added integer columns: `doy` (1-366), `month`
(1-12), and `year` (4-digit calendar year).
