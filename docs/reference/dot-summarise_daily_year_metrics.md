# Summarise daily hazard output to annual metrics

Produces one row per `(location, sim_year)` with the annual metrics used
by the query and stress-selection workflows.

## Usage

``` r
.summarise_daily_year_metrics(daily_tbl)
```

## Arguments

- daily_tbl:

  Flat tibble from
  [`.resolve_daily_tbl()`](https://tanerumit.github.io/ipdcstorm/reference/dot-resolve_daily_tbl.md).

## Value

Tibble with annual hazard and impact summaries.
