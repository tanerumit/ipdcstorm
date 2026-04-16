# Extract one annual metric from shared daily-year summaries

Extract one annual metric from shared daily-year summaries

## Usage

``` r
.annual_metric_tbl(daily_tbl, metric)
```

## Arguments

- daily_tbl:

  Flat tibble from
  [`.resolve_daily_tbl()`](https://tanerumit.github.io/ipdcstorm/reference/dot-resolve_daily_tbl.md).

- metric:

  Character scalar naming the annual metric to extract.

## Value

Tibble with `location`, `sim_year`, and `annual_metric`.
