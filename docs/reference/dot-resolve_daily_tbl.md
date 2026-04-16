# Coerce daily hazard output to a flat tibble

Accepts either a named list of tibbles from
[`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)
or a single tibble and returns one flat tibble, optionally filtered to
the requested locations.

## Usage

``` r
.resolve_daily_tbl(daily, location = NULL)
```

## Arguments

- daily:

  Named list of tibbles or single tibble.

- location:

  Character vector of locations to retain, or `NULL`.

## Value

Tibble with the resolved daily rows.
