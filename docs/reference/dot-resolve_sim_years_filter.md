# Resolve sim_years filter to a per-location list

Accepts either (a) an integer vector applied uniformly across all
locations, or (b) a tibble with at least `location` and `sim_year`
columns (the output of
[`query_storm_track_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_storm_track_years.md)
or
[`query_impact_years()`](https://tanerumit.github.io/ipdcstorm/reference/query_impact_years.md))
and returns a named list mapping each location to its allowed `sim_year`
values.

## Usage

``` r
.resolve_sim_years_filter(sim_years, locations)
```

## Arguments

- sim_years:

  Integer vector, tibble with `location`/`sim_year` columns, or `NULL`
  (no filter).

- locations:

  Character vector of all locations present in `daily`.

## Value

Named list: `list(loc = integer_vector, ...)` or `NULL`.
