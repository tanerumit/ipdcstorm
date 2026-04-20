# Build a shared event pool from multiple location event libraries

Unions all events from a named list of location event libraries, keeping
one row per SID (the row with the highest site-level peak wind across
all locations). Only `"TS"` and `"HUR"` class events are retained.

## Usage

``` r
.build_shared_event_pool(libs)
```

## Arguments

- libs:

  Named list of event libraries as returned by
  [`build_event_library_from_out()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library_from_out.md),
  one element per location.

## Value

Tibble with one row per unique SID and a `storm_class` column alongside
any shared event attributes carried from the library tables.
