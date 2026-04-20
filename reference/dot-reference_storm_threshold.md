# Derive an impact threshold for a reference storm

Resolves a numeric threshold for a named metric from either:

1.  The historical event record in `out$events` (for
    `metric = "peak_wind_kt"`).

2.  The rows of `daily` where that storm appears as the sampled event
    (for damage metrics), taking the median annual metric across those
    occurrences.

## Usage

``` r
.reference_storm_threshold(daily, out, storm_id, metric, location = NULL)
```

## Arguments

- daily:

  Flat tibble (already resolved by
  [`.resolve_daily_tbl()`](https://tanerumit.github.io/ipdcstorm/reference/dot-resolve_daily_tbl.md)).

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- storm_id:

  Character scalar IBTrACS SID of the reference storm.

- metric:

  Character scalar metric name; one of `"peak_wind_kt"`, `"cum_damage"`,
  or `"max_damage_rate"`.

- location:

  Character vector of locations to consider, or `NULL` for all.

## Value

Named numeric vector of thresholds, one entry per location.
