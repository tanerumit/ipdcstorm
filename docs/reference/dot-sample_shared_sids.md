# Sample SIDs from the shared event pool

Sample SIDs from the shared event pool

## Usage

``` r
.sample_shared_sids(pool, storm_class, n)
```

## Arguments

- pool:

  Tibble from
  [`.build_shared_event_pool()`](https://tanerumit.github.io/ipdcstorm/reference/dot-build_shared_event_pool.md).

- storm_class:

  Character scalar; one of `"TS"` or `"HUR"`.

- n:

  Integer scalar number of SIDs to draw (with replacement).

## Value

Character vector of length `n`, or `character(0)` when `n == 0` or the
class sub-pool is empty.
