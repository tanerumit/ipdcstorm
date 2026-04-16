# Precompute spatial event lookup tables for fast SID resolution

Precompute spatial event lookup tables for fast SID resolution

## Usage

``` r
.prepare_spatial_event_lookup(lib)
```

## Arguments

- lib:

  Event library from
  [`build_event_library_from_out()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library_from_out.md).

## Value

Input library with a cached `spatial_lookup` entry.
