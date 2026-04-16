# Resolve Holland outer cutoff radius.

Computes the outer cutoff radius in kilometres used to taper winds
beyond the resolved 34-kt radius.

## Usage

``` r
.resolve_holland_outer_cutoff_km(R34_km, R34_is_fallback = FALSE)
```

## Arguments

- R34_km:

  Numeric vector of 34-kt wind radii in kilometres.

- R34_is_fallback:

  Logical vector indicating whether each `R34_km` value comes from
  climatological fallback rather than observations.

## Value

Numeric vector of outer cutoff radii in kilometres.
