# Enforce monotone wind radii.

Adjusts directional radii so that `R64_km <= R50_km <= R34_km` wherever
paired values are finite.

## Usage

``` r
.enforce_monotone_radii(R34_km, R50_km, R64_km)
```

## Arguments

- R34_km:

  Numeric vector of 34-kt radii in kilometres.

- R50_km:

  Numeric vector of 50-kt radii in kilometres.

- R64_km:

  Numeric vector of 64-kt radii in kilometres.

## Value

List with elements `R34_km`, `R50_km`, and `R64_km` after monotonicity
enforcement.
