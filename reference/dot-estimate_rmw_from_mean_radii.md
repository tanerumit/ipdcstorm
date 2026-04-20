# Estimate radius of maximum wind from mean wind radii.

Infers the radius of maximum wind from storm-wide mean 64-, 50-, or
34-kt radii using fixed regression coefficients.

## Usage

``` r
.estimate_rmw_from_mean_radii(R64_mean_km, R50_mean_km, R34_mean_km)
```

## Arguments

- R64_mean_km:

  Numeric vector of mean 64-kt radii in kilometres.

- R50_mean_km:

  Numeric vector of mean 50-kt radii in kilometres.

- R34_mean_km:

  Numeric vector of mean 34-kt radii in kilometres.

## Value

Numeric vector of inferred radius of maximum wind values in kilometres.
