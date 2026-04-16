# Cap inferred radius of maximum wind.

Applies intensity-dependent lower and upper bounds to inferred
radius-of-maximum-wind values.

## Usage

``` r
.cap_inferred_rmw_km(rmw_km, Vmax_kt)
```

## Arguments

- rmw_km:

  Numeric vector of inferred radius of maximum wind values in
  kilometres.

- Vmax_kt:

  Numeric vector of maximum sustained wind speeds in knots.

## Value

Numeric vector of bounded radius of maximum wind values in kilometres.
