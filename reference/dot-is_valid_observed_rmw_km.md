# Validate observed radius of maximum wind.

Flags observed radius-of-maximum-wind values that fall inside the
accepted physical range.

## Usage

``` r
.is_valid_observed_rmw_km(rmw_km)
```

## Arguments

- rmw_km:

  Numeric vector of observed radius of maximum wind values in
  kilometres.

## Value

Logical vector indicating whether each value is finite and within the
accepted bounds.
