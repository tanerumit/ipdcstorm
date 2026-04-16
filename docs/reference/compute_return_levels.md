# Compute empirical return levels from a vector of annual maxima

Extracts return levels from a sample of annual maxima using the
plotting-position approach: sort values, assign return periods T =
(n+1)/rank, interpolate to requested return periods.

## Usage

``` r
compute_return_levels(annual_max, return_periods = c(5, 10, 25, 50))
```

## Arguments

- annual_max:

  Numeric vector of annual maximum values (e.g., peak_wind_kt).

- return_periods:

  Numeric vector of return periods (years) to compute.

## Value

Named numeric vector of return levels at requested periods.
