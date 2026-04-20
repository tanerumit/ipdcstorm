# Compute return levels using a hurdle-GEV model

Accounts for years with zero impact (no events) separately from the
intensity distribution. The model is: P(annual_max \<= v) = p0 + (1 -
p0) \* F_GEV(v) where p0 is the probability of a zero-event year, and
F_GEV is the GEV CDF fitted to nonzero annual maxima.

## Usage

``` r
compute_return_levels_gev(
  annual_max,
  return_periods = c(5, 10, 25, 50),
  xi_bounds = c(-0.3, 0.4)
)
```

## Arguments

- annual_max:

  Numeric vector of annual maxima (including zeros).

- return_periods:

  Numeric vector of return periods (years).

- xi_bounds:

  Bounds on GEV shape parameter (default: c(-0.3, 0.4) for TSs).

## Value

A list with:

- return_levels:

  Named numeric vector of return levels.

- gev_fit:

  GEV fit object for nonzero maxima.

- p_zero:

  Fraction of zero-event years.

- n_total, n_nonzero:

  Sample sizes.
