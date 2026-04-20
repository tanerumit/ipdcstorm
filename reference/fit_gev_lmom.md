# Fit GEV distribution using L-moments (Hosking 1990)

Estimates GEV parameters (location \\\mu\\, scale \\\sigma\\, shape
\\\xi\\) using the method of L-moments. This is more robust than MLE for
small samples (n \< 50) and requires no optimization. The sign
convention follows the standard meteorological form: \\\xi \> 0\\ =
heavy tail (Frechet), \\\xi \< 0\\ = bounded tail (Weibull), \\\xi = 0\\
= Gumbel.

## Usage

``` r
fit_gev_lmom(x, xi_bounds = c(-0.4, 0.5))
```

## Arguments

- x:

  Numeric vector of block maxima (e.g., annual maxima).

- xi_bounds:

  Numeric vector of length 2; allowed range for shape parameter. Default
  c(-0.5, 0.5) prevents extreme tail behavior.

## Value

A list with elements:

- mu:

  Location parameter.

- sigma:

  Scale parameter.

- xi:

  Shape parameter.

- n:

  Sample size.

- l_moments:

  Named vector of L1, L2, L3, tau3.

- converged:

  Logical; whether estimation produced valid parameters.
