# Background wind configuration for correlated Weibull generation

Specifies monthly Weibull marginal parameters per location and an
optional Gaussian copula correlation matrix for generating spatially
correlated background wind on all days. Background winds are combined
with storm pulses via `pmax`, so the storm signal always dominates on
active days.

The Gaussian copula workflow per simulated year:

1.  Simulate a \\K\\-variate AR(1) process in standardised normal space,
    using the Cholesky factor of `cor_matrix` to impose spatial
    correlation and `ar1` for day-to-day persistence.

2.  Map each margin through
    [`pnorm()`](https://rdrr.io/r/stats/Normal.html) to obtain uniform
    scores.

3.  Transform each score through its site- and month-specific Weibull
    quantile function to produce background wind speeds in kt.

## Usage

``` r
make_background_wind_cfg(weibull_params, cor_matrix = NULL, ar1 = 0)
```

## Arguments

- weibull_params:

  Named list of data frames, one per location. Each data frame must have
  columns `month` (integer 1-12), `shape` (Weibull shape \> 0), and
  `scale` (Weibull scale \> 0).

- cor_matrix:

  Numeric Pearson correlation matrix for the Gaussian copula. Row and
  column names must match `names(weibull_params)`. Must be symmetric
  positive-definite with unit diagonal. `NULL` (default) treats
  locations as independent.

- ar1:

  Numeric scalar in `[0, 1)`; AR(1) coefficient for day-to-day
  persistence in the normal domain. `0` (default) gives independent
  daily draws.

## Value

A list of class `"background_wind_cfg"`.
