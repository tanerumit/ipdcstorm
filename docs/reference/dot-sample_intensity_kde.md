# Sample from a fitted intensity KDE

Sample from a fitted intensity KDE

## Usage

``` r
.sample_intensity_kde(fit, n)
```

## Arguments

- fit:

  List from `.fit_intensity_kde()`.

- n:

  Integer; number of draws.

## Value

Numeric vector of n intensity draws within the interval `fit$lower` to
`fit$upper`.
