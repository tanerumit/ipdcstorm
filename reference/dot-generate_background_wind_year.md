# Generate one year of correlated background wind for multiple locations

Simulates a \\K\\-variate AR(1) in standardised normal space with
spatial covariance given by the Gaussian copula, then transforms each
margin through its site- and month-specific Weibull quantile function.

## Usage

``` r
.generate_background_wind_year(year, location, cfg)
```

## Arguments

- year:

  Integer calendar year (determines leap year and month mapping).

- location:

  Character vector of location names; must be a subset of
  `cfg$locations`.

- cfg:

  A `background_wind_cfg` object.

## Value

Named list of numeric vectors (one per location), each of length 365 or
366, containing background wind speeds in kt.
