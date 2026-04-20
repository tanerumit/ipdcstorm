# Estimate I3 (intensity shift coefficient) from historical data

Fits a logistic regression of annual hurricane fraction on SST anomaly
to estimate the intensification trend coefficient I3 in:

`p_HUR(t) = p_HUR_base * (1 + gamma * sst_anomaly_t)`

The coefficient I3 captures how the fraction of storms reaching
hurricane intensity shifts with SST warming. Literature estimates
suggest roughly 5-8% increase in Cat 4-5 fraction per degC of tropical
SST warming (Knutson et al. 2020; Kossin et al. 2020).

Uses a binomial GLM: cbind(n_HUR, n_TS) ~ sst_anomaly, then converts the
logistic coefficient to the linear I3 parameterization.

## Usage

``` r
estimate_gamma_intensity(
  annual_counts,
  sst_df,
  min_year = 1970L,
  gamma_prior = 0.065,
  verbose = TRUE
)
```

## Arguments

- annual_counts:

  Tibble with columns: year, storm_class, n_events.

- sst_df:

  Tibble with columns: year, sst_anomaly.

- min_year:

  Integer; earliest year to include.

- gamma_prior:

  Optional numeric; prior value for shrinkage (default: 0.065, i.e. 6.5%
  increase in HUR fraction per degC).

- verbose:

  Logical.

## Value

A list with:

- gamma:

  Estimated (or shrunk) I3 coefficient.

- gamma_mle:

  Raw MLE from logistic regression (converted to linear scale).

- gamma_se:

  Approximate standard error on I3.

- p_hur_base:

  Baseline hurricane fraction.

- method:

  Character; estimation method used.

- n_years:

  Number of years in fit.

- fit_data:

  Tibble used for fitting.
