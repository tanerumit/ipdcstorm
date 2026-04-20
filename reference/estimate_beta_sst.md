# Estimate the SST-activity scaling coefficient I2_SST

Fits a Poisson (or negative binomial) GLM of annual TS counts on MDR SST
anomaly to estimate I2_SST in:

`E[N_t] = exp(alpha + beta_SST * sst_anomaly_t)`

The coefficient I2_SST represents the log-linear sensitivity of annual
TS activity to SST anomalies. Typical values from the literature are
around 0.4-0.8 per degC for the North Atlantic (Villarini et al. 2011;
Vecchi et al. 2021).

If the `MASS` package is available, a negative binomial GLM is preferred
(accounts for overdispersion). Otherwise, a quasi-Poisson GLM is used.

## Usage

``` r
estimate_beta_sst(
  annual_counts,
  sst_df,
  min_year = 1970L,
  beta_prior = NULL,
  verbose = TRUE
)
```

## Arguments

- annual_counts:

  Tibble from
  [`compute_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md)
  with columns: year, storm_class, n_events.

- sst_df:

  Tibble with columns: year, sst_anomaly.

- min_year:

  Integer; earliest year to include (default: cfg\$start_year).

- beta_prior:

  Optional numeric; if provided, shrinks the estimate toward this prior
  value using a simple Bayesian-inspired weighted average: beta_final =
  w \* beta_mle + (1 - w) \* beta_prior, where w = min(1, n_years/50).
  This stabilizes estimates when the historical record is short.
  Climate-rate calibration must use basin-consistent annual counts in
  storms/year; passing target-conditioned counts with a `location`
  column is rejected because that duplicates annual activity across
  targets.

- verbose:

  Logical; print diagnostic output.

## Value

A list with:

- beta_sst:

  Estimated (or shrunk) I2_SST coefficient.

- beta_se:

  Standard error of the MLE I2_SST.

- beta_mle:

  Raw MLE estimate before shrinkage.

- alpha:

  Intercept (log baseline rate).

- method:

  Character; "negbin", "quasipoisson", or "literature_fallback".

- n_years:

  Number of years used in estimation.

- r_squared_dev:

  Deviance-based pseudo-RA2 (proportion of deviance explained).

- aic:

  AIC of the fitted model (NA for quasi-Poisson).

- annual_count_series:

  Tibble of the basin-consistent annual total count series used for
  calibration, in storms/year.

- annual_count_source:

  Character label describing the provenance of `annual_count_series`.

- guardrail:

  List describing any regularization or fallback applied to keep
  `beta_sst` within documented plausibility bounds.

- fit_data:

  Tibble of the joined data used for fitting.

- glm_fit:

  The fitted GLM object.

## Examples

``` r
# Small example workflow:
sst_df <- get_mdr_sst_builtin() |>
  compute_sst_anomaly()
# annual_counts should come from compute_annual_counts()
# beta_info <- estimate_beta_sst(annual_counts, sst_df)
```
