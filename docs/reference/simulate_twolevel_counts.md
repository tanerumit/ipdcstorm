# Simulate annual counts with climate adjustments

Non-stationary extension of the Poisson-Gamma annual count model with a
rate effect and an intensity effect applied to a stationary time slice:

**Rate effect:** Each year's activity factor is modulated by a scalar
SST shift: `A_t = activity_factor * exp(beta_SST * delta_sst)`

**Intensity effect:** The storm-class split varies with the same scalar
SST shift:
`p_HUR(t) = clamp(p_HUR_base * (1 + gamma * delta_sst), 0.01, 0.99)`
`N_total_t ~ Poisson(lambda_total * A_t)`
`n_HUR_t ~ Binomial(N_total_t, p_HUR(t))` `n_TS_t = N_total_t - n_HUR_t`

When gamma_intensity is 0, the storm-class split is constant. When both
beta_sst and gamma_intensity are 0, reduces to stationary model.

## Usage

``` r
simulate_twolevel_counts(
  lambda_table,
  k_hat,
  n_years_sim,
  delta_sst = 0,
  beta_sst = 0,
  gamma_intensity = 0,
  p_hur_base = NULL,
  .sst_abs_max = 10,
  .sst_scale_max = 1000,
  .mu_total_max = 1e+06
)
```

## Arguments

- lambda_table:

  Tibble from
  [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md).

- k_hat:

  Numeric; Gamma shape parameter for overdispersion.

- n_years_sim:

  Integer; number of years to simulate.

- delta_sst:

  Numeric scalar SST shift (degC) applied uniformly to all simulation
  years.

- beta_sst:

  Numeric; SST rate-effect coefficient.

- gamma_intensity:

  Numeric; SST intensity-effect coefficient. Represents fractional
  change in p_HUR per degC of SST warming.

- p_hur_base:

  Optional numeric; baseline hurricane fraction. If NULL, computed from
  lambda_table.

- .sst_abs_max:

  Numeric guardrail for absolute SST anomaly magnitude.

- .sst_scale_max:

  Numeric guardrail for `exp(beta_sst * delta_sst)`.

- .mu_total_max:

  Numeric guardrail for annual Poisson mean.

## Value

Tibble with columns: sim_year, activity_factor, climate_scale,
activity_combined, p_hurricane, n_total, n_ts, n_hur. The applied
`delta_sst` is recorded as an attribute on the returned tibble.
