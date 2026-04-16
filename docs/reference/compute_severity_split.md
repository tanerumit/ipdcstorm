# Compute time-varying storm-class split

Computes per-year hurricane probability and corresponding TS/HUR rates:

`p_HUR(t) = clamp(p_HUR_base * (1 + gamma * sst_anomaly_t), 0.01, 0.99)`
`lambda_HUR(t) = lambda_total * p_HUR(t)`
`lambda_TS(t) = lambda_total * (1 - p_HUR(t))`

## Usage

``` r
compute_severity_split(lambda_table, sst_anomaly, gamma = 0, p_hur_base = NULL)
```

## Arguments

- lambda_table:

  Tibble from
  [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md).

- sst_anomaly:

  Numeric vector of SST anomalies per simulation year (degC).

- gamma:

  Numeric; intensity shift coefficient.

- p_hur_base:

  Optional numeric; if NULL, computed from lambda_table.

## Value

Tibble with columns: sim_year, sst_anomaly, p_hur, lam_TS, lam_HUR.
