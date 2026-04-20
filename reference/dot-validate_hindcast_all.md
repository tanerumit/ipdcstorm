# Run hindcast validation across all locations in a hazard model output

Run hindcast validation across all locations in a hazard model output

## Usage

``` r
.validate_hindcast_all(
  out,
  holdout_years = 10,
  n_sim = 5000,
  return_periods = c(5, 10, 25, 50),
  conf_level = 0.95,
  seed = 42,
  beta_sst = 0,
  gamma_intensity = 0,
  use_raw_rates = TRUE,
  xi_bounds = c(-0.3, 0.4),
  n_boot = 500
)
```
