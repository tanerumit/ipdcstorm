# Perturb sampled event properties for climate sensitivity

Adjusts individual storm properties in a sampled-event tibble to reflect
projected changes in storm structure under warmer SSTs.

Perturbations applied (all multiplicative, proportional to delta_SST):
V_peak: scaled by (1 + v_scale \* delta_SST) RMW_mean_km: scaled by (1 +
r_scale \* delta_SST) dur_days: scaled by 1 / (1 + speed_scale \*
delta_SST / 2) precip_scaling: new column = 1 + precip_scale \*
delta_SST

At delta_SST = 0 all factors equal 1 (identity property for validation).

## Usage

``` r
perturb_event(events, delta_sst, cc_params = NULL)
```

## Arguments

- events:

  Tibble of sampled events with at least V_peak, dur_days.

- delta_sst:

  Numeric scalar; SST anomaly (degC) for the simulation year.

- cc_params:

  Named list of per-degC scaling factors (NULL = defaults).

## Value

The input tibble with perturbed columns plus precip_scaling and
delta_sst columns.
