# Validate wind field estimates against station observations

For each reference observation, finds the matching storm in the model's
trackpoint data, extracts the model's peak site wind for that storm at
the target location, and compares against the observed wind. Each
observation carries a quality tier (A/B/C) that determines the
acceptable bias threshold.

## Usage

``` r
validate_wind_field(out, obs_table = NULL)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- obs_table:

  Optional tibble of observations (default:
  [`get_wind_observations()`](https://tanerumit.github.io/ipdcstorm/reference/get_wind_observations.md)).

## Value

Tibble with model vs observed comparison per storm-station pair,
including `obs_quality`, `bias_threshold_pct`, and `bias_ok` columns.
