# Build multiplicative lambda scalers from a rate-check table

Build multiplicative lambda scalers from a rate-check table

## Usage

``` r
.lambda_scalers_from_rate_check(
  rate_tbl,
  scale_min = 0.25,
  scale_max = 4,
  scaling_mode = "target"
)
```

## Arguments

- rate_tbl:

  Tibble with location, storm_class, lambda_model or lambda_model_raw,
  lambda_ref, and optional expected_ratio.

- scale_min, scale_max:

  Numeric scalar clamp bounds for lambda_scale.

- scaling_mode:

  Character scalar. `"target"` preserves historical behavior;
  `"down_only"` prevents lambda upscaling above modeled rates.

## Value

Tibble keyed by location and storm_class with raw lambda, target lambda,
lambda_scale, adjusted lambda, and scaler status.
