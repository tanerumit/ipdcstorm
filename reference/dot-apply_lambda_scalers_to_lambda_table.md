# Apply site/class lambda scalers to a TS/HUR lambda table

Apply site/class lambda scalers to a TS/HUR lambda table

## Usage

``` r
.apply_lambda_scalers_to_lambda_table(
  lambda_table,
  location,
  lambda_scalers = NULL
)
```

## Arguments

- lambda_table:

  Tibble from compute_lambda_table().

- location:

  Character scalar; site name.

- lambda_scalers:

  Optional output from .lambda_scalers_from_rate_check().

## Value

Tibble with adjusted lambda values. Existing columns are preserved;
lambda_raw, lambda_scale, and lambda_adj are added for traceability.
