# Extract TS/HUR lambdas from a lambda table

Extract TS/HUR lambdas from a lambda table

## Usage

``` r
.extract_lambdas(lambda_table)
```

## Arguments

- lambda_table:

  Tibble from
  [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)
  with severities "TS" and "HUR".

## Value

List with elements `ts`, `hur`, `total`.
