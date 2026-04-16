# Estimate the historical hurricane fraction from annual counts

Computes p_HUR_base = I\>\_HUR / (I\>\_TS + I\>\_HUR) from the
historical record. This is the baseline probability that an event
reaching TS intensity or above becomes a hurricane (\>=64 kt).

## Usage

``` r
compute_p_hur_base(lambda_table)
```

## Arguments

- lambda_table:

  Tibble from
  [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)
  with severities "TS" and "HUR".

## Value

Numeric scalar: baseline hurricane fraction p_HUR_base.
