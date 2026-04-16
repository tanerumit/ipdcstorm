# Reference HURDAT2/literature annual rates for Leeward Islands region

Returns a tibble of published annual TS passage rates from literature,
for comparison against model-fitted lambdas.

## Usage

``` r
get_reference_rates()
```

## Value

Tibble with columns: region, storm_class, lambda_ref, source,
gate_approx_nm, period.
