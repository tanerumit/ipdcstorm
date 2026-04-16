# Compare model-fitted rates against reference climatologies

Takes the lambda table from the hazard model output and compares
per-location per-storm_class rates against published references.

## Usage

``` r
validate_rates(out, ref_rates = NULL)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- ref_rates:

  Optional tibble of reference rates (default:
  [`get_reference_rates()`](https://tanerumit.github.io/ipdcstorm/reference/get_reference_rates.md)).

## Value

Tibble with model vs reference rate comparison.
