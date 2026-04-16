# Derive annual counts from model output.

Builds location-specific annual storm counts from `out$events` and
zero-fills missing year-class combinations within each location.

## Usage

``` r
get_annual_counts(out)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
  containing an `events` component.

## Value

Tibble with columns `location`, `year`, `storm_class`, and `n_events`.

## See also

[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md),
[`compute_lambda_table`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)

## Examples

``` r
out <- list(events = data.frame(
  location = c("Saba", "Saba"),
  year = c(2000, 2001),
  storm_class = c("TS", "HUR"),
  storm_id = c("A", "B")
))
get_annual_counts(out)
#> # A tibble: 4 × 4
#>   location  year storm_class n_events
#>   <chr>    <dbl> <chr>          <int>
#> 1 Saba      2000 HUR                0
#> 2 Saba      2000 TS                 1
#> 3 Saba      2001 HUR                1
#> 4 Saba      2001 TS                 0
```
