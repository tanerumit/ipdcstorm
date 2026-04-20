# Compute Poisson rate table by storm class.

Summarises annual event counts into storm-class-specific Poisson rates
and annual exceedance probabilities.

## Usage

``` r
compute_lambda_table(annual_counts)
```

## Arguments

- annual_counts:

  Tibble returned by
  [`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md).

## Value

Tibble with `storm_class`, `lambda`, `n_years`, `prob_annual`, and
`prob_none` columns.

## See also

[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md),
[`estimate_k_hat`](https://tanerumit.github.io/ipdcstorm/reference/estimate_k_hat.md)

## Examples

``` r
annual_counts <- data.frame(
  year = c(2000, 2001, 2000, 2001),
  storm_class = c("TS", "TS", "HUR", "HUR"),
  n_events = c(1, 0, 0, 1)
)
compute_lambda_table(annual_counts)
#> # A tibble: 2 × 5
#>   storm_class lambda n_years prob_annual prob_none
#>   <chr>        <dbl>   <int>       <dbl>     <dbl>
#> 1 HUR            0.5       2       0.393     0.607
#> 2 TS             0.5       2       0.393     0.607
```
