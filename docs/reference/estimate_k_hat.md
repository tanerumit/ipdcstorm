# Estimate annual-count overdispersion.

Estimates the negative-binomial overdispersion parameter from total
annual storm counts using the identity `Var(N) = mu + mu^2 / k`.

## Usage

``` r
estimate_k_hat(annual_counts)
```

## Arguments

- annual_counts:

  Tibble returned by
  [`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md).

## Value

List with elements `k_hat`, `annual_total`, `mu`, and `var`.

## See also

[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md),
[`compute_lambda_table`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)

## Examples

``` r
annual_counts <- data.frame(
  year = c(2000, 2000, 2001, 2001),
  storm_class = c("TS", "HUR", "TS", "HUR"),
  n_events = c(1, 0, 2, 1)
)
estimate_k_hat(annual_counts)
#> $k_hat
#> [1] 1e+06
#> 
#> $annual_total
#> # A tibble: 2 × 2
#>    year     N
#>   <dbl> <dbl>
#> 1  2000     1
#> 2  2001     3
#> 
#> $mu
#> [1] 2
#> 
#> $var
#> [1] 2
#> 
```
