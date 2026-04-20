# Compute annual storm-event counts.

Counts unique storm events per year and storm class, then completes the
series so missing year-class combinations are filled with zeros.

## Usage

``` r
compute_annual_counts(events, storm_classes = c("TS", "HUR"))
```

## Arguments

- events:

  Tibble with at least `year`, `storm_class`, and `storm_id` columns.

- storm_classes:

  Character vector of storm classes to include.

## Value

Tibble with columns `year`, `storm_class`, and `n_events`, completed
over all represented years and classes.

## See also

[`compute_lambda_table`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md),
[`get_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/get_annual_counts.md)

## Examples

``` r
events <- data.frame(
  year = c(2000, 2000, 2001),
  storm_class = c("TS", "TS", "HUR"),
  storm_id = c("A", "B", "C")
)
compute_annual_counts(events)
#> # A tibble: 4 × 3
#>    year storm_class n_events
#>   <dbl> <chr>          <int>
#> 1  2000 HUR                0
#> 2  2000 TS                 2
#> 3  2001 HUR                1
#> 4  2001 TS                 0
```
