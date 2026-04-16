# Plot empirical wind return levels

Computes empirical return-period points using either annual block maxima
or peaks-over-threshold exceedances and plots them on a log-scaled
return-period axis.

## Usage

``` r
plot_return_levels(daily, block_maxima = TRUE, threshold = NULL)
```

## Arguments

- daily:

  A data frame or tibble with columns `sim_year` and `wind_kt` (numeric,
  knots).

- block_maxima:

  Logical scalar; if `TRUE`, use annual maxima, otherwise use
  exceedances over `threshold`.

- threshold:

  Optional numeric scalar (knots); required when `block_maxima = FALSE`.

## Value

A `ggplot` object.

## Details

For `block_maxima = TRUE`, one annual maximum is taken per `sim_year`.
For `block_maxima = FALSE`, `threshold` must be non-`NULL`, and only
rows with `wind_kt > threshold` are used. Wind summaries use
`na.rm = TRUE`.

## Examples

``` r
d <- data.frame(
  sim_year = rep(2001:2005, each = 3),
  wind_kt = c(20, 35, 50, 15, 30, 70, 10, 25, 45, 18, 40, 55, 22, 60, 65)
)
plot_return_levels(d, block_maxima = TRUE)
```
