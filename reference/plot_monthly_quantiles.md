# Plot monthly wind quantile summary

Computes monthly wind quantiles and plots a median-to-upper ribbon with
upper and extreme quantile lines.

## Usage

``` r
plot_monthly_quantiles(
  daily,
  probs = c(0.5, 0.95, 0.99),
  log_scale = FALSE,
  thr_tc = 34,
  thr_hur = 64
)
```

## Arguments

- daily:

  A data frame or tibble with columns `date` (`Date`) and `wind_kt`
  (numeric, knots; `NA` allowed and ignored in quantiles).

- probs:

  Numeric vector of length 3 with probabilities in `[0, 1]`, used in
  order `probs[1]`, `probs[2]`, `probs[3]`.

- log_scale:

  Logical scalar; if `TRUE`, use a log10 y-axis.

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

## Value

A `ggplot` object.

## Details

`probs` is indexed as lower/median display (`probs[1]`), ribbon upper
bound (`probs[2]`), and extreme line (`probs[3]`). If `log_scale` is
`TRUE`, quantile values are floored at 0.5 before
[`scale_y_log10()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
to avoid non-positive values.

## Examples

``` r
d <- data.frame(
  date = as.Date("2001-01-01") + 0:59,
  wind_kt = pmax(1, 20 + rnorm(60, sd = 5))
)
plot_monthly_quantiles(d)
```
