# Plot distribution of daily wind speed

Plots daily wind speed as either a histogram or kernel density estimate
with TS/HUR threshold reference lines.

## Usage

``` r
plot_wind_distribution(
  daily,
  type = c("histogram", "density"),
  log_scale = FALSE,
  thr_tc = 34,
  thr_hur = 64
)
```

## Arguments

- daily:

  A data frame or tibble with column `wind_kt` (numeric, knots; `NA`
  values are ignored by ggplot stat layers).

- type:

  Character scalar; one of `"histogram"` or `"density"`.

- log_scale:

  Logical scalar; if `TRUE`, use a log10 x-axis.

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

## Value

A `ggplot` object.

## Details

`type` is matched with
[`match.arg()`](https://rdrr.io/r/base/match.arg.html). If
`log_scale = TRUE`,
[`scale_x_log10()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
is applied; non-positive `wind_kt` values are therefore not displayable
on the transformed axis.

## Examples

``` r
d <- data.frame(wind_kt = c(5, 10, 15, 20, 35, 40, 65, 70))
plot_wind_distribution(d, type = "histogram")
```
