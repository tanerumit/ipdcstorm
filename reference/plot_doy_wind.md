# Plot mean and maximum wind by day of year

Aggregates daily wind by day-of-year and plots mean and maximum wind
speed curves with optional loess smoothing.

## Usage

``` r
plot_doy_wind(daily, smooth = TRUE, span = 0.15, thr_tc = 34, thr_hur = 64)
```

## Arguments

- daily:

  A data frame or tibble with columns `date` (`Date`) and `wind_kt`
  (numeric, knots).

- smooth:

  Logical scalar; if `TRUE`, use `geom_smooth(method = "loess")`,
  otherwise plot unsmoothed lines.

- span:

  Numeric scalar in `(0, 1]`; loess span used when `smooth = TRUE`.

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

## Value

A `ggplot` object.

## Details

Summary statistics use `na.rm = TRUE`; all-`NA` day groups propagate
non-finite summaries. Threshold lines are always drawn at `thr_tc` and
`thr_hur`.

## Examples

``` r
d <- data.frame(
  date = as.Date("2001-01-01") + 0:29,
  wind_kt = 15 + sin((1:30) / 5) * 10
)
plot_doy_wind(d, smooth = FALSE)
```
