# Classify storm class from peak site wind.

Classifies peak site wind into tropical depression, tropical storm,
hurricane, or unknown.

## Usage

``` r
classify_severity(
  V_site_max_kt,
  ts_threshold_kt = 34,
  hurricane_threshold_kt = 64
)
```

## Arguments

- V_site_max_kt:

  Numeric vector of peak site wind speeds in knots.

- ts_threshold_kt:

  Numeric scalar tropical-storm threshold in knots.

- hurricane_threshold_kt:

  Numeric scalar hurricane threshold in knots.

## Value

Character vector containing `"TD"`, `"TS"`, `"HUR"`, or `"unknown"`.

## See also

[`make_storm_events`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md),
[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md)

## Examples

``` r
classify_severity(c(20, 40, 80, NA_real_))
#> [1] "TD"      "TS"      "HUR"     "unknown"
```
