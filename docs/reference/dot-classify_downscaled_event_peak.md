# Classify realized downscaled event peaks

Classifies realized downscaled event peaks using the canonical
[`classify_severity()`](https://tanerumit.github.io/ipdcstorm/reference/classify_severity.md)
thresholds when available and otherwise falls back to the local simple
classifier.

## Usage

``` r
.classify_downscaled_event_peak(peak_kt)
```

## Arguments

- peak_kt:

  Numeric vector of realized event peak winds (kt).

## Value

Character vector of storm classes with `NA` for unknown peaks.
