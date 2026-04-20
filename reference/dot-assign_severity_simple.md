# Assign simple severity class from peak wind

Applies fixed tropical-depression, tropical-storm, and hurricane wind
thresholds to a realized peak wind value.

## Usage

``` r
.assign_severity_simple(wind_max_kt)
```

## Arguments

- wind_max_kt:

  Numeric; peak wind (kt).

## Value

Character scalar in `c("TD", "TS", "HUR")` (or `NA`).
