# Build an empirical event library for resampling

Builds a resampling library from historical storm events for one target
location. The library stores empirical day-of-year samples by storm
class, stratified historical event bins, and sampler closures used by
the daily downscaling workflow.

## Usage

``` r
build_event_library(
  track_df,
  event_df,
  storm_classes = c("TD", "TS", "HUR"),
  bins = list(wind = c(0, 34, 64, 83, 96, 113, Inf), Pc = c(850, 900, 940, 970, 1000,
    1050), RMW = c(0, 10, 20, 30, 40, 60, Inf)),
  seed = NULL,
  resampling_method = c("stratified", "copula_nn"),
  copula_min_n = 30L,
  copula_k = 1L,
  copula_robust_scale = TRUE
)
```

## Arguments

- track_df:

  Tibble/data.frame of track points for one location; must contain
  columns `SID` and `iso_time`.

- event_df:

  Tibble/data.frame of storm events with one row per storm; must contain
  `SID` or `storm_id`.

- storm_classes:

  Character vector of storm classes retained in the library.

- bins:

  Named list of numeric break vectors for wind, pressure, and RMW
  stratification.

- seed:

  Optional integer seed for deterministic library construction.

- resampling_method:

  Character scalar; one of `"stratified"` or `"copula_nn"`.

- copula_min_n:

  Integer scalar; minimum complete-case sample size needed to fit a
  class-specific copula sampler.

- copula_k:

  Integer scalar; number of nearest neighbours sampled from a copula
  proposal.

- copula_robust_scale:

  Logical scalar; if `TRUE`, use median/MAD scaling in nearest-neighbour
  distance calculations.

## Value

List with empirical day-of-year samples, stratification bins, the
processed event table, and sampler functions.

## See also

[`build_event_library_from_out`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library_from_out.md)

## Examples

``` r
track_df <- tibble::tibble(
  SID = c("A", "A", "B", "B"),
  iso_time = as.POSIXct(c(
    "2000-08-01 00:00:00", "2000-08-01 06:00:00",
    "2000-09-10 00:00:00", "2000-09-10 06:00:00"
  ), tz = "UTC")
)
event_df <- tibble::tibble(
  SID = c("A", "B"),
  peak_wind_kt = c(45, 80),
  storm_intensity_kt = c(45, 80),
  min_pressure_hpa = c(995, 970),
  rmw_mean_km = c(40, 25),
  start_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-09-10 00:00:00"), tz = "UTC")
)
lib <- build_event_library(track_df, event_df, storm_classes = c("TS", "HUR"))
names(lib)
#> [1] "doy_by_severity" "strat_bins"      "events"          "sample_doy"     
#> [5] "sample_event"   
```
