# Plot event intensity versus duration

Creates a scatter plot of event duration (days) against event peak wind
speed, colored by event class.

## Usage

``` r
plot_intensity_duration(daily = NULL, events = NULL, thr_tc = 34, thr_hur = 64)
```

## Arguments

- daily:

  Optional data frame/tibble of daily rows used only when
  `events = NULL`; must contain columns required by
  [`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md).

- events:

  Optional data frame/tibble with columns `dur_days` (integer days),
  `max_wind_kt` (numeric, knots), and `event_class` (`"TS"`/`"HUR"`).

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

## Value

A `ggplot` object.

## Details

If `events` is `NULL`, `daily` must be supplied and is summarized via
[`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md).
The function errors if both are `NULL`.

## Examples

``` r
d <- data.frame(
  date = as.Date("2001-08-01") + 0:9,
  wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 70, 65),
  event_id = c(NA, NA, 1, 1, 1, NA, NA, NA, 2, 2),
  event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, NA, "HUR", "HUR"),
  location = "A",
  sim_year = 2001
)
plot_intensity_duration(daily = d)
```
