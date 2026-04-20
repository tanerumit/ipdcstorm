# Plot daily wind speed with event-duration overlays

Draws a daily wind time series and overlays each event as a horizontal
segment spanning event start-to-end dates at that event's peak wind
speed.

## Usage

``` r
plot_wind_timeseries(
  daily,
  events = NULL,
  thr_tc = 34,
  thr_hur = 64,
  show_thresholds = TRUE,
  title = "Daily Wind with TS/Hurricane Events"
)
```

## Arguments

- daily:

  A data frame or tibble of daily records with columns `date` (`Date`),
  `wind_kt` (numeric, knots), and `event_id` (event identifier; `NA`
  means no event).

- events:

  Optional precomputed event table (data frame/tibble) with columns
  `start_date`, `end_date` (`Date`), `max_wind_kt` (numeric, knots), and
  `event_class` (`"TS"`/`"HUR"`). If `NULL`, it is derived from `daily`.

- thr_tc:

  Numeric scalar; tropical-cyclone threshold in knots.

- thr_hur:

  Numeric scalar; hurricane threshold in knots.

- show_thresholds:

  Logical scalar; if `TRUE`, dashed horizontal threshold lines are
  added.

- title:

  Character scalar used as plot title.

## Value

A `ggplot` object.

## Details

If `events` is `NULL`, event summaries are computed from `daily`. Event
class is derived from each event's peak daily wind (`HUR` for peaks at
or above 64 kt, otherwise `TS`). `NA` values in `event_id` are ignored
when computing events.

## Examples

``` r
d <- data.frame(
  date = as.Date("2001-08-01") + 0:9,
  wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 12, 8),
  event_id = c(NA, NA, 1, 1, 1, NA, NA, NA, NA, NA),
  event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, NA, NA, NA),
  location = "A",
  sim_year = 2001
)
plot_wind_timeseries(d)
```
