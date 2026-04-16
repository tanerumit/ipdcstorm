# Plot day-of-year distribution of event timing

Plots a day-of-year histogram for either event-active days or event
starts, optionally faceted by event class.

## Usage

``` r
plot_seasonality_doy(
  daily,
  metric = c("event_days", "starts"),
  facet_class = TRUE,
  binwidth = 7
)
```

## Arguments

- daily:

  A data frame or tibble with at least `date` (`Date`), `event_id`,
  `event_class`, plus fields required by
  [`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md)
  when `metric = "starts"` (`location`, `sim_year`, `wind_kt`).

- metric:

  Character scalar; one of `"event_days"` or `"starts"`.

- facet_class:

  Logical scalar; if `TRUE`, create one panel per class.

- binwidth:

  Numeric scalar \> 0; histogram bin width in day-of-year units (days).

## Value

A `ggplot` object.

## Details

For `metric = "event_days"`, each day with non-`NA` `event_id`
contributes one count. For `metric = "starts"`, counts are based on
event start dates derived from grouped events. Event classes are
collapsed to `TS`/`HUR` using the same rule as
[`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md).

## Examples

``` r
d <- data.frame(
  date = as.Date("2001-08-01") + 0:9,
  wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 12, 8),
  event_id = c(NA, NA, 1, 1, 1, NA, NA, 2, 2, NA),
  event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, "HUR", "HUR", NA),
  location = "A",
  sim_year = 2001
)
plot_seasonality_doy(d, metric = "event_days")
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_bar()`).
```
