# Plot monthly counts of event starts

Summarizes event starts by calendar month and plots grouped bars by
event class, optionally normalized to an annual rate.

## Usage

``` r
plot_monthly_events(daily, normalize = FALSE)
```

## Arguments

- daily:

  A data frame or tibble with columns required by
  [`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md)
  (`sim_year`, `event_id`, `date`, `wind_kt`, and optional `location`).
  Rows with `NA` `event_id` are excluded during event extraction.

- normalize:

  Logical scalar; if `TRUE`, plot events per year.

## Value

A `ggplot` object.

## Details

Event starts are computed from
[`prep_events()`](https://tanerumit.github.io/ipdcstorm/reference/prep_events.md).
Missing month/class combinations are filled with zero. When
`normalize = TRUE`, counts are divided by `n_distinct(daily$sim_year)`.

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
plot_monthly_events(d)
```
