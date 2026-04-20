# Summarize event-level features from daily rows

Summarize event-level features from daily rows

## Usage

``` r
prep_events(daily)
```

## Arguments

- daily:

  A data frame or tibble containing daily rows with columns `sim_year`,
  `event_id`, `date`, and `wind_kt`. A `location` column is optional for
  single-location inputs and is reconstructed from tibble metadata when
  absent. Rows with `NA` in `event_id` are excluded.

## Value

A tibble with one row per unique (`location`, `sim_year`, `event_id`)
and columns `event_class`, `start_date`, `end_date`, `dur_days`,
`max_wind_kt`, `start_doy`, and `start_month`.
