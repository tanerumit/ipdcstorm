# Generate a daily wind + dominant-event attribute series for a single calendar year

Converts sampled storm events for one calendar year into a daily
sustained wind series while tracking the dominant event and key event
attributes for each day.

## Usage

``` r
generate_daily_year_extended(
  year,
  sampled_events,
  pulse_shape = "cosine",
  sim_year = NA_integer_,
  location = NA_character_,
  scenario = NA_character_
)
```

## Arguments

- year:

  Integer scalar calendar year.

- sampled_events:

  Tibble from
  [`sample_events_for_year_extended()`](https://tanerumit.github.io/ipdcstorm/reference/sample_events_for_year_extended.md).

- pulse_shape:

  Character scalar pulse shape passed to
  [`event_pulse()`](https://tanerumit.github.io/ipdcstorm/reference/event_pulse.md).

## Value

Tibble with one row per calendar day and daily hazard attributes.

## See also

[`sample_events_for_year_extended`](https://tanerumit.github.io/ipdcstorm/reference/sample_events_for_year_extended.md),
[`event_pulse`](https://tanerumit.github.io/ipdcstorm/reference/event_pulse.md)

## Examples

``` r
sampled_events <- tibble::tibble(
  start_date = as.Date("2000-08-01"),
  dur_days = 3L,
  V_peak = 60,
  event_id = "evt_1",
  event_class = "TS",
  Pc_min_hPa = 995,
  dP_max_hPa = 18,
  RMW_mean_km = 40
)
generate_daily_year_extended(2000, sampled_events)
#> Error in generate_daily_year_extended(2000, sampled_events): could not find function "generate_daily_year_extended"
```
