# Generate a parametric wind pulse for a storm event

Creates a deterministic within-event daily wind profile from an event
duration and peak sustained wind.

## Usage

``` r
event_pulse(dur_days, V_peak, shape = c("cosine", "triangle"))
```

## Arguments

- dur_days:

  Integer scalar event duration in days.

- V_peak:

  Numeric scalar peak sustained wind in kt.

- shape:

  Character scalar pulse shape; one of `"cosine"` or `"triangle"`.

## Value

Numeric vector of daily sustained wind values in kt.

## See also

[`generate_daily_year_extended`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_year_extended.md)

## Examples

``` r
event_pulse(dur_days = 5, V_peak = 60)
#> Error in event_pulse(dur_days = 5, V_peak = 60): could not find function "event_pulse"
```
