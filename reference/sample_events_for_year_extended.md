# Sample synthetic storm events for a year with extended attributes

Samples synthetic storm events for one calendar year from an event
library and carries forward the atmospheric attributes and identifiers
needed by
[`generate_daily_year_extended()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_year_extended.md).

For each sampled event, extracts from the event library row:

- `event_id`: unique identifier (SID or generated).

- `event_class`: "TS" or "HUR" (for the daily dominant-event tracker).

- `Pc_min_hPa`: minimum central pressure.

- `dP_max_hPa`: maximum pressure deficit.

- `RMW_mean_km`: mean radius of maximum wind.

## Usage

``` r
sample_events_for_year_extended(lib, year, n_ts, n_hur, seed = NULL)
```

## Arguments

- lib:

  List event library produced by
  [`build_event_library()`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library.md).

- year:

  Integer scalar calendar year.

- n_ts:

  Integer scalar number of tropical storms to sample.

- n_hur:

  Integer scalar number of hurricanes to sample.

- seed:

  Optional integer seed for deterministic sampling.

## Value

Tibble with one row per sampled event and event-level attributes used by
the daily downscaling workflow.

## See also

[`generate_daily_year_extended`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_year_extended.md)

## Examples

``` r
sample_lib <- list(
  sample_doy = function(sev) if (sev == "TS") 220L else 250L,
  sample_event = function(sev) {
    tibble::tibble(
      SID = paste0(sev, "_1"),
      dur_days = if (sev == "TS") 3L else 4L,
      V_site_max_kt = if (sev == "TS") 45 else 80,
      Pc_min_hPa = if (sev == "TS") 995 else 970,
      dP_max_hPa = if (sev == "TS") 18 else 40,
      RMW_mean_km = if (sev == "TS") 40 else 25
    )
  }
)
sample_events_for_year_extended(sample_lib, year = 2000, n_ts = 1, n_hur = 1)
#> Error in sample_events_for_year_extended(sample_lib, year = 2000, n_ts = 1,     n_hur = 1): could not find function "sample_events_for_year_extended"
```
