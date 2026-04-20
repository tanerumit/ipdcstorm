# Reference observations for wind field validation

Returns a tibble of known station/buoy observations during
well-documented storms for comparison against model-estimated site wind.
Each observation includes a quality tier (`obs_quality`) that determines
the acceptable bias threshold:

- **A** (Direct): Station measurement, instrument survived (\\\pm\\15%).

- **B** (Converted): 10-min reading with averaging-period conversion
  (\\\pm\\25%).

- **C** (Estimated): NHC best-track or indirect estimate (\\\pm\\35%).

## Usage

``` r
get_wind_observations()
```

## Value

Tibble with columns: storm_sid, storm_name, year, target_island,
station, obs_wind_kt, obs_type, obs_quality, obs_source, notes.
