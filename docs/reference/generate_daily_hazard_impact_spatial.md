# Generate spatially coherent daily synthetic hazard and impact time series

Drop-in replacement for `generate_daily_hazard_impact_spatial()` that
enforces spatial coherence across locations by sampling storms once at
the basin level per simulated year and assigning each drawn storm to
every location whose event library contains it (Option 1 — shared event
pool).

The key difference from the independent-sampling baseline is that the
annual storm draw happens *before* the location loop rather than inside
it. Consequently, if Hurricane Irma (`"2017242N16333"`) is drawn in
simulated year 47, it will appear simultaneously at every location that
has an Irma entry in its event library (e.g. both St. Martin and Saba),
each retaining its own site-level wind profile, duration, and pressure
attributes. Storms whose tracks never came within the search radius of a
location are absent from that location's library and are therefore
skipped silently for that location, so per-location event counts may be
lower than the basin-level draw.

Basin-level annual counts are resolved as `max(n_ts)` and `max(n_hur)`
across all requested locations for each simulation year, using the
counts already present in `out$sim`.

All downstream processing (climate perturbation, pulse generation,
damage forcing) is identical to the independent-sampling variant.

## Usage

``` r
generate_daily_hazard_impact_spatial(
  out,
  location,
  sim_years = 1:1000,
  year0 = 2000,
  gust_factor = 1,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = NA_character_,
  seed = NULL,
  pinned_sids = NULL,
  background_wind = NULL,
  verbose = FALSE
)
```

## Arguments

- out:

  List returned by
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).

- location:

  Character vector of one or more target location names.

- sim_years:

  Integer vector of simulation-year indices to generate.

- year0:

  Integer scalar base calendar year for `sim_year == 1`.

- gust_factor:

  Numeric scalar gust multiplier applied to `wind_kt`.

- damage:

  Named list defining the daily damage model; same specification as
  `generate_daily_hazard_impact_spatial()`.

- pulse_shape:

  Character scalar pulse shape; `"cosine"` or `"triangle"`.

- scenario:

  Optional character scalar scenario label carried into output.

- seed:

  Optional integer scalar seed. Defaults to `out$run_metadata$seed` or
  `1L`. Per-location library seeds are offset by location index for
  reproducibility.

## Value

Named list of tibbles — one per requested location — with the same
column schema as `generate_daily_hazard_impact_spatial()`.

## See also

`generate_daily_hazard_impact_spatial`,
[`build_event_library_from_out`](https://tanerumit.github.io/ipdcstorm/reference/build_event_library_from_out.md),
[`generate_daily_year_extended`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_year_extended.md)

## Examples

``` r
if (FALSE) { # \dontrun{
daily_spatial <- generate_daily_hazard_impact_spatial(
  out         = hazard_out,
  location    = c("Saba", "St_Martin", "Statia"),
  sim_years   = 1:2000,
  year0       = 2025,
  gust_factor = 1.3,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "stationary"
)
# Verify spatial coherence: Irma should appear in the same sim_years
# at every location that has Irma in its event library.
lapply(daily_spatial, function(loc_tbl) {
  loc_tbl |>
    dplyr::filter(grepl("^2017242N16333", event_id)) |>
    dplyr::distinct(sim_year) |>
    dplyr::pull(sim_year)
})
} # }
```
