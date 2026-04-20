# Create a hazard model configuration

Creates a user-facing configuration object for
[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md).
Most users only need to set the IBTrACS path and simulation horizon.

The `preset = "default"` setting applies standard threshold and cap
values. Expert tuning knobs can be provided via `advanced`, but are
optional.

## Usage

``` r
make_hazard_cfg(
  data_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
  search_radius_km = 800,
  historical_start_year = 1970L,
  simulation_years = 1000L,
  preset = "default",
  climate = make_climate_cfg(scenario = "stationary"),
  advanced = NULL,
  background_wind = NULL
)
```

## Arguments

- data_path:

  Character; path to IBTrACS CSV input data.

- search_radius_km:

  Numeric; maximum distance from each target used to include track
  points.

- historical_start_year:

  Integer; first historical year used to fit rates.

- simulation_years:

  Integer; number of synthetic years to simulate.

- preset:

  Character; currently only `"default"`.

- climate:

  Climate configuration object created by
  [`make_climate_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_climate_cfg.md).
  Defaults to the canonical stationary baseline
  `make_climate_cfg(scenario = "stationary")`.

- advanced:

  Optional named list of expert parameters. Most users should leave this
  as `NULL`. Supported names: `ts_threshold_kt`,
  `strong_storm_threshold_kt`, `hurricane_threshold_kt`, `r34_cap_nm`,
  `r50_cap_nm`, `r64_cap_nm`, `lambda_scaling_mode` (`"target"` or
  `"down_only"`).

- background_wind:

  Optional `background_wind_cfg` object created by
  [`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md).
  When supplied, non-storm days in the daily downscaling step receive
  wind speeds sampled from Weibull marginals with optional spatial
  copula correlation and AR(1) persistence, rather than being set to
  zero. Defaults to `NULL` (zero background wind).

## Value

A list with class `c("hazard_cfg", "list")`.
