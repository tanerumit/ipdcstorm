# Run the stochastic tropical-cyclone hazard model for one or more target locations

Runs the end-to-end site-level hurricane hazard workflow using cleaned
IBTrACS North Atlantic track data and a user-supplied set of target
locations.

The function reads and filters the IBTrACS archive, selects storm track
points within a specified search radius of each target, computes
location-specific wind exposure metrics from the passing storms,
aggregates those track points into storm-event summaries, estimates
historical annual occurrence rates by storm class, calibrates dispersion
for annual counts, applies climate-change adjustments defined in
`cfg$climate`, and then simulates synthetic annual storm counts for each
location over `cfg$n_sim` years.

This is the main user-level entry point for the hazard model. It returns
both the historical intermediate products used for calibration
(`events`, `trackpoints`, `rates`, `fit`) and the final stochastic
simulation output (`sim`), along with resolved configuration and run
metadata for reproducibility.

## Usage

``` r
run_hazard_model(
  cfg,
  targets,
  storm_classes = c("TS", "HUR"),
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- cfg:

  A `hazard_cfg` object created by
  [`make_hazard_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md).
  This controls the IBTrACS source file, target search radius,
  simulation length, start year for the calibration sample, and advanced
  storm/wind/rate settings.

- targets:

  A data frame or tibble of target locations with required columns:

  `name`

  :   Character location name used as the location key in output tables
      and lists.

  `lat`

  :   Numeric latitude in decimal degrees.

  `lon`

  :   Numeric longitude in decimal degrees.

  Each row defines one site for which a separate hazard calibration and
  stochastic simulation will be performed.

- storm_classes:

  Character vector of storm classes to model. These are normalized
  internally and currently intended to include `"TS"` and `"HUR"`. Only
  the requested classes are carried into the annual-count and rate
  workflow.

- seed:

  Optional integer scalar random seed. If `NULL`, a seed is generated
  internally, set at function entry, and recorded in
  `run_metadata$seed`. All stochastic simulation in the run inherits
  from this seed.

- verbose:

  Logical; if `TRUE` (default), print progress, summary statistics,
  climate settings, and run metadata to the console.

## Value

A named list with the following elements:

- `sim`:

  Tibble of synthetic annual storm counts for all simulated years and
  locations. Contains one row per simulation year and location.

- `events`:

  Tibble of historical storm-event summaries across all locations after
  filtering, event construction, and classification.

- `trackpoints`:

  Named list of tibbles, one per location, containing filtered IBTrACS
  track points and computed site-wind diagnostics.

- `rates`:

  Tibble of calibrated mean annual rates (`lambda`) by location and
  storm class.

- `lambda_scalers`:

  Tibble of location/class-specific lambda scaling factors used to
  adjust the calibrated rate table before simulation.

- `lambda_scaler_id`:

  Character identifier summarizing the applied lambda-scaler
  configuration.

- `fit`:

  Tibble of fitted dispersion and summary parameters by location, with
  climate-related attributes attached.

- `cfg`:

  Resolved hazard configuration used for the run, including normalized
  data path, embedded climate config, and resolved climate metadata.

- `config`:

  Duplicate of `cfg`, returned for compatibility.

- `run_metadata`:

  List with reproducibility metadata including `seed`, `ibtracs_file`,
  `ibtracs_rows`, `ibtracs_data_id`, `parameter_id`, and
  `lambda_scaling_mode`, plus a nested `climate` list with the resolved
  climate diagnostics used for the run.

## Details

The workflow consists of the following major steps:

1.  **IBTrACS data loading and filtering.** The function resolves the
    IBTrACS input path from `cfg$data_path`, reads the cleaned North
    Atlantic archive with
    [`read_ibtracs_clean()`](https://tanerumit.github.io/ipdcstorm/reference/read_ibtracs_clean.md),
    and keeps all available fields needed by downstream wind-field and
    event-processing functions.

2.  **Target-based spatial gating.** For each target location, all track
    points are assigned a great-circle distance to the site using
    [`dist_to_target()`](https://tanerumit.github.io/ipdcstorm/reference/dist_to_target.md).
    Only points within `cfg$search_radius_km` are retained for that
    location. This is a site-centred filtering step: the same storm may
    be retained for one island and excluded for another, depending on
    distance.

3.  **Site wind estimation.** The retained track points are passed to
    [`compute_site_winds_full()`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md),
    which estimates the wind at the target site from storm intensity,
    geometry, and wind-radii information. The returned track-point table
    includes location-specific wind metrics such as asymmetric site
    wind, symmetric site wind, storm maximum wind, and bearing to
    target.

4.  **Storm-event construction and classification.** Filtered track
    points are aggregated into per-storm event summaries using
    [`make_storm_events()`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md).
    Each event is then classified into the requested storm classes with
    [`classify_severity()`](https://tanerumit.github.io/ipdcstorm/reference/classify_severity.md)
    using the thresholds stored in `cfg$advanced`, currently
    tropical-storm and hurricane thresholds in knots.

5.  **Historical annual rates and dispersion calibration.** The event
    catalogue is reduced to annual counts with
    [`compute_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md).
    Mean annual rates by storm class are estimated with
    [`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md),
    and over-dispersion of annual total activity is estimated with
    [`estimate_k_hat()`](https://tanerumit.github.io/ipdcstorm/reference/estimate_k_hat.md).
    These provide the frequency inputs for the stochastic simulation.

6.  **Rate-scaling adjustment.** Before simulation, the function
    computes location/class-specific lambda scalers from the calibrated
    rate table using the configured `lambda_scaling_mode`. This can
    downscale or, depending on mode, upscale location/class rates
    relative to the reference calibration target. The applied scalers
    and their identifier are returned.

7.  **Climate adjustments.** The embedded climate configuration
    `cfg$climate` is resolved with
    [`resolve_climate_inputs()`](https://tanerumit.github.io/ipdcstorm/reference/resolve_climate_inputs.md).
    This produces a scalar `delta_sst`, a frequency sensitivity
    `beta_sst`, a hurricane-fraction sensitivity `gamma`, the applied
    SST-driven count multiplier `rate_scale`, and optional storm
    perturbation settings. Historical climate sensitivity is calibrated
    from basin-consistent annual counts derived from de-duplicated
    storm-year records, independent of the target set. These affect the
    stochastic simulation but do not rewrite the historical event
    catalogue itself.

8.  **Simulation of synthetic annual counts.** For each location, the
    rate table is adjusted, the island-specific base hurricane fraction
    is computed, and
    [`simulate_twolevel_counts()`](https://tanerumit.github.io/ipdcstorm/reference/simulate_twolevel_counts.md)
    is used to generate `cfg$n_sim` synthetic years of annual storm
    activity. The simulation produces counts such as total storms,
    tropical storms, and hurricanes by year and location.

Notes:

- **Historical sample versus synthetic output.** The historical
  component is built from IBTrACS North Atlantic storms gated around
  each site. The synthetic component does not resample full tracks; it
  simulates annual storm counts conditional on the calibrated local rate
  structure.

- **Storm classes.** The model currently uses normalized storm-class
  labels such as `"TS"` and `"HUR"`. These are assigned from site-level
  peak wind at the target location, not basin-wide best-track status
  labels.

- **Climate interpretation.** Baseline runs are represented by
  `make_climate_cfg(scenario = "stationary")` embedded inside `cfg`.
  Baseline and future runs use the same climate-resolution pathway and
  return the same metadata structure; the stationary baseline simply
  resolves `delta_sst = 0`.

- **Reproducibility metadata.** The returned object includes
  `run_metadata` with the seed, IBTrACS file identifier, row count,
  derived data identifier, parameter hash, lambda-scaling mode, and
  resolved climate scalars.

## See also

[`make_hazard_cfg`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md),
[`make_climate_cfg`](https://tanerumit.github.io/ipdcstorm/reference/make_climate_cfg.md),
[`read_ibtracs_clean`](https://tanerumit.github.io/ipdcstorm/reference/read_ibtracs_clean.md),
[`compute_site_winds_full`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md),
[`make_storm_events`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md),
[`compute_annual_counts`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md),
[`compute_lambda_table`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md),
[`estimate_k_hat`](https://tanerumit.github.io/ipdcstorm/reference/estimate_k_hat.md),
[`simulate_twolevel_counts`](https://tanerumit.github.io/ipdcstorm/reference/simulate_twolevel_counts.md)

## Examples

``` r
if (FALSE) { # \dontrun{
targets <- tibble::tribble(
  ~name,        ~lat,     ~lon,
  "Saba",        17.63,  -63.23,
  "St_Martin",   18.07,  -63.05
)

cfg <- make_hazard_cfg(
  data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  historical_start_year = 1980,
  simulation_years = 1000
)

# Historical/stationary baseline run through the common climate resolver
out <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  storm_classes = c("TS", "HUR"),
  seed = 123
)

out$sim
out$rates
out$run_metadata

# Future climate run with embedded climate config
cfg_future <- make_hazard_cfg(
  data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  historical_start_year = 1980,
  simulation_years = 1000,
  climate = make_climate_cfg(scenario = "ssp245", target_year = 2050)
)

out_climate <- run_hazard_model(
  cfg = cfg_future,
  targets = targets,
  seed = 123
)
} # }
```
