# Tutorial 1 — Introduction, Setup & Baseline Workflow

## 1. Overview

`ipdcstorm` is a parametric tropical-cyclone wind hazard model. It
combines NOAA’s IBTrACS historical track database with stochastic
simulation to generate long synthetic records of daily site-level wind
exposure. The model was developed for climate stress testing of
supply-chain infrastructure on the Dutch Caribbean islands (St. Martin,
Saba, Statia), but can be adapted for any point location in the North Atlantic
basin.

The workflow has three stages:

1.  **Historical calibration.** Read IBTrACS, filter to track points
    within `search_radius_km` of each site, and compute site-level winds
    with a Holland (1980) parametric profile that uses directional
    R34/R50/R64 wind radii (with Knaff-style climatological fallbacks)
    and a forward-motion asymmetry term. Site events are aggregated and
    classified by peak site wind into TS (≥ 34 kt) and HUR (≥ 64 kt).
    Annual class rates and a Gamma-Poisson overdispersion parameter
    `k_hat` are estimated from the historical record.
2.  **Stochastic simulation.** Draw synthetic seasons from a two-level
    process: a Gamma-distributed annual activity multiplier `A`
    reproduces the observed burstiness of seasons, then total counts and
    the TS / HUR split follow Poisson and binomial draws conditional on
    `A`. Climate forcing — when enabled in Tutorial 3 — enters via
    SST-driven changes to annual rates and the hurricane fraction, with
    optional event-level perturbation of intensity, size, and
    translation speed.
3.  **Temporal downscaling.** For each simulated year, sample events
    from a severity-stratified empirical event library, place them on
    the calendar with empirical day-of-year distributions, spread each
    event over its duration with a smooth pulse, and combine overlapping
    events into a daily wind (kt) series. A damage proxy and supporting
    daily fields (gust, surge, dominant-event metadata) are computed
    from the daily wind.

## 2. Prerequisites

### 2.1. Installation

``` r
# install.packages("devtools")
# devtools::install_github("tanerumit/ipdcstorm")
```

``` r
library(ipdcstorm)
library(dplyr)
```

### 2.2. IBTrACS data

Calibration uses NOAA’s International Best Track Archive for Climate
Stewardship (IBTrACS). The model uses the **North Atlantic subset**
(`ibtracs.NA.list.v04r01.csv`), with 6-hourly track points for every
recorded Atlantic storm.

| IBTrACS field(s) | Used for |
|----|----|
| `SID`, `SEASON`, `BASIN` | Storm identification, year, basin filter |
| `USA_LAT`, `USA_LON`, `ISO_TIME` | 6-hourly track position |
| `USA_WIND` | Maximum sustained wind (kt) |
| `USA_R34_NE/SE/SW/NW`, `USA_R50_*`, `USA_R64_*` | Directional wind radii (nm) at 34, 50, 64 kt thresholds |
| `USA_PRES`, `USA_POCI` | Central pressure and outermost closed-isobar pressure |
| `USA_RMW`, `USA_ROCI` | Radius of maximum wind and outermost closed isobar |

Wind radii drive the spatial extent of the wind field. When they are
missing (mostly pre-2004), the model substitutes Knaff et al. (2015)
climatological estimates.

**Download:** [IBTrACS
NCEI](https://www.ncei.noaa.gov/products/international-best-track-archive)


## 3. Configuration

``` r
SEED             <- 42L
N_SIM            <- 2000L
HIST_START_YEAR  <- 1970L
SEARCH_RADIUS_KM <- 800

# Calendar year mapped to sim_year = 1 in the daily output.
YEAR0 <- 2000L

GUST_FACTOR <- 1.0
# 1.0 → 1-min sustained wind (use this for validation against IBTrACS).
# 1.3 → standard WMO gust factor (use for damage / infrastructure impact).
```

``` r
targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740
)
```

``` r
DATA_PATH <- system.file("extdata", "ibtracs_1970.csv", package = "ipdcstorm")

hazard_cfg <- make_hazard_cfg(
  data_path             = DATA_PATH,
  search_radius_km      = SEARCH_RADIUS_KM,
  historical_start_year = HIST_START_YEAR,
  simulation_years      = N_SIM,
  climate               = make_climate_cfg(scenario = "stationary")
)

hazard_cfg
```

    #> Hazard configuration (preset: "default")
    #>   IBTrACS data  : C:/Users/taner/AppData/Local/R/win-library/4.5/ipdcstorm/extdata/ibtracs_1970.csv
    #>   Study period  : 1970 - present
    #>   Search radius : 800 km
    #>   Simulation    : 2,000 synthetic years
    #>   Thresholds    : WMO standard (TS >= 34 kt, Hurricane >= 64 kt) [preset]
    #>   Lambda scaling: target

> **What the configuration parameters do**
>
> | Parameter | Role |
> |----|----|
> | `data_path` | Path to the IBTrACS NA CSV (or packaged subset). |
> | `historical_start_year` | First IBTrACS year used to fit storm rates. 1970 marks the start of reliable satellite-era observations in the Atlantic. |
> | `search_radius_km` | Track-point gate radius per site. Keep at 800 km unless you have a specific reason to deviate; a storm 700 km away can still produce TS-strength winds at the site through its extended wind field. |
> | `simulation_years` | Number of synthetic years to generate. |
> | `climate` | A `climate_cfg` object. `make_climate_cfg(scenario = "stationary")` switches off all SST forcing. |
> | `preset` / `advanced` | Hooks for tuning class thresholds, wind-radii caps, and the lambda-scaling mode. Defaults follow WMO conventions (TS ≥ 34 kt, HUR ≥ 64 kt). |
> | `background_wind` | Optional `background_wind_cfg` object (see [Section 5.3](#sec-background-wind)). Fills non-storm days with trade-wind speeds sampled from site-specific Weibull distributions. `NULL` (default) sets all non-storm days to 0 kt. |
>
> See
> [`?make_hazard_cfg`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md)
> for the full parameter list.

## 4. Model Calibration and Simulation

[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
executes the historical calibration and stochastic simulation
(Stage 2) in one call. With a stationary `climate` config, no SST
forcing is applied.

``` r
out_base <- run_hazard_model(
  cfg     = hazard_cfg,
  targets = targets,
  seed    = SEED,
  verbose = FALSE
)
```

### 4.1. What happens inside `run_hazard_model()`

**Step 1 — Read and clean IBTrACS**
([`read_ibtracs_clean()`](https://tanerumit.github.io/ipdcstorm/reference/read_ibtracs_clean.md)).
Standardize columns, retain the USA best-track wind/pressure fields,
convert wind radii to km, and filter to the North Atlantic.

**Step 2 — Gate track points.** Compute great-circle distance from every
6-hour track point to each site and keep only points within
`search_radius_km`.

**Step 3 — Compute site-level winds**
([`compute_site_winds_full()`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md)).
Determine the site’s quadrant relative to the storm center, pull the
matching directional wind radii (R34 / R50 / R64), estimate RMW
([`estimate_RMW_knaff()`](https://tanerumit.github.io/ipdcstorm/reference/estimate_RMW_knaff.md)),
evaluate the parametric wind profile, and apply a forward-motion
asymmetry adjustment. Output: `site_wind_kt` per track point.

**Step 4 — Aggregate into events**
([`make_storm_events()`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md)).
Group track points by storm and summarise into one event per (storm ×
site) with peak site wind, intensity and pressure summaries, mean RMW,
timing, and the gated-point count. Classify by peak site wind into TS
(34–63 kt) or HUR (≥ 64 kt).

**Step 5 — Estimate annual rates**
([`compute_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md)
→
[`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)).
Convert per-class annual counts into mean rates `lambda` and estimate
the Gamma overdispersion `k_hat` that controls year-to-year clustering.

**Step 6 — Simulate annual counts**
([`simulate_twolevel_counts()`](https://tanerumit.github.io/ipdcstorm/reference/simulate_twolevel_counts.md)).
For `simulation_years` synthetic years, draw the activity multiplier `A`
from `Gamma(k_hat, k_hat)`, then total impactful storm counts from a
Poisson mean conditional on `A`, then split into TS vs HUR via a
binomial draw at `p_hurricane_base` (in stationary mode).

> **Inspect the calibrated rates first**
>
> The first thing to look at after a run is `out_base$rates`. If a
> site’s `lambda` is far from a published rate, the rest of the
> workflow will inherit that bias — fix it before downstream analysis.

``` r
out_base$rates |>
  select(location, storm_class, lambda) |>
  mutate(lambda = round(lambda, 3)) |>
  arrange(location, storm_class)
```

    #> # A tibble: 6 × 3
    #>   location  storm_class lambda
    #>   <chr>     <chr>        <dbl>
    #> 1 Saba      HUR          0.064
    #> 2 Saba      TS           0.489
    #> 3 St_Martin HUR          0.125
    #> 4 St_Martin TS           0.393
    #> 5 Statia    HUR          0.087
    #> 6 Statia    TS           0.5

### 4.2. Inspecting `out_base`

The return value is a named list with one tidy table per stage.

#### 4.2.1. Rates

``` r
out_base$rates
```

    #> # A tibble: 6 × 6
    #>   location  storm_class lambda n_years prob_annual prob_none
    #>   <chr>     <chr>        <dbl>   <int>       <dbl>     <dbl>
    #> 1 St_Martin HUR         0.125       56      0.118      0.882
    #> 2 St_Martin TS          0.393       56      0.325      0.675
    #> 3 Saba      HUR         0.0638      47      0.0618     0.938
    #> 4 Saba      TS          0.489       47      0.387      0.613
    #> 5 Statia    HUR         0.0870      46      0.0833     0.917
    #> 6 Statia    TS          0.5         46      0.393      0.607

| Column        | Description                                              |
|---------------|----------------------------------------------------------|
| `location`    | Target site.                                             |
| `storm_class` | `TS` (34–63 kt site wind) or `HUR` (≥ 64 kt).            |
| `lambda`      | Mean annual rate of class events.                        |
| `n_years`     | Historical years contributing to the estimate.           |
| `prob_annual` | Probability of ≥ 1 event per year, $`1 - e^{-\lambda}`$. |
| `prob_none`   | Probability of zero events, $`e^{-\lambda}`$.            |

A `lambda` of 1.2 for TS means roughly 1.2 tropical-storm impacts per
year on average; the actual count varies because of the
Gamma-distributed activity multiplier.

#### 4.2.2. Fit

``` r
out_base$fit
```

    #> # A tibble: 3 × 7
    #>   location     k_hat mean_annual_total var_annual_total beta_sst gamma_intensity
    #>   <chr>        <dbl>             <dbl>            <dbl>    <dbl>           <dbl>
    #> 1 St_Martin   1   e6             0.518            0.509      0.6           0.244
    #> 2 Saba        8.27e1             0.553            0.557      0.6           0.244
    #> 3 Statia      2.10e1             0.587            0.603      0.6           0.244
    #> # ℹ 1 more variable: p_hurricane_base <dbl>

| Column | Description |
|----|----|
| `location` | Target site. |
| `k_hat` | Gamma-Poisson overdispersion. Smaller values (~ 2–5) imply burstier seasons; larger values approach pure Poisson. |
| `mean_annual_total` / `var_annual_total` | First two moments of the total annual count. |
| `beta_sst` | SST rate-scaling coefficient (0 in stationary mode). |
| `gamma_intensity` | Intensity / class-share shift coefficient (0 in stationary mode). |
| `p_hurricane_base` | Baseline probability that an impactful event is HUR-class. |

#### 4.2.3. Sim

``` r
head(out_base$sim)
```

    #> # A tibble: 6 × 9
    #>   location  sim_year activity_factor climate_scale activity_combined p_hurricane
    #>   <chr>        <int>           <dbl>         <dbl>             <dbl>       <dbl>
    #> 1 St_Martin        1           1.00              1             1.00        0.190
    #> 2 St_Martin        2           0.999             1             0.999       0.190
    #> 3 St_Martin        3           1.00              1             1.00        0.190
    #> 4 St_Martin        4           0.999             1             0.999       0.190
    #> 5 St_Martin        5           1.000             1             1.000       0.190
    #> 6 St_Martin        6           0.999             1             0.999       0.190
    #> # ℹ 3 more variables: n_total <int>, n_ts <int>, n_hur <int>

| Column | Description |
|----|----|
| `location` | Target site. |
| `sim_year` | Synthetic year index, 1 to `N_SIM`. |
| `n_total` / `n_tc` | Total impactful storms in the year. |
| `n_ts` | TS-class events. |
| `n_hur` | HUR-class events. |
| `A` | Activity multiplier from `Gamma(k_hat, k_hat)`; \> 1 = hyperactive season, \< 1 = quiet. |
| `p_hur` | Hurricane share applied (constant in stationary mode). |

This table drives the daily downscaling and all downstream return-period
and exceedance work.

#### 4.2.4. Events

``` r
head(out_base$events)
```

    #> # A tibble: 6 × 13
    #>   location  storm_id      start_time          end_time            n_points
    #>   <chr>     <chr>         <dttm>              <dttm>                 <int>
    #> 1 St_Martin 1970220N14342 1970-08-13 03:00:00 1970-08-14 15:00:00       12
    #> 2 St_Martin 1970223N19340 1970-08-17 03:00:00 1970-08-18 18:00:00       13
    #> 3 St_Martin 1970230N12316 1970-08-20 03:00:00 1970-08-21 21:00:00       15
    #> 4 St_Martin 1970247N17316 1970-09-10 06:00:00 1970-09-13 09:00:00       23
    #> 5 St_Martin 1970274N15306 1970-10-01 03:00:00 1970-10-10 03:00:00       52
    #> 6 St_Martin 1971231N13308 1971-08-19 03:00:00 1971-08-20 21:00:00       14
    #> # ℹ 8 more variables: peak_wind_kt <dbl>, storm_intensity_kt <dbl>,
    #> #   min_pressure_hpa <dbl>, pressure_deficit_hpa <dbl>, rmw_min_km <dbl>,
    #> #   rmw_mean_km <dbl>, year <dbl>, storm_class <chr>

One row per (historical storm × site) with site-attenuated wind:

| Column | Description |
|----|----|
| `storm_id` | IBTrACS SID (e.g. `2017242N16333` for Irma). |
| `year` | Season year. |
| `peak_wind_kt` | **Peak wind at the site** — the model-computed value, accounting for distance and wind-field structure. |
| `storm_intensity_kt` | The storm’s lifetime peak intensity, regardless of location. |
| `min_pressure_hpa` / `pressure_deficit_hpa` | Central pressure and maximum environment-minus-centre deficit. |
| `rmw_mean_km` | Mean radius of maximum wind. |
| `storm_class` | Site classification (`TS` or `HUR`). |
| `location` | Target site. |

The distinction between `peak_wind_kt` and `storm_intensity_kt` matters:
a Cat-5 storm (≈ 155 kt) passing 400 km from a site might only produce
Cat-1 winds (~ 70 kt) there. The Holland profile captures that
attenuation.

#### 4.2.5. Trackpoints

``` r
head(out_base$trackpoints[["Saba"]])
```

    #> # A tibble: 6 × 54
    #>   SID       SEASON BASIN iso_time            storm_name storm_status   lat   lon
    #>   <chr>      <dbl> <chr> <dttm>              <chr>      <chr>        <dbl> <dbl>
    #> 1 1970220N…   1970 NA    1970-08-13 03:00:00 UNNAMED    TS            14.7 -57.3
    #> 2 1970220N…   1970 NA    1970-08-13 06:00:00 UNNAMED    TS            15.3 -58  
    #> 3 1970220N…   1970 NA    1970-08-13 09:00:00 UNNAMED    TS            16   -58.7
    #> 4 1970220N…   1970 NA    1970-08-13 12:00:00 UNNAMED    TS            16.7 -59.5
    #> 5 1970220N…   1970 NA    1970-08-13 15:00:00 UNNAMED    TS            17.5 -60.4
    #> 6 1970220N…   1970 NA    1970-08-13 18:00:00 UNNAMED    TS            18.2 -61.3
    #> # ℹ 46 more variables: wind_kt <dbl>, pres_hpa <dbl>, pres_source <chr>,
    #> #   usa_pres_hpa <dbl>, wmo_pres_hpa <dbl>, poci_hpa <dbl>, roci_km <dbl>,
    #> #   rmw_km <dbl>, r34_ne_nm <dbl>, r34_se_nm <dbl>, r34_sw_nm <dbl>,
    #> #   r34_nw_nm <dbl>, r50_ne_nm <dbl>, r50_se_nm <dbl>, r50_sw_nm <dbl>,
    #> #   r50_nw_nm <dbl>, r64_ne_nm <dbl>, r64_se_nm <dbl>, r64_sw_nm <dbl>,
    #> #   r64_nw_nm <dbl>, storm_speed_kt <dbl>, storm_dir_deg <dbl>, dist_km <dbl>,
    #> #   bearing_deg <dbl>, quadrant <chr>, nq34 <dbl>, R34_nm_dir <dbl>, …

The lowest-level output: every gated 6-hourly IBTrACS track point per
site, with all computed wind-field columns (`site_wind_kt`, `dist_km`,
`bearing_deg`, `RMW_km`, etc.). Useful for diagnostics and spot-checking
individual storms.

## 5. Daily Wind Time Series

[`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)
converts the per-year storm counts into a daily site-level wind series
for every requested simulated year and location.

``` r
daily_by_loc <- generate_daily_hazard_impact_spatial(
  out         = out_base,
  location    = targets$name,
  sim_years   = seq_len(N_SIM),
  year0       = YEAR0,
  gust_factor = GUST_FACTOR,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "stationary",
  seed        = SEED
)
```

The returned object is a named list of daily tibbles — one per location.

### 5.1. What the function does

**Step 1 — Build event libraries.** For each location, the model builds
a stratified empirical library of historical events with their physical
attributes paired (peak wind, duration, pressure, RMW). A shared
basin-level event pool is constructed across locations to keep
multi-site sampling consistent.

**Step 2 — Sample events for each simulated year.** For each `sim_year`,
draw the prescribed number of TS and HUR events from the library.
Attribute pairing is preserved so wind, duration, pressure, and RMW
remain physically consistent.

**Step 3 — Place events on the calendar.** Each sampled event gets a
start date drawn from the empirical day-of-year distribution,
concentrating events in the August–October peak.

**Step 4 — Spread each event into a daily wind pulse.** The `cosine`
pulse shape gives a smooth bell-shaped envelope rising to the event peak
and decaying over the event duration. `pulse_shape = "uniform"` is also
available for a flat box.

**Step 5 — Combine overlapping events.** Multiple events that overlap in
calendar time are reduced to the daily maximum sustained wind. This
avoids double-counting while still capturing compound seasons.

**Step 6 — Derived daily fields.** Apply the damage function
(`damage = list(method = "intensity")` returns the normalized
intensity-based damage proxy; `"powerlaw"` is also available), estimate
`surge_m`, and forward dominant-event metadata (`event_id`,
`pressure_hpa`, `rmw_km`).

> **Common-random-numbers comparisons**
>
> Pass the same `seed` and the same `sim_years` vector across baseline
> and scenario runs, and any difference in daily output isolates the
> climate forcing. 

> **Pinned events and per-year jitter**
>
> `pinned_sids` is a named list (keys = `sim_year` as character) that
> forces a specific historical SID to keep a focal storm
> (e.g. Irma) in every synthetic year. 

### 5.2. Daily output schema

``` r
head(daily_by_loc[["Saba"]])
```

    #> # A tibble: 6 × 11
    #>   sim_year   doy wind_kt surge_m event_id pressure_hpa pressure_deficit_hpa
    #>      <int> <int>   <dbl>   <dbl> <chr>           <dbl>                <dbl>
    #> 1        1     1       0      NA <NA>               NA                   NA
    #> 2        1     2       0      NA <NA>               NA                   NA
    #> 3        1     3       0      NA <NA>               NA                   NA
    #> 4        1     4       0      NA <NA>               NA                   NA
    #> 5        1     5       0      NA <NA>               NA                   NA
    #> 6        1     6       0      NA <NA>               NA                   NA
    #> # ℹ 4 more variables: rmw_km <dbl>, damage_intensity <dbl>, damage_rate <dbl>,
    #> #   cum_damage <dbl>

| Column | Description |
|----|----|
| `sim_year` | Simulated year index (serial 1..N). |
| `doy` | Day of year, integer 1–366. Synthetic `sim_year` values are serial and not tied to real calendar years; `doy` is the only temporal axis in daily output. |
| `wind_kt` | Daily site-level sustained wind (kt). |
| `surge_m` | Estimated storm surge (m), derived from event properties. |
| `event_id` | Identifier of the dominant event on that day. |
| `pressure_hpa` / `pressure_deficit_hpa` | Central pressure and deficit of the dominant event. |
| `rmw_km` | Radius of maximum wind of the dominant event. |
| `damage_intensity` | Normalised hazard intensity on a 0–1 scale. |
| `damage_rate` | Daily fractional damage rate from the damage function. |
| `cum_damage` | Cumulative damage within each simulated year. |

The location name is stored as the list-element name and as a tibble
attribute; it is not a column of the returned tibble. To build a flat
long tibble for cross-location analysis, use
`dplyr::bind_rows(daily_by_loc, .id = "location")`.

### 5.3. Background Wind

By default, every day without a TC event has `wind_kt = 0`. This is
intentional for hazard-return-period work (zeros represent TC-free days
cleanly), but it makes the daily series unsuitable for applications that
need a realistic non-storm wind envelope.

[`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
specifies Weibull marginal distributions for the trade-wind background,
with an optional Gaussian copula for spatial correlation across sites
and an AR(1) coefficient for day-to-day persistence. On each simulated
day, a background wind speed is drawn from the site- and month-specific
Weibull distribution.

> **TC-day winds are unaffected**
>
> Background wind is applied via `pmax(bg, tc_pulse)`. On days when the
> TC pulse exceeds the background draw, the TC value is kept unchanged.
> Return-level statistics for TS and HUR thresholds are therefore
> unaffected by the background wind setting.

## 6. Annual Activity Summary

Collapse the daily series into per-site annual statistics.

``` r
daily_all <- bind_rows(daily_by_loc, .id = "location")

annual_summary <- daily_all |>
  group_by(location, sim_year) |>
  summarise(
    peak_wind_kt  = max(wind_kt,    na.rm = TRUE),
    n_ts_days     = sum(wind_kt >= 34, na.rm = TRUE),
    n_hur_days    = sum(wind_kt >= 64, na.rm = TRUE),
    annual_damage = sum(damage_rate, na.rm = TRUE),
    .groups = "drop"
  )

activity_means <- annual_summary |>
  group_by(location) |>
  summarise(
    mean_peak_wind_kt  = mean(peak_wind_kt),
    mean_ts_days       = mean(n_ts_days),
    mean_hur_days      = mean(n_hur_days),
    mean_annual_damage = mean(annual_damage),
    p_hur_year         = mean(n_hur_days > 0),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

activity_means
```

    #> # A tibble: 3 × 6
    #>   location  mean_peak_wind_kt mean_ts_days mean_hur_days mean_annual_damage
    #>   <chr>                 <dbl>        <dbl>         <dbl>              <dbl>
    #> 1 Saba                   55.2         2.62         0.382              0.003
    #> 2 St_Martin              68.0         3.19         0.961              0.007
    #> 3 Statia                 53.8         2.67         0.38               0.002
    #> # ℹ 1 more variable: p_hur_year <dbl>

`p_hur_year` is the model-implied probability of at least one hurricane
day in a season — a useful quick check that the simulation reproduces
the expected basin behaviour.

## 7. Key Parameters — What to Tune and Why

- **`simulation_years`.** Controls Monte Carlo noise. Use ≥ 1,000;
  2,000+ for return-period work above ~ 50 yr.
- **`search_radius_km`.** Track-point gate radius per site. 800 km is
  conservative and usually captures all relevant storms. Reducing it can
  starve the rate fit; increasing it pulls in tracks that contribute
  negligible site wind.
- **`historical_start_year`.** Start of the calibration window. 1970 is
  the practical minimum for reliable Atlantic best tracks; later starts
  shorten the record and inflate sampling uncertainty in `lambda` and
  `k_hat`.
- **`gust_factor`.** Sustained-to-gust conversion. Use 1.0 when
  comparing model output to IBTrACS sustained wind; ~1.3 when the
  downstream model wants gust forcing.
- **`damage`.** `list(method = "intensity")` returns a normalised
  wind-intensity damage proxy; `"powerlaw"` exposes a configurable
  wind-cubed-style damage curve.
- **`pulse_shape`.** `"cosine"` gives a smooth bell envelope;
  `"uniform"` is a flat-top alternative useful for sensitivity testing.
- **`background_wind`.** Pass a
  [`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
  object to fill non-storm days with realistic trade-wind speeds. Keep
  `NULL` for return-level work (zeros represent TC-free days
  unambiguously); use a non-`NULL` value for energy yield, chronic
  exposure, or port disruption analyses that need a continuous wind
  record. See [Section 5.3](#sec-background-wind) for the full workflow.

