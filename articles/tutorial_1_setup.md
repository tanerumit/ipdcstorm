# Tutorial 1 - Introduction, Setup & Baseline Workflow

## Tutorial 1 — Introduction, Setup & Baseline Workflow

### 1. Overview

`ipdcstorm` simulates hurricane winds at specific coastal sites. It uses
NOAA’s IBTrACS historical track database as the backbone, then rolls
thousands of plausible seasons to produce long daily records of the wind
each site would experience. The package was built for climate
stress-testing of supply-chain infrastructure on the Dutch Caribbean
islands (St. Martin, Saba, Statia), but works for any point location in
the North Atlantic basin.

The workflow has three stages:

1.  **Learn from history.** For every site, the model finds every
    historical storm that passed within range, estimates how strong the
    wind was right at the site (using Holland’s 1980 wind-profile
    formula with directional R34 / R50 / R64 wind radii and a
    forward-motion correction), and labels each visit as a tropical
    storm (≥ 34 kt) or hurricane (≥ 64 kt). From the record it learns
    how often each class hits and how much seasons cluster — captured by
    the overdispersion parameter `k_hat`.
2.  **Simulate plausible seasons.** Draw many synthetic years. For each
    year, first roll an activity multiplier `A` (some years are quiet,
    some busy), then draw the total storm count and the storm /
    hurricane split. When climate forcing is switched on (Tutorial 3),
    warmer sea-surface temperatures nudge both the count and the
    hurricane share, and can optionally shift each event’s intensity,
    size, and forward speed.
3.  **Turn counts into daily winds.** For each simulated season, the
    model picks matching real storms from its catalogue, drops them on
    the calendar at realistic dates, and spreads each storm’s peak wind
    over the days it lasts. Overlapping storms are merged by taking the
    highest wind per day. A damage estimate, storm surge, and the
    dominant event’s details are added alongside.

### 2. Prerequisites

#### 2.1. Installation

``` r
# install.packages("devtools")
# devtools::install_github("tanerumit/ipdcstorm")
```

``` r
library(ipdcstorm)
library(dplyr)
```

#### 2.2. IBTrACS data

The model learns from NOAA’s International Best Track Archive for
Climate Stewardship (IBTrACS) — specifically the **North Atlantic
subset** (`ibtracs.NA.list.v04r01.csv`), which records the position and
strength of every Atlantic storm every six hours.

| IBTrACS field(s) | Used for |
|----|----|
| `SID`, `SEASON`, `BASIN` | Storm identification, year, basin filter |
| `USA_LAT`, `USA_LON`, `ISO_TIME` | 6-hourly track position |
| `USA_WIND` | Maximum sustained wind (kt) |
| `USA_R34_NE/SE/SW/NW`, `USA_R50_*`, `USA_R64_*` | Directional wind radii (nm) at 34, 50, 64 kt thresholds |
| `USA_PRES`, `USA_POCI` | Central pressure and outermost closed-isobar pressure |
| `USA_RMW`, `USA_ROCI` | Radius of maximum wind and outermost closed isobar |

Wind radii tell the model how far the strong winds reach from the
storm’s centre. When they are missing (mostly in pre-2004 records), the
model falls back to climatological estimates from Knaff et al. (2015).

**Download:** [IBTrACS
NCEI](https://www.ncei.noaa.gov/products/international-best-track-archive)

### 3. Configuration

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
> | `search_radius_km` | How close a storm must pass to count. Keep 800 km — even storms 700 km away can push tropical-storm winds to a site through their outer wind field. |
> | `simulation_years` | How many synthetic years to generate. |
> | `climate` | A climate scenario. `make_climate_cfg(scenario = "stationary")` means “no warming” — present-day conditions. |
> | `preset` / `advanced` | Expert knobs for class thresholds and rate-scaling. Defaults follow WMO conventions (TS ≥ 34 kt, HUR ≥ 64 kt). |
> | `background_wind` | Optional `background_wind_cfg` object (see [Section 5.3](#sec-background-wind)). Fills non-storm days with trade-wind speeds sampled from site-specific Weibull distributions. `NULL` (default) sets all non-storm days to 0 kt. |
>
> See
> [`?make_hazard_cfg`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md)
> for the full parameter list.

### 4. Model Calibration and Simulation

[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
runs Stages 1 and 2 end to end in a single call: it learns from the
historical record, then simulates `N_SIM` synthetic seasons. With a
stationary `climate` config, no warming is applied.

``` r
out_base <- run_hazard_model(
  cfg     = hazard_cfg,
  targets = targets,
  seed    = SEED,
  verbose = FALSE
)
```

#### 4.1. What happens inside `run_hazard_model()`

**Step 1 — Read and tidy the IBTrACS file**
([`read_ibtracs_clean()`](https://tanerumit.github.io/ipdcstorm/reference/read_ibtracs_clean.md)).
Standardise column names, keep the wind and pressure fields, convert
wind-radii units to km, and restrict to the Atlantic basin.

**Step 2 — Keep only nearby track points.** For every six-hour position
in the record, measure the great-circle distance to each site and
discard points beyond `search_radius_km`.

**Step 3 — Estimate the wind at each site**
([`compute_site_winds_full()`](https://tanerumit.github.io/ipdcstorm/reference/compute_site_winds_full.md)).
Work out which side of the storm the site sits on, look up the
directional wind radii, estimate the eye radius
([`estimate_RMW_knaff()`](https://tanerumit.github.io/ipdcstorm/reference/estimate_RMW_knaff.md)),
and evaluate the Holland wind-profile formula — with a small adjustment
for how fast the storm is moving. The output is the estimated wind
`site_wind_kt` at every track point.

**Step 4 — Group track points into events**
([`make_storm_events()`](https://tanerumit.github.io/ipdcstorm/reference/make_storm_events.md)).
Combine the 6-hourly points into one summary row per (storm, site), with
the peak site wind, mean eye radius, timing, and pressure. Each event is
then labelled TS (34–63 kt) or HUR (≥ 64 kt).

**Step 5 — Measure how often storms hit**
([`compute_annual_counts()`](https://tanerumit.github.io/ipdcstorm/reference/compute_annual_counts.md)
→
[`compute_lambda_table()`](https://tanerumit.github.io/ipdcstorm/reference/compute_lambda_table.md)).
Turn the per-year counts into average rates (`lambda`) and measure how
clustered busy and quiet seasons are (`k_hat`).

**Step 6 — Roll synthetic seasons**
([`simulate_twolevel_counts()`](https://tanerumit.github.io/ipdcstorm/reference/simulate_twolevel_counts.md)).
For each synthetic year: first draw an activity multiplier `A` (is this
a quiet or hyperactive season?), then a total storm count consistent
with `A`, and finally split that total into tropical storms and
hurricanes.

> **Always check the rates first**
>
> The first thing to look at after a run is `out_base$rates`. If a
> site’s `lambda` looks very different from the rate you expect (say,
> from a published climatology), every downstream result will be biased
> too — fix that before running anything else.

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

#### 4.2. Inspecting `out_base`

The return value is a named list with one tidy table per stage.

##### 4.2.1. Rates

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

| Column | Description |
|----|----|
| `location` | Target site. |
| `storm_class` | `TS` (34–63 kt site wind) or `HUR` (≥ 64 kt). |
| `lambda` | Average number of this kind of event per year. |
| `n_years` | How many historical years the estimate is based on. |
| `prob_annual` | Chance of at least one event in a year, \$\`1 - e^{-\lambda}\`\$. |
| `prob_none` | Chance of no events, \$\`e^{-\lambda}\`\$. |

A `lambda` of 1.2 for TS means a site averages about 1.2 tropical storms
per year — but some years are quiet and others see several, because
seasons naturally cluster.

##### 4.2.2. Fit

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
| `k_hat` | Season-to-season clustering. Small values (~ 2–5) mean very bursty seasons (some quiet, some hyperactive); large values mean seasons are fairly uniform. |
| `mean_annual_total` / `var_annual_total` | Mean and spread of the yearly storm count. |
| `beta_sst` | How strongly warming changes storm rates (zero for present-day). |
| `gamma_intensity` | How strongly warming shifts the TS / hurricane mix (zero for present-day). |
| `p_hurricane_base` | Baseline chance that a given storm reaches hurricane strength. |

##### 4.2.3. Sim

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
| `sim_year` | Which synthetic year this row describes (1 to `N_SIM`). |
| `n_total` / `n_tc` | Total storms of any class that hit the site this year. |
| `n_ts` | Tropical-storm-class visits. |
| `n_hur` | Hurricane-class visits. |
| `A` | Activity roll for the year: \> 1 means busier than normal, \< 1 means quieter. |
| `p_hur` | Share of storms that escalate to hurricane strength. |

This table feeds the daily downscaling and every return-period
calculation that follows.

##### 4.2.4. Events

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

One row per (historical storm × site), showing the wind the site
actually felt:

| Column | Description |
|----|----|
| `storm_id` | IBTrACS ID (e.g. `2017242N16333` for Irma). |
| `year` | Season year. |
| `peak_wind_kt` | **Peak wind at the site** — after accounting for how far away the storm passed. |
| `storm_intensity_kt` | The storm’s strongest wind anywhere in its life — often much higher. |
| `min_pressure_hpa` / `pressure_deficit_hpa` | Lowest central pressure and its gap below the surrounding environment. |
| `rmw_mean_km` | Average radius of the eyewall. |
| `storm_class` | Classification at the site (`TS` or `HUR`). |
| `location` | Target site. |

The gap between `peak_wind_kt` and `storm_intensity_kt` matters: a Cat-5
storm (≈ 155 kt) passing 400 km away might only produce Cat-1 winds (~
70 kt) at the site. The Holland wind profile captures that drop-off.

##### 4.2.5. Trackpoints

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

This is the rawest output: every 6-hour IBTrACS track point near each
site, plus all the computed wind-field values (`site_wind_kt`,
`dist_km`, `bearing_deg`, `RMW_km`, and more). Handy when you want to
spot-check a specific storm.

### 5. Daily Wind Time Series

[`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)
turns the yearly storm counts into a day-by-day wind record for every
simulated year and every site.

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

#### 5.1. What the function does

**Step 1 — Build a catalogue of historical events.** For each site,
collect every past event into a library that keeps peak wind, duration,
pressure and eye radius bundled together. A shared pool across sites
keeps multi-site sampling consistent.

**Step 2 — Pick events for each simulated year.** For every `sim_year`,
draw the required number of tropical storms and hurricanes from the
catalogue. Each storm keeps its real physical profile (wind, duration,
pressure, eye radius all come from the same event).

**Step 3 — Drop each event onto the calendar.** Each event gets a
realistic start date drawn from the historical day-of-year pattern,
which puts most activity in August–October.

**Step 4 — Spread each event’s wind over its duration.** The `cosine`
pulse shape builds a smooth bell curve: winds rise, peak, and fall over
the event’s life. Use `pulse_shape = "uniform"` for a flat-top
alternative.

**Step 5 — Handle overlaps.** When two storms hit the same day, take the
higher wind. No double-counting, but compound seasons are still visible.

**Step 6 — Add derived daily fields.** Compute a damage estimate
(`damage = list(method = "intensity")` gives a normalised wind-based
proxy; `"powerlaw"` is the alternative), estimate `surge_m`, and carry
through the dominant event’s details (`event_id`, `pressure_hpa`,
`rmw_km`).

> **Compare two runs side by side**
>
> Use the same `seed` and the same `sim_years` in both a baseline run
> and a climate-change run. Because the random draws are identical, any
> difference in the daily output comes from the climate signal alone.

> **Force a specific storm into every year**
>
> The `pinned_sids` option lets you lock a historical storm (e.g. Irma)
> into every synthetic year. Used by Tutorial 4’s stress-test workflow.

#### 5.2. Daily output schema

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

#### 5.3. Background Wind

By default, any day without a storm has `wind_kt = 0`. That is clean for
return-period analysis (zeros clearly mean “no storm”), but it is
unrealistic for applications that need a continuous wind record — for
example, wind-energy assessments or chronic-exposure studies.

[`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
fills the quiet days with realistic trade-wind speeds. You supply a
monthly wind distribution per site; the function also supports spatial
correlation across sites and a small day-to-day persistence so that one
windy day tends to follow another.

> **Storm-day winds stay the same**
>
> On any day where the storm wind is stronger than the background trade
> wind, the storm value wins. Return-period statistics at TS and
> hurricane thresholds are therefore unchanged whether you use
> background wind or not.

### 6. Annual Activity Summary

Roll the daily series up into one row per site per year.

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

`p_hur_year` is the model’s estimate of the chance a site sees at least
one hurricane day in a given year — a quick sanity check that the
simulation matches what you’d expect from the region.

### 7. Key Parameters — What to Tune and Why

- **`simulation_years`.** How many synthetic seasons to roll. More
  seasons mean smoother statistics. Use at least 1,000; 2,000+ for
  return periods beyond ~ 50 years.
- **`search_radius_km`.** How wide to cast the net for historical
  storms. 800 km is a safe default that captures all relevant hits.
  Smaller values risk missing real cases; larger values add tracks that
  barely touch the site.
- **`historical_start_year`.** First year of the record used to learn
  rates. 1970 is the safe minimum — satellite coverage before that is
  patchy. Later starts shorten the record and make `lambda` and `k_hat`
  less reliable.
- **`gust_factor`.** Converts sustained wind to gust. Keep at 1.0 when
  comparing with IBTrACS; bump to ~ 1.3 when the downstream analysis
  needs peak gusts (typical for damage and infrastructure studies).
- **`damage`.** `list(method = "intensity")` gives a simple wind-based
  damage score; `"powerlaw"` uses a tuneable wind-cubed damage curve.
- **`pulse_shape`.** `"cosine"` draws a smooth bell curve of wind over
  the event; `"uniform"` gives every storm day the same wind — useful
  for sensitivity tests.
- **`background_wind`.** Pass a
  [`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
  object to fill quiet days with realistic trade winds.

### 8. References

Holland, G. J. (1980). An analytic model of the wind and pressure
profiles in hurricanes. *Monthly Weather Review*, *108*(8), 1212–1218.
https://doi.org/10.1175/1520-0493(1980)108%3C1212:AAMOTW%3E2.0.CO;2

Knaff, J. A., Longmore, S. P., DeMaria, R. T., & Molenar, D. A. (2015).
Improved tropical-cyclone flight-level wind estimates using routine
infrared satellite reconnaissance. *Journal of Applied Meteorology and
Climatology*, *54*(2), 463–478. https://doi.org/10.1175/JAMC-D-14-0112.1

Knaff, J. A., & Zehr, R. M. (2007). Reexamination of tropical cyclone
wind–pressure relationships. *Weather and Forecasting*, *22*(1), 71–88.
https://doi.org/10.1175/WAF965.1

Knapp, K. R., Diamond, H. J., Kossin, J. P., Kruk, M. C., & Schreck, C.
J. (2018). *International Best Track Archive for Climate Stewardship
(IBTrACS) project, version 4* \[Data set\]. NOAA National Centers for
Environmental Information. https://doi.org/10.25921/82ty-9e16

Knapp, K. R., Kruk, M. C., Levinson, D. H., Diamond, H. J., & Neumann,
C. J. (2010). The International Best Track Archive for Climate
Stewardship (IBTrACS): Unifying tropical cyclone data. *Bulletin of the
American Meteorological Society*, *91*(3), 363–376.
https://doi.org/10.1175/2009BAMS2755.1

World Meteorological Organization. (2017). *Global guide to tropical
cyclone forecasting* (WMO-No. 1194). World Meteorological Organization.
https://library.wmo.int/idurl/4/56236
