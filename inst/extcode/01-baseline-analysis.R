# =============================================================================
# Baseline hazard analysis pipeline
#
# Purpose:
#   Set up and run a stationary baseline tropical-cyclone hazard model for
#   a set of Caribbean target locations, then produce daily hazard-impact
#   time series for every simulated year.
#
#   The baseline is a stationary Poisson/NB model calibrated to the
#   IBTrACS historical record (default: 1970–present).  No climate forcing
#   is applied; SST is treated as fixed at its historical mean.  Use
#   02-model-validation.R to check model skill before drawing conclusions,
#   and 03-climate-scenarios.R to layer KNMI'23 warming scenarios on top.
#
# Key config switches (all in the Config block below):
#   SEED             — global RNG seed for full reproducibility
#   N_SIM            — number of synthetic years (≥ 1 000 recommended)
#   HIST_START_YEAR  — first IBTrACS year used to fit storm rates
#   SEARCH_RADIUS_KM — storm-track proximity radius per site (km)
#   YEAR0            — calendar year assigned to sim_year = 1 in daily output
#   GUST_FACTOR      — multiplicative gust-to-sustained-wind conversion
#
# Outputs:
#   hazard_cfg      — hazard_cfg object (model parameters)
#   out_base        — run_hazard_model() output: events, fit, rates, sim
#   daily_by_loc    — named list[location] of daily tibbles (all N_SIM years)
#   annual_summary  — tibble: per-location annual storm activity statistics
#
# Data path:
#   Uses the packaged demo IBTrACS subset by default.
#   For production, replace `DATA_PATH` with the full IBTrACS NA CSV path.
# =============================================================================

library(dplyr)
library(ipdcstorm)

# =============================================================================
# Config
# =============================================================================

SEED             <- 42L
N_SIM            <- 2000L
HIST_START_YEAR  <- 1970L
SEARCH_RADIUS_KM <- 800

# Calendar year mapped to the first synthetic year in daily output.
# Setting YEAR0 = 2000 means sim_year = 1 → calendar year 2000,
# sim_year = 2 → 2001, and so on.  Adjust to align with your analysis window.
YEAR0 <- 2000L

# Multiplicative conversion from 1-minute sustained wind to 3-second gust.
# 1.0 → no conversion (sustained wind, useful for model validation vs IBTrACS).
# 1.3 → standard WMO gust factor (use for damage or infrastructure impact).
GUST_FACTOR <- 1.0

# Target locations (Dutch Caribbean islands)
targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740
)

# IBTrACS dataset CSV for the analysis
DATA_PATH <- system.file("extdata", "ibtracs_1970.csv", package = "ipdcstorm")
#DATA_PATH <- system.file("extdata", "ibtracs_demo.csv", package = "ipdcstorm")


# =============================================================================
# Part 1: Configure and run the stationary baseline model
# =============================================================================

message("Part 1: Configuring stationary baseline model")
message("  Seed            : ", SEED)
message("  Simulation years: ", N_SIM)
message("  Historical start: ", HIST_START_YEAR)
message("  Search radius   : ", SEARCH_RADIUS_KM, " km")
message("  Locations       : ", paste(targets$name, collapse = ", "))

hazard_cfg <- make_hazard_cfg(
  data_path             = DATA_PATH,
  search_radius_km      = SEARCH_RADIUS_KM,
  historical_start_year = HIST_START_YEAR,
  simulation_years      = N_SIM,
  climate               = make_climate_cfg(scenario = "stationary")
)

message("\n  Running hazard model ...")

out_base <- run_hazard_model(
  cfg     = hazard_cfg,
  targets = targets,
  seed    = SEED,
  verbose = TRUE
)

# --- Inspect calibrated storm rates ---
#
# out_base$rates contains per-location, per-class (TS / HUR) annual rates
# (lambda) fitted from the historical record over the selected search radius.
# These are the Poisson rates driving synthetic storm counts in simulation.

message("\nCalibrated annual storm rates:")
print(
  out_base$rates |>
    select(location, storm_class, lambda) |>
    mutate(lambda = round(lambda, 3)) |>
    arrange(location, storm_class)
)




# =============================================================================
# Part 2: Generate daily hazard-impact series
#
# Converts the discrete event catalogue (out_base$sim) into a continuous
# daily wind time series for each synthetic year and location.
#
# Each day carries:
#   wind_kt      — site-level 1-min sustained wind (kt)
#   event_class  — storm classification at that day (TS / HUR / NA)
#   damage_rate  — fractional damage proxy [0, 1] derived from wind intensity
#
# pulse_shape = "cosine" applies a bell-shaped intensity envelope over the
# passage window of each event, which is more physically realistic than a
# uniform box profile.
# =============================================================================

message("\nPart 2: Generating daily hazard-impact series (", N_SIM, " years per location)")

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

# =============================================================================
# Part 3: Annual activity summary
#
# Collapse the daily series into per-location annual statistics.
# `peak_wind_kt`  — highest single-day wind in each simulated year (kt)
# `n_ts_days`     — annual count of days with wind >= 34 kt (tropical storm)
# `n_hur_days`    — annual count of days with wind >= 64 kt (hurricane)
# `annual_damage` — integrated damage rate summed over all days of the year
# =============================================================================

message("\nPart 3: Computing annual activity summary ...")

daily_all <- dplyr::bind_rows(daily_by_loc)

annual_summary <- daily_all |>
  group_by(location, sim_year) |>
  summarise(
    peak_wind_kt  = max(wind_kt,    na.rm = TRUE),
    n_ts_days     = sum(wind_kt >= 34, na.rm = TRUE),
    n_hur_days    = sum(wind_kt >= 64, na.rm = TRUE),
    annual_damage = sum(damage_rate, na.rm = TRUE),
    .groups = "drop"
  )

# Per-location means across all simulated years
activity_means <- annual_summary |>
  group_by(location) |>
  summarise(
    mean_peak_wind_kt  = mean(peak_wind_kt),
    mean_ts_days       = mean(n_ts_days),
    mean_hur_days      = mean(n_hur_days),
    mean_annual_damage = mean(annual_damage),
    p_hur_year         = mean(n_hur_days > 0),
    .groups = "drop"
  )

message("\nMean annual hazard statistics across ", N_SIM, " simulated years:")
print(
  activity_means |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# =============================================================================
# Summary
# =============================================================================

message("\n--- Baseline analysis complete ---")
message("Seed              : ", SEED)
message("Simulation years  : ", N_SIM)
message("Historical start  : ", HIST_START_YEAR)
message("Search radius     : ", SEARCH_RADIUS_KM, " km")
message("Gust factor       : ", GUST_FACTOR)
message("Calendar year0    : ", YEAR0)
message("Locations         : ", paste(targets$name, collapse = ", "))
message("")
message("Objects in workspace:")
message("  hazard_cfg      — hazard_cfg: model configuration (pass to validation or scenario scripts)")
message("  out_base        — list: run_hazard_model() output (events, fit, rates, sim, config)")
message("  daily_by_loc    — list[location]: daily tibbles for all ", N_SIM, " simulated years")
message("  annual_summary  — tibble: annual storm statistics per location and sim_year")
message("  activity_means  — tibble: mean annual statistics aggregated across all simulated years")
