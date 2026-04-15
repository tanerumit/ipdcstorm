# =============================================================================
# Climate scenario analysis pipeline
#
# Purpose:
#   Compare a stationary baseline against KNMI'23 climate scenarios to
#   quantify how projected SST warming affects tropical-cyclone hazard at
#   Dutch Caribbean target locations.
#
#   Part 1 — Run the stationary baseline (same configuration as
#             01-baseline-analysis.R) to establish the reference climatology.
#
#   Part 2 — Run each KNMI'23 scenario under the same RNG seed using the
#             common random numbers (CRN) approach.  Differences between
#             scenario and baseline outputs therefore isolate the climate
#             forcing signal rather than sampling noise.
#
#   Part 3 — Compare annual storm activity (counts, hurricane fraction) and
#             daily hazard metrics (TC days, hurricane days, damage) across
#             the baseline and all scenarios.
#
# Key config switches (all in the Config block below):
#   SEED           — RNG seed; must be identical across baseline and scenarios
#   N_SIM          — synthetic years (≥ 1 000 recommended; 500 for quick runs)
#   TARGET_YEAR    — 2050 (near-term) or 2100 (end-of-century) for delta_sst
#   SCENARIOS      — named subset of KNMI'23 scenario IDs to run
#   PERTURB_EVENTS — TRUE applies perturb_event() intensity shift per storm
#
# Outputs:
#   out_base       — run_hazard_model() output for the stationary baseline
#   out_cc         — named list[scenario]: run_hazard_model() output per scenario
#   daily_base     — named list[location]: daily tibbles for baseline
#   daily_cc       — named list[scenario][location]: daily tibbles per scenario
#   sim_compare    — combined sim table with scenario column for plotting
#   daily_compare  — combined daily table with scenario column for plotting
#
# Data path:
#   Uses the packaged demo IBTrACS subset by default.
#   For production, replace `DATA_PATH` with the full IBTrACS NA CSV path.
# =============================================================================

library(dplyr)
library(ggplot2)
library(ipdcstorm)

# =============================================================================
# Config
# =============================================================================

SEED <- 42L
N_SIM <- 2000L
HIST_START_YEAR  <- 1970L
SEARCH_RADIUS_KM <- 800

# Future target year for KNMI'23 delta_sst lookup.
# 2050 → near-term forcing (moderate SST increase for all scenarios)
# 2100 → end-of-century forcing (knmi_Hd reaches ~2.6 °C above baseline)
TARGET_YEAR <- 2050

# KNMI'23 scenario IDs to run.
# Available: "knmi_Ld" (low-end dry), "knmi_Ln" (low-end neutral),
#            "knmi_Hd" (high-end dry), "knmi_Hn" (high-end neutral).
# Running knmi_Ld + knmi_Hd captures the low-to-high plausible spread.
SCENARIOS <- c("knmi_Ld", "knmi_Hd")

# Apply per-storm intensity perturbation (perturb_event()).
# TRUE  — each sampled event's intensity and duration are shifted upward by
#         the scenario's delta_sst, reflecting the projected intensification
#         of individual storms under warming (thermodynamic pathway).
# FALSE — only storm frequency (Poisson rate) is adjusted; no per-storm shift.
PERTURB_EVENTS <- TRUE

# Target locations (Dutch Caribbean islands)
targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740
)

# Data: packaged demo subset — replace with full IBTrACS CSV for production runs
DATA_PATH <- system.file("extdata", "ibtracs_demo.csv", package = "ipdcstorm")
# DATA_PATH <- "/path/to/ibtracs.NA.list.v04r01.csv"

# =============================================================================
# Part 1: Stationary baseline
#
# Run the baseline under the same seed that will be used for scenarios.
# The CRN approach ensures that seed-driven sampling differences between
# baseline and scenario outputs reflect climate forcing only, not random
# variation in storm draws.
# =============================================================================

message("Part 1: Stationary baseline (", N_SIM, " years, seed = ", SEED, ")")

cfg_base <- make_hazard_cfg(
  data_path             = DATA_PATH,
  search_radius_km      = SEARCH_RADIUS_KM,
  historical_start_year = HIST_START_YEAR,
  simulation_years      = N_SIM,
  climate               = make_climate_cfg(scenario = "stationary")
)

out_base <- run_hazard_model(cfg_base, targets = targets, seed = SEED, verbose = FALSE)

message("  Generating baseline daily series ...")

daily_base <- generate_daily_hazard_impact_spatial(
  out         = out_base,
  location    = targets$name,
  sim_years   = seq_len(N_SIM),
  year0       = 2000L,
  gust_factor = 1.0,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "stationary",
  seed        = SEED
)

message("  Done.")

# =============================================================================
# Part 2: KNMI'23 climate scenarios
#
# Each scenario is run with the identical seed, so random draws follow the
# same sequence.  The make_climate_cfg() call resolves the KNMI'23 delta_sst
# for the requested TARGET_YEAR and passes it to run_hazard_model(), where
# it adjusts both the Poisson/NB storm rate (beta_sst pathway) and,
# optionally, individual event properties (perturb_event() pathway).
# =============================================================================

message("\nPart 2: KNMI'23 climate scenarios at target year ", TARGET_YEAR)
message("  Scenarios       : ", paste(SCENARIOS, collapse = ", "))
message("  Perturb events  : ", PERTURB_EVENTS)

out_cc    <- vector("list", length(SCENARIOS))
daily_cc  <- vector("list", length(SCENARIOS))
names(out_cc)   <- SCENARIOS
names(daily_cc) <- SCENARIOS

for (scen in SCENARIOS) {

  message("  [", scen, "] Running model ...")

  cfg_cc <- make_hazard_cfg(
    data_path             = DATA_PATH,
    search_radius_km      = SEARCH_RADIUS_KM,
    historical_start_year = HIST_START_YEAR,
    simulation_years      = N_SIM,
    climate               = make_climate_cfg(
      scenario    = scen,
      target_year = TARGET_YEAR,
      perturb     = PERTURB_EVENTS
    )
  )

  # Same seed → CRN; all between-scenario differences trace to climate forcing.
  out_cc[[scen]] <- run_hazard_model(cfg_cc, targets = targets, seed = SEED, verbose = FALSE)

  message("  [", scen, "] Generating daily series ...")

  daily_cc[[scen]] <- generate_daily_hazard_impact_spatial(
    out         = out_cc[[scen]],
    location    = targets$name,
    sim_years   = seq_len(N_SIM),
    year0       = 2000L,
    gust_factor = 1.0,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = scen,
    seed        = SEED
  )

  message("  [", scen, "] Done.")
}

# =============================================================================
# Part 3: Compare annual storm activity across scenarios
#
# Combine the sim tables from all runs into a single tibble.  Each row is one
# (location, sim_year, scenario) combination; `n_total`, `n_ts`, and `n_hur`
# are the synthetic event counts for that year drawn from the calibrated model.
#
# `p_hurricane` is the parameterised probability that a storm reaching the
# site escalates to hurricane intensity under the given SST scenario; it is
# derived from the gamma parameter in run_hazard_model().
# =============================================================================

message("\nPart 3: Comparing annual storm activity across scenarios ...")

sim_compare <- bind_rows(
  out_base$sim |> mutate(scenario = "stationary"),
  lapply(SCENARIOS, function(s) out_cc[[s]]$sim |> mutate(scenario = s)) |>
    bind_rows()
)

activity_summary <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    mean_n_total   = mean(n_total, na.rm = TRUE),
    mean_n_hur     = mean(n_hur,   na.rm = TRUE),
    mean_p_hur     = mean(p_hurricane, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

cat("\nMean annual storm activity by scenario and location:\n")
print(
  activity_summary |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# Delta relative to stationary baseline (percentage change)
activity_delta <- activity_summary |>
  left_join(
    activity_summary |>
      filter(scenario == "stationary") |>
      select(location, base_n_total = mean_n_total, base_n_hur = mean_n_hur),
    by = "location"
  ) |>
  filter(scenario != "stationary") |>
  mutate(
    delta_n_total_pct = 100 * (mean_n_total - base_n_total) / base_n_total,
    delta_n_hur_pct   = 100 * (mean_n_hur   - base_n_hur)   / pmax(base_n_hur, 1e-9)
  ) |>
  select(scenario, location, delta_n_total_pct, delta_n_hur_pct)

cat("\nClimate change signal (% change vs stationary baseline):\n")
print(
  activity_delta |>
    mutate(across(where(is.numeric), ~ round(.x, 1)))
)

# =============================================================================
# Part 4: Compare daily hazard metrics across scenarios
#
# Aggregate the daily series to annual totals and compare distributions.
# Annual damage is the sum of damage_rate across all days of the year;
# it integrates both the intensity (per-day damage) and duration (active days).
# =============================================================================

message("\nPart 4: Comparing daily hazard metrics across scenarios ...")

# Combine all daily series into one tibble with a scenario label
daily_compare <- bind_rows(
  dplyr::bind_rows(daily_base) |> mutate(scenario = "stationary"),
  lapply(SCENARIOS, function(s) {
    dplyr::bind_rows(daily_cc[[s]]) |> mutate(scenario = s)
  }) |>
    bind_rows()
)

daily_annual <- daily_compare |>
  group_by(scenario, location, sim_year) |>
  summarise(
    n_tc_days    = sum(wind_kt > 0,    na.rm = TRUE),
    n_ts_days    = sum(wind_kt >= 34,  na.rm = TRUE),
    n_hur_days   = sum(wind_kt >= 64,  na.rm = TRUE),
    annual_damage = sum(damage_rate,   na.rm = TRUE),
    .groups = "drop"
  )

daily_summary <- daily_annual |>
  group_by(scenario, location) |>
  summarise(
    mean_tc_days     = mean(n_tc_days),
    mean_hur_days    = mean(n_hur_days),
    mean_ann_damage  = mean(annual_damage),
    p90_ann_damage   = quantile(annual_damage, 0.90),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

cat("\nMean annual daily hazard metrics by scenario:\n")
print(
  daily_summary |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# =============================================================================
# Part 5: Quick comparison plot — mean annual storm count by scenario
# =============================================================================

p_activity <- activity_summary |>
  mutate(
    scenario = factor(scenario, levels = c("stationary", SCENARIOS))
  ) |>
  ggplot(aes(x = location, y = mean_n_total, fill = scenario)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c(
      stationary = "#7f8c8d",
      knmi_Ld    = "#3498db",
      knmi_Ln    = "#2980b9",
      knmi_Hd    = "#e74c3c",
      knmi_Hn    = "#c0392b"
    ),
    name = "Scenario"
  ) +
  labs(
    x     = NULL,
    y     = "Mean annual storm count",
    title = sprintf("Annual storm activity — baseline vs KNMI'23 (%d)", TARGET_YEAR)
  ) +
  theme_light(base_size = 11)

print(p_activity)

# =============================================================================
# Summary
# =============================================================================

message("\n--- Climate scenario analysis complete ---")
message("Seed              : ", SEED)
message("Simulation years  : ", N_SIM)
message("Historical start  : ", HIST_START_YEAR)
message("Search radius     : ", SEARCH_RADIUS_KM, " km")
message("Target year       : ", TARGET_YEAR)
message("Scenarios         : ", paste(SCENARIOS, collapse = ", "))
message("Perturb events    : ", PERTURB_EVENTS)
message("")
message("Objects in workspace:")
message("  out_base       — list: baseline run_hazard_model() output")
message("  out_cc         — list[scenario]: scenario run_hazard_model() outputs")
message("  daily_base     — list[location]: baseline daily tibbles (all ", N_SIM, " years)")
message("  daily_cc       — list[scenario][location]: scenario daily tibbles")
message("  sim_compare    — tibble: combined sim table with scenario column")
message("  daily_compare  — tibble: combined daily table with scenario column")
message("  activity_summary — tibble: mean annual storm counts per scenario/location")
message("  activity_delta   — tibble: % change in storm counts vs stationary baseline")
message("  daily_summary    — tibble: mean annual hazard metrics per scenario/location")
message("  p_activity       — ggplot: bar chart of mean annual storm count by scenario")
