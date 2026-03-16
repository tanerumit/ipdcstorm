# =============================================================================
# Simple climate-change assessment pipeline
#
# Purpose:
#   Illustrate a clean baseline-vs-climate comparison using the current
#   stationary time-slice climate API.
#
# What this example does:
#   1) Define the hazard model setup and target locations
#   2) Run a stationary baseline
#   3) Build two climate scenarios (SSP2-4.5 and SSP5-8.5)
#   4) Compare annual activity and hurricane fraction across scenarios
#   5) Generate a daily series example for one island
#
# Scientific note:
#   This script uses the simplified climate workflow: scenario `delta_sst` is
#   resolved from the built-in scenario table and historical sensitivities are
#   calibrated inside `run_hazard_model()`. It does not build transient or
#   nested climate-response objects inside the pipeline.
# =============================================================================

library(dplyr)
library(ggplot2)
library(ipdcstorm)

# =============================================================================
# 1) Core model setup
# =============================================================================

seed <- 123L
simulation_years <-500
ibtracs_path <- "/inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"


targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)



# =============================================================================
# 2) Stationary baseline
# =============================================================================


hazard_cfg <- make_hazard_cfg(
  data_path = ibtracs_path,
  search_radius_km = 800,
  historical_start_year = 1970,
  simulation_years = simulation_years,
  climate = make_climate_cfg(scenario = "stationary")
)

out_baseline <- run_hazard_model(
  cfg = hazard_cfg,
  targets = targets,
  seed = seed,
  verbose = FALSE
)

# =============================================================================
# 3) Climate scenario
# =============================================================================

# 1) Set hypothetical delta_sst directly for the target year.
delta_sst_value <- 1.5

# 2) Make the climate configuration file
climate_hyp <- make_climate_cfg(
  delta_sst = delta_sst_value,
  sensitivity_mode = "fixed")

print(climate_hyp)

scenario_cfg <- make_hazard_cfg(
  data_path = ibtracs_path,
  search_radius_km = 800,
  historical_start_year = 1970,
  simulation_years = simulation_years,
  climate = climate_hyp)

out_scenario <- run_hazard_model(
  cfg = scenario_cfg,
  targets = targets,
  seed = seed,
  verbose = TRUE)

# =============================================================================
# 5) Compare simulated annual activity across scenarios
# =============================================================================

sim_compare <- bind_rows(
  out_baseline$sim |> mutate(scenario = "Baseline"),
  out_scenario$sim |> mutate(scenario = "future_scenario"))

activity_summary <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_ts = mean(n_ts, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    mean_p_hur = mean(n_hur / pmax(1, n_total), na.rm = TRUE),
    .groups = "drop") |>
  arrange(location, scenario)

cat("\n--- Mean annual activity by scenario ---\n")
print(
  activity_summary |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# =============================================================================
# 6) Compare hurricane fraction across scenarios
# =============================================================================

split_summary <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    p_hur_sim = mean(n_hur / pmax(1, n_total), na.rm = TRUE),
    mean_p_hur_param = mean(p_hurricane, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

cat("\n--- Hurricane fraction by scenario ---\n")
print(
  split_summary |>
    mutate(across(where(is.numeric), ~ round(.x, 4)))
)

# =============================================================================
# 8) Daily hazard illustration for one island
# =============================================================================

daily_baseline <- generate_daily_hazard_impact(
  out = out_baseline,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = "baseline",
  seed = seed
)$Saba

daily_scenario <- generate_daily_hazard_impact(
  out = out_scenario,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = "future",
  seed = seed
)$Saba

daily_compare <- bind_rows(daily_baseline, daily_scenario)

daily_summary <- daily_compare |>
  group_by(scenario) |>
  summarise(
    n_days = n(),
    tc_days = sum(wind_kt > 0, na.rm = TRUE),
    ts_days = sum(wind_kt >= 34, na.rm = TRUE),
    hur_days = sum(wind_kt >= 64, na.rm = TRUE),
    mean_annual_damage = sum(damage_rate, na.rm = TRUE) / 200,
    .groups = "drop"
  )

cat("\n--- Daily impact comparison for Saba (200 simulated years) ---\n")
print(
  daily_summary |>
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

# =============================================================================
# 9) Optional quick plot
# =============================================================================

p_activity <- activity_summary |>
  ggplot(aes(x = location, y = mean_total, fill = scenario)) +
  geom_col(position = "dodge") +
  labs(
    x = NULL,
    y = "Mean annual storm count",
    title = "Simulated annual storm activity by scenario"
  ) +
  theme_light(base_size = 11)

print(p_activity)


