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

ibtracs_path <- "/inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"

cfg <- make_hazard_cfg(
  data_path = ibtracs_path,
  search_radius_km = 800,
  historical_start_year = 1970,
  simulation_years = 500
)

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

out_baseline <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  climate = NULL,
  seed = seed,
  verbose = TRUE
)

# =============================================================================
# 3) Climate scenario
# =============================================================================

future_period <- 2035:2065

scenario_table <- tibble::tibble(
  scenario = "ssp585",
  delta_sst = get_scenario_delta("ssp585", future_period = future_period)
)

clim_in <- make_climate_cfg(
  scenario = scenario_table$scenario[[1]],
  start_year = min(future_period),
  sensitivity_mode = c("linear_shifted")
)

out_585 <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  climate = clim_in,
  seed = seed,
  verbose = FALSE)

# =============================================================================
# 5) Compare simulated annual activity across scenarios
# =============================================================================

sim_compare <- bind_rows(
  out_baseline$sim |> mutate(scenario = "Baseline"),
  out_585$sim |> mutate(scenario = "SSP5-8.5"))

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
# 7) Percentage change relative to baseline
# =============================================================================

change_585 <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(scenario %in% c("Baseline", "SSP5-8.5"))

baseline_change <- change_585 |>
  filter(scenario == "Baseline") |>
  transmute(
    location,
    mean_total_baseline = mean_total,
    mean_hur_baseline = mean_hur
  ) |>
  arrange(location)

ssp585_change <- change_585 |>
  filter(scenario == "SSP5-8.5") |>
  transmute(
    location,
    mean_total_585 = mean_total,
    mean_hur_585 = mean_hur
  ) |>
  arrange(location)

change_585 <- baseline_change |>
  left_join(ssp585_change, by = "location") |>
  mutate(
    pct_total_585 = 100 * (mean_total_585 - mean_total_baseline) / mean_total_baseline,
    pct_hur_585 = 100 * (mean_hur_585 - mean_hur_baseline) / mean_hur_baseline
  ) |>
  select(location, pct_total_585, pct_hur_585)

cat("\n--- SSP5-8.5 change relative to baseline (%) ---\n")
print(
  change_585 |>
    mutate(across(where(is.numeric), ~ round(.x, 1)))
)

# =============================================================================
# 8) Daily hazard illustration for one island
# =============================================================================

daily_baseline <- generate_daily_hazard_impact(
  out = out_baseline,
  location = "Saba",
  sim_years = 1:200,
  year0 = min(future_period),
  gust_factor = 1.25,
  damage_method = "powerlaw",
  scenario = "baseline",
  seed = seed
)$Saba

daily_585 <- generate_daily_hazard_impact(
  out = out_585,
  location = "Saba",
  sim_years = 1:200,
  year0 = min(future_period),
  gust_factor = 1.25,
  damage_method = "powerlaw",
  scenario = "ssp585",
  seed = seed
)$Saba

daily_compare <- bind_rows(daily_baseline, daily_585)

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
