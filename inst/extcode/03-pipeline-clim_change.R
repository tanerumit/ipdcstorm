# =============================================================================
# Simple climate-change assessment pipeline
#
# Purpose:
#   Illustrate a clean baseline-vs-climate comparison using the refactored API.
#
# What this example does:
#   1) Define the hazard model setup and target locations
#   2) Run a stationary baseline
#   3) Run two climate scenarios (SSP2-4.5 and SSP5-8.5)
#   4) Compare annual activity and hurricane fraction across scenarios
#   5) Generate a daily series example for one island

# =============================================================================

rm(list = ls())

# =============================================================================
# Packages
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ipdcstorm)

# =============================================================================
# 1) Core model setup
# =============================================================================
# Define the historical hazard model configuration and the target locations.

cfg <- make_hazard_cfg(
  data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  search_radius_km = 800,
  start_year = 1970,
  n_sim = 1000
)

targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

# Randomization seed (use for deterministic results)
seed <- 123

# =============================================================================
# 2) Stationary baseline
# =============================================================================

# Run the model without climate adjustments. This is the reference case.
out_baseline <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  seed = seed,
  climate_cfg = NULL

)

# =============================================================================
# 3) Climate scenarios
# =============================================================================

# Define climate scenarios via rate effects (beta_sst) and hurricane-fraction/intensity effects (gamma)
# Storm perturbation is optional and off unless explicitly supplied.

clim_245 <- make_climate_cfg(
  enabled = TRUE,
  scenario = "ssp245",             # future climate pathway used to generate projected SST anomalies
  sst_source = "builtin",          # source of historical MDR SST data used to compute anomalies
  baseline_years = 1991L:2020L,    # user-specified
  start_year = 2025L,              # user-specified
  beta_sst = NULL,                 # estimate from data
  gamma = NULL,                    # estimate from data
  perturb = NULL                   # no storm perturbation in this simple illustration
)

scen_245 <- generate_sst_scenario(
  n_years = 200,
  mode = clim_245$scenario,
  start_year = clim_245$start_year,
  baseline_years = clim_245$baseline_years
)

head(scen_245)
plot(scen_245$calendar_year, scen_245$sst_anomaly, type = "l")





clim_585 <- make_climate_cfg(
  enabled = TRUE,
  scenario = "ssp585",
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  start_year = 2025L,
  beta_sst = NULL,
  gamma = NULL,
  perturb = NULL
)

# Optional variant:
# If you want to include storm perturbation in a simple way, use:
#
# clim_585 <- make_climate_cfg(
#   enabled = TRUE,
#   scenario = "ssp585",
#   sst_source = "builtin",
#   baseline_years = 1991L:2020L,
#   start_year = 2025L,
#   beta_sst = NULL,
#   gamma = NULL,
#   perturb = list(
#     v_scale = 0.05,
#     r_scale = 0.08,
#     speed_scale = -0.10,
#     precip_scale = 0.07
#   )
# )

# =============================================================================
# 4) Climate runs
# =============================================================================
# Run the same hazard model under two climate scenarios.

set.seed(seed)
out_245 <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  seed = seed,
  climate_cfg = clim_245
)

set.seed(seed)
out_585 <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  seed = seed,
  climate_cfg = clim_585
)

# =============================================================================
# 5) Compare simulated annual activity across scenarios
# =============================================================================
# Combine model outputs and summarise mean annual storm activity.

sim_compare <- bind_rows(
  out_baseline$sim %>% mutate(scenario = "Baseline"),
  out_245$sim      %>% mutate(scenario = "SSP2-4.5"),
  out_585$sim      %>% mutate(scenario = "SSP5-8.5")
)

activity_summary <- sim_compare %>%
  group_by(scenario, location) %>%
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_ts = mean(n_ts, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    mean_p_hur = mean(n_hur / pmax(1, n_total), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(location, scenario)

cat("\n--- Mean annual activity by scenario ---\n")
print(
  activity_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# =============================================================================
# 6) Compare hurricane fraction across scenarios
# =============================================================================
# This focuses on the modeled hurricane share rather than total activity alone.

split_summary <- sim_compare %>%
  group_by(scenario, location) %>%
  summarise(
    p_hur_sim = mean(n_hur / pmax(1, n_total), na.rm = TRUE),
    mean_p_hur_param = mean(p_hurricane, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(location, scenario)

cat("\n--- Hurricane fraction by scenario ---\n")
print(
  split_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
)

# =============================================================================
# 7) Percentage change relative to baseline
# =============================================================================
# Calculate simple percentage changes for SSP5-8.5 relative to the baseline.

change_table <- sim_compare %>%
  group_by(scenario, location) %>%
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    .groups = "drop"
  )

base_table <- change_table %>%
  filter(scenario == "Baseline") %>%
  transmute(
    location,
    mean_total_baseline = mean_total,
    mean_hur_baseline = mean_hur
  )

ssp585_table <- change_table %>%
  filter(scenario == "SSP5-8.5") %>%
  transmute(
    location,
    mean_total_585 = mean_total,
    mean_hur_585 = mean_hur
  )

change_585 <- base_table %>%
  left_join(ssp585_table, by = "location") %>%
  mutate(
    pct_total_585 = 100 * (mean_total_585 - mean_total_baseline) / mean_total_baseline,
    pct_hur_585 = 100 * (mean_hur_585 - mean_hur_baseline) / mean_hur_baseline
  ) %>%
  select(location, pct_total_585, pct_hur_585)

cat("\n--- SSP5-8.5 change relative to baseline (%) ---\n")
print(
  change_585 %>%
    mutate(across(where(is.numeric), ~ round(.x, 1)))
)

# =============================================================================
# 8) Daily hazard illustration for one island
# =============================================================================
# Generate a daily synthetic series for Saba.
#
# Important API change:
# generate_daily_hazard_impact() now:
#   - uses `location` instead of `island`
#   - uses `scenario` instead of `cc_scenario`
#   - returns a named list, so extract the tibble with $Saba

daily_baseline <- generate_daily_hazard_impact(
  out = out_baseline,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage_method = "powerlaw",
  scenario = "baseline",
  seed = seed
)$Saba

daily_585 <- generate_daily_hazard_impact(
  out = out_585,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage_method = "powerlaw",
  scenario = "ssp585",
  seed = seed
)$Saba

daily_compare <- bind_rows(daily_baseline, daily_585)

daily_summary <- daily_compare %>%
  group_by(scenario) %>%
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
  daily_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

# =============================================================================
# 9) Optional quick plot
# =============================================================================
# Example plot of total annual storm activity by scenario.

p_activity <- activity_summary %>%
  ggplot(aes(x = location, y = mean_total, fill = scenario)) +
  geom_col(position = "dodge") +
  labs(
    x = NULL,
    y = "Mean annual storm count",
    title = "Simulated annual storm activity by scenario"
  ) +
  theme_light(base_size = 11)

print(p_activity)
