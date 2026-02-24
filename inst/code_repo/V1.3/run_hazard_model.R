# run_hazard_model.R
# Config + per-island thresholds + run pipeline

rm(list = ls())

library(dplyr)
library(tibble)

source("R/V1.3/functions_hazard.R")

set.seed(123)

# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

cfg <- list(
  # Data
  ibtracs_path = "./data/ibtracs.NA.list.v04r01.csv",
  min_year = 1970,
  
  # Computational gate (NOT impact definition)
  gate_km = 800,
  
  # Wind thresholds (kt) used for hazard severity classification
  thr_ts = 34,
  thr_50 = 50,
  thr_hur = 64,
  
  # Radii parsing caps (nautical miles)
  cap_r34_nm = 600,
  cap_r50_nm = 400,
  cap_r64_nm = 250,
  
  # Simulation
  n_years_sim = 1000
)

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

# -----------------------------------------------------------------------------
# Per-island config: thresholds for SD coupling (ports, critical infra, etc.)
# (All units: knots sustained at site)
# -----------------------------------------------------------------------------
# You can extend each list with additional parameters (e.g. "thr_power", "thr_road", ...)
# Anything not provided defaults to NA and produces NA disruption flags.

per_target_cfg <- list(
  Saba = list(
    thr_port  = 40,  # example: port operational disruption threshold
    thr_infra = 55   # example: critical infrastructure disruption threshold
  ),
  Statia = list(
    thr_port  = 38,
    thr_infra = 52
  ),
  St_Martin = list(
    thr_port  = 45,
    thr_infra = 60
  )
  # Puerto_Rico / Miami can be added similarly
)

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

out <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  severities = c("TS", "HUR")
)

# -----------------------------------------------------------------------------
# Primary tidy output for SD coupling:
# one row per (island, storm SID, season) with site wind metrics + disruption flags
# -----------------------------------------------------------------------------

events_all <- out$events_all

# Example: inspect columns
print(dplyr::glimpse(events_all))

# Optional: write tidy table
# readr::write_csv(events_all, "outputs/hazard_events_all.csv")

# Optional: island-level lambdas
print(out$lambda_all)

# Optional: two-level simulation results per island
# readr::write_csv(out$sim_all, "outputs/sim_twolevel_all.csv")

ev <- out$events_all  # or out$events_by_island[["Saba"]]

# 1) What year window is actually present for the severities you count?
ev %>%
  filter(severity %in% c("TS", "HUR")) %>%
  summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    n_distinct_years = n_distinct(year)
  )

# 2) Inspect annual counts object directly
ac <- out$annual_counts_all %>% filter(island == "Saba")
ac %>% summarise(min_year = min(year), max_year = max(year), n_years = n())

# 3) Cross-check: mean counts from ac should reproduce lambda
ac %>%
  group_by(severity) %>%
  summarise(lambda_check = mean(n_events), n_years = n(), .groups = "drop")



