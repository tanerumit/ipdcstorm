
# =============================================================================
# Script overview: runnable workflow for hazard modeling and daily impacts.
# - No function definitions in this script.
# - Uses functions sourced from hazard_core.R, hazard_ibtracs.R, hazard_run.R,
#   hazard_downscale.R, and hazard_utils.R to run the model, summarize results,
#   and generate daily synthetic series.
# =============================================================================

### RUN HAZARD MODEL AND ESTIMATE DAILY DAMAGE

rm(list = ls())

library(tidyr)
library(stringr)
library(lubridate)
library(geosphere)
library(readr)
library(dplyr)
library(tibble)

#source("functions_hazard.R")
source("R/hazard_core.R")
source("R/hazard_downscale.R")
source("R/hazard_ibtracs.R")
source("R/hazard_utils.R")
source("R/hazard_run.R")



set.seed(123)



# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

cfg <- list(
  # Data
  ibtracs_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
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
  severities = c("TC", "HUR")
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
  filter(severity %in% c("TC", "HUR")) %>%
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





# Temporal Downscaling ##########################################################

# Simple duration first:
#June 1 – Nov 30 (DOY 152–334)

#TS: 1–3 days
#HUR: 1–2 days (site exceedance is often brief),

# Run for single location
saba_target <- targets %>% filter(name == "Saba")

out <- run_hazard_model(
  cfg = cfg,
  targets = saba_target,
  per_target_cfg = per_target_cfg,
  severities = c("TC", "HUR")
)

# Lamda table
out$lambda_all
out$events_all
View(out$events_all)
knitr::kable(out$events_all, format = "pipe")

daily_saba <- generate_daily_from_hazard(
  out = out,
  island = "Saba",
  sim_years = 1:200, # 200 synthetic years
  year0 = 2000,
  thr_port = 40,
  thr_infra = 55,
  seed = 42
)

# Example: apply to your daily_saba
daily_saba2 <- daily_saba |>
  dplyr::mutate(
    damage_rate = damage_rate_from_wind(
      wind_kt,
      thr = 34, V_ref = 80, d_ref = 0.03, p = 3, d_max = 0.10
    ),
    cum_damage = cumsum(dplyr::coalesce(damage_rate, 0))
)




