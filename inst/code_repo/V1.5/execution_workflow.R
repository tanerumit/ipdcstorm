################################################################################
# Workflow: Saba (end-to-end) — hazard -> annual simulation -> daily winds -> damage
################################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(tibble)
library(geosphere)

# Set the path of ibtracs dataset
ibtracs_path <- "./data/ibtracs/ibtracs.NA.list.v04r01.csv"

# Load functions
source("R/V1.5/functions_hazard.R")

# set randomization seed 
set.seed(123)

# -----------------------------------------------------------------------------
# 1) Configuration and target
# -----------------------------------------------------------------------------

# Wind thresholds for Tropical Cyclones (TC) and Hurricanes, kt
thr_tc  <- 34
thr_hur <- 64


cfg <- list(
  ibtracs_path = ibtracs_path,
  min_year = 1970,
  gate_km = 800,
  thr_ts = 34,
  thr_hur = 64,
  n_years_sim = 500  # keep small for an example
)

targets <- tibble::tibble(
  name = "Saba",
  lat  = 17.6350,
  lon  = -63.2300
)

per_target_cfg <- list(
  Saba = list(
    thr_port  = 40,
    thr_infra = 55
  )
)

# -----------------------------------------------------------------------------
# 2) Run hazard model + build event library
# -----------------------------------------------------------------------------

out <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  severities = c("TS", "HUR64plus")
)

# Key outputs
events_saba <- out$events_by_island[["Saba"]]
ac_saba     <- out$annual_counts_all %>% filter(island == "Saba")
lambda_saba <- out$lambda_all %>% filter(island == "Saba")
sim_saba    <- out$sim_all %>% filter(island == "Saba")

# Inspect outputs
glimpse(events_saba)
print(lambda_saba)
print(out$kinfo_all)

# Quick sanity checks
events_saba %>%
  summarise(
    n_events = n(),
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    share_with_pressure = mean(is.finite(Pc_min_hPa)),
    share_with_rmw      = mean(is.finite(RMW_mean_km)))

events_saba %>% count(severity)

lib_saba <- build_event_library_from_out(out, island = "Saba", seed = 1)

# What DOY samples look like by class (TD/TS/HUR)
lib_saba$doy_by_severity

# -----------------------------------------------------------------------------
# 4) Generate synthetic daily wind+impact series from simulated annual counts
# -----------------------------------------------------------------------------

daily_saba_impact <- generate_daily_hazard_impact(
  out = out,
  island = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  thr_port = per_target_cfg$Saba$thr_port,
  thr_infra = per_target_cfg$Saba$thr_infra,
  damage_method = "intensity",     # or "powerlaw"
  damage_params = list(V0 = 34, V1 = 120, p = 3, dmax = 0.02),
  seed = 42
)

write_csv(daily_saba_impact, "output/saba_daily_impact.csv")


library(dplyr)

weekly_saba_impact <- daily_saba_impact %>%
  mutate(
    iso_weekday = as.integer(format(date, "%u")),          # 1 = Mon, ..., 7 = Sun
    week_start  = as.Date(date) - (iso_weekday - 1L),
    week_end    = week_start + 6L,
    week_num    = as.integer(format(date, "%V")),          # ISO week number
    week_year   = as.integer(format(date, "%G"))           # ISO week-based year
  ) %>%
  group_by(island, sim_year, week_year, week_num, week_start, week_end) %>%
  summarise(
    # hazard
    wind_max_kt  = max(wind_kt, na.rm = TRUE),
    wind_p95_kt  = quantile(wind_kt, 0.95, na.rm = TRUE),
    
    # events
    n_event_days  = sum(!is.na(event_id_day)),
    n_events_week = n_distinct(event_id_day[!is.na(event_id_day)]),
    any_tc_week   = any(is_tc_event_day %in% TRUE),
    any_hur_week  = any(is_hur_event_day %in% TRUE),
    
    # disruptions
    port_disrupt_any  = any(port_disrupt %in% TRUE),
    infra_disrupt_any = any(infra_disrupt %in% TRUE),
    
    # damage
    damage_rate_sum = sum(damage_rate, na.rm = TRUE),
    damage_inc_sum  = sum(damage_increment, na.rm = TRUE),
    cum_damage_end  = last(cum_damage),
    
    .groups = "drop"
  )%>% select(-c(port_disrupt_any, infra_disrupt_any, week_year, week_start, week_end))

write_csv(weekly_saba_impact, "output/saba_weekly_impact.csv")

glimpse(daily_saba_impact)


# Optional: write outputs
# write_csv(events_saba, "outputs/saba_events.csv")
# write_csv(sim_saba, "outputs/saba_sim_annual.csv")





 
