

# Temporal Downscaling ######################################

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
  severities = c("TS", "HUR")
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




#Histogram of yday(date) where wind_kt >= 34 should resemble historical.