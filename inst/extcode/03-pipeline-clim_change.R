
# =============================================================================
# Pipeline assessment with Level 1+2+3 Climate Change Modifications
#
# Level 1: SST-conditioned rate scaling (beta_SST)
# Level 2: Intensity distribution shift (gamma)
# Level 3: Storm characteristic perturbation (V_peak, RMW, duration, precip)
#
# Demonstrates:
#   1. Stationary baseline (no climate mods)
#   2. Level 1 only (rate scaling, constant severity split)
#   3. Level 1+2 (rate scaling + intensity shift) under SSP scenarios
#   4. Severity split comparison across scenarios
#   5. Daily impact series with cc_scenario label
#   6. Level 1+2+3 (full perturbation) comparison
# =============================================================================

rm(list = ls())

# =============================================================================
# Configuration
# =============================================================================

# library(tidyr)
# library(stringr)
# library(lubridate)
# library(geosphere)
# library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ipdcstorm)

cfg <- make_hazard_cfg(
  data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  search_radius_km = 800,
  start_year = 1970,
  n_sim_years = 1000
)

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918)

per_target_cfg <- list(
  Saba    = list(thr_port = 40, thr_infra = 55),
  Statia  = list(thr_port = 38, thr_infra = 52),
  St_Martin = list(thr_port = 45, thr_infra = 60))


# =============================================================================
# 1) STATIONARY BASELINE
# =============================================================================

seed_compare <- 123
set.seed(seed_compare)
out_stat <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = NULL)

#print(out_stat$lambda_all)

# =============================================================================
# 2) Rate scaling + intensity shift under SSP scenarios
# =============================================================================

# Rate scaling + intensity shift under SSP scenarios

sst_245 <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  scenario = "ssp245",
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,
    beta_prior = 0.6,
    gamma_intensity = NULL,   # estimate from data
    gamma_prior = 0.065,      # 6.5% HUR fraction increase per C
    cc_params = list(
      v_scale     =  NULL,     #0.05,  # +5% peak intensity per degC
      r_scale     =  NULL,     #0.08,  # +8% radii expansion per degC
      speed_scale =  NULL,     #-0.10, # -10% translation speed per degC
      precip_scale = NULL      #0.07   # +7% rainfall per degC (CC scaling)
    )
  )
)

set.seed(seed_compare)
out_245 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_245)

sst_585 <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  scenario = "ssp585",
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,
    beta_prior = 0.6,
    gamma_intensity = NULL,
    gamma_prior = 0.065
  )
)

set.seed(seed_compare)
out_585 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_585
)

# =============================================================================
# 4) COMPARE ACROSS SCENARIOS
# =============================================================================

sim_compare <- bind_rows(
  out_stat$sim |> mutate(scenario = "Stationary"),
  out_245$sim  |> mutate(scenario = "SSP2-4.5"),
  out_585$sim  |> mutate(scenario = "SSP5-8.5"))

# --- Activity summary ---
activity_summary <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    mean_total = mean(n_total),
    mean_ts = mean(n_ts),
    mean_hur = mean(n_hur),
    mean_p_hur = mean(n_hur / pmax(1, n_total)),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

cat("\n--- Mean Annual Activity by Scenario ---\n")
print(knitr::kable(activity_summary |> mutate(across(where(is.numeric), ~round(., 3))),
                   format = "pipe"))

# --- Severity split comparison  ---
split_summary <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    p_hur_sim = mean(n_hur / pmax(1, n_total)),
    mean_p_hur_param = mean(p_hurricane),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

cat("\n--- Hurricane Fraction by Scenario (L2 Effect) ---\n")
print(knitr::kable(split_summary |> mutate(across(where(is.numeric), ~round(., 4))),
                   format = "pipe"))

# --- Percentage changes vs stationary ---
change_table <- sim_compare |>
  group_by(scenario, location) |>
  summarise(
    mean_total = mean(n_total),
    mean_hur = mean(n_hur),
    .groups = "drop"
  )

base_table <- change_table |>
  filter(scenario == "Stationary") |>
  transmute(location, mean_total_stationary = mean_total, mean_hur_stationary = mean_hur)

ssp585_table <- change_table |>
  filter(scenario == "SSP5-8.5") |>
  transmute(location, mean_total_585 = mean_total, mean_hur_585 = mean_hur)

change_table <- base_table |>
  left_join(ssp585_table, by = "location") |>
  mutate(
    pct_total_585 = 100 * (mean_total_585 - mean_total_stationary) / mean_total_stationary,
    pct_hur_585 = 100 * (mean_hur_585 - mean_hur_stationary) / mean_hur_stationary
  )

print(knitr::kable(
  change_table |>
    select(location, pct_total_585, pct_hur_585) |>
    mutate(across(where(is.numeric), ~round(., 1))),
  format = "pipe"))


# =============================================================================
# 5) DAILY IMPACT SERIES WITH cc_scenario LABEL
# =============================================================================

message("\n", strrep("=", 72))
message("  5) DAILY IMPACT SERIES (Saba)")
message(strrep("=", 72))

daily_stat <- generate_daily_hazard_impact(
  out = out_stat, island = "Saba", sim_years = 1:200, year0 = 2025,
  thr_port = 40, thr_infra = 55, gust_factor = 1.25,
  damage_method = "powerlaw", seed = 42,
  cc_scenario = "stationary"
)

daily_585 <- generate_daily_hazard_impact(
  out = out_585, island = "Saba", sim_years = 1:200, year0 = 2025,
  thr_port = 40, thr_infra = 55, gust_factor = 1.25,
  damage_method = "powerlaw", seed = 42,
  cc_scenario = "ssp585"
)

# Combine for comparison
daily_both <- bind_rows(daily_stat, daily_585)

# Summary by scenario
daily_summary <- daily_both |>
  group_by(cc_scenario) |>
  summarise(
    n_days = n(),
    tc_days = sum(wind_kt > 0),
    ts_days = sum(wind_kt >= 34),
    hur_days = sum(wind_kt >= 64),
    port_disrupt_days = sum(port_disrupt, na.rm = TRUE),
    mean_annual_damage = sum(damage_rate, na.rm = TRUE) / 200,
    .groups = "drop"
  )

cat("\nDaily impact comparison (Saba, 200 years):\n")
print(knitr::kable(daily_summary |> mutate(across(where(is.numeric), ~round(., 2))),
                   format = "pipe"))


# =============================================================================
# 6) DIAGNOSTIC PLOTS
# =============================================================================

if (nrow(out_585$fit) > 0) {

  dir.create("output/climate", recursive = TRUE, showWarnings = FALSE)

  # --- Scenario trajectories ---
  scenarios <- bind_rows(
    generate_sst_scenario(76, mode = "stationary", start_year = 2025),
    generate_sst_scenario(76, mode = "ssp245", start_year = 2025),
    generate_sst_scenario(76, mode = "ssp585", start_year = 2025)
  )

  fit_row <- out_585$fit |>
    select(beta_sst, gamma_intensity, p_hurricane_base) |>
    distinct() |>
    slice(1)

  gamma <- fit_row$gamma_intensity
  p_base <- fit_row$p_hurricane_base
  scenarios <- scenarios |>
    mutate(p_hur_t = pmin(0.99, pmax(0.01, p_base * (1 + gamma * sst_anomaly))))

  p3 <- ggplot(scenarios, aes(x = calendar_year)) +
    geom_line(aes(y = p_hur_t, color = scenario), linewidth = 0.8) +
    geom_hline(yintercept = p_base, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      stationary = "grey40", ssp245 = "orange", ssp585 = "red"
    )) +
    labs(
      x = "Year", y = "p_HUR (hurricane fraction)",
      title = "Projected Hurricane Fraction Under SSP Scenarios",
      subtitle = sprintf("gamma=%.4f, p_HUR_base=%.3f", gamma, p_base),
      color = "Scenario"
    ) +
    theme_light(base_size = 11)

  ggsave("output/climate/L2_p_hur_projections.png", p3, width = 8, height = 5, dpi = 150)
  message("Saved: output/climate/L2_p_hur_projections.png")
}


# =============================================================================
# 7) LEVEL 3: STORM PERTURBATION COMPARISON
# =============================================================================

message("\n", strrep("=", 72))
message("  7) LEVEL 1+2+3 (SSP5-8.5 with storm perturbation)")
message(strrep("=", 72))

sst_585_L3 <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  scenario = "ssp585",
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,
    beta_prior = 0.6,
    gamma_intensity = NULL,
    gamma_prior = 0.065,
    cc_params = list(
      v_scale     =  0.05,   # +5% peak intensity per degC
      r_scale     =  0.08,   # +8% radii expansion per degC
      speed_scale = -0.10,   # -10% translation speed per degC
      precip_scale =  0.07   # +7% rainfall per degC (CC scaling)
    )
  )
)

set.seed(seed_compare)
out_585_L3 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_585_L3
)

cat("\nL3 cc_params:\n")
print(out_585_L3$cc_params)

# --- Daily comparison: L1+L2 vs L1+L2+L3 ---
message("\n  Generating daily series: L1+L2 vs L1+L2+L3 (Saba, 200 years)")

daily_L12 <- generate_daily_hazard_impact(
  out = out_585, island = "Saba", sim_years = 1:200, year0 = 2025,
  thr_port = 40, thr_infra = 55, gust_factor = 1.25,
  damage_method = "powerlaw", seed = 42,
  cc_scenario = "L1+L2_ssp585"
)

daily_L123 <- generate_daily_hazard_impact(
  out = out_585_L3, island = "Saba", sim_years = 1:200, year0 = 2025,
  thr_port = 40, thr_infra = 55, gust_factor = 1.25,
  damage_method = "powerlaw", seed = 42,
  cc_scenario = "L1+L2+L3_ssp585"
)

daily_L3_compare <- bind_rows(daily_stat, daily_L12, daily_L123)

L3_summary <- daily_L3_compare |>
  group_by(cc_scenario) |>
  summarise(
    n_days = n(),
    ts_days = sum(wind_kt >= 34),
    hur_days = sum(wind_kt >= 64),
    cat3plus_days = sum(wind_kt >= 96),
    mean_peak_wind = mean(peak_wind_year_kt, na.rm = TRUE),
    mean_annual_damage = sum(damage_rate, na.rm = TRUE) / 200,
    .groups = "drop"
  )

cat("\nLevel 3 impact comparison (Saba, 200 years):\n")
print(knitr::kable(L3_summary |> mutate(across(where(is.numeric), ~round(., 2))),
                   format = "pipe"))


# =============================================================================
# 8) LEVEL 3 SENSITIVITY ANALYSIS: delta_SST = {0, +1, +2, +3}
# =============================================================================

message("\n", strrep("=", 72))
message("  8) L3 SENSITIVITY: monotonicity check at fixed delta_SST")
message(strrep("=", 72))

# Quick validation: run sampling + perturbation at fixed delta_SST values
# and verify monotonic increase in intensity metrics

lib_saba <- build_event_library_from_out(out_stat, island = "Saba", seed = 42)

set.seed(99)
sens_results <- list()

for (dsst in c(0, 1, 2, 3)) {
  events_list_sens <- list()
  for (rep in 1:500) {
    sampled <- sample_events_for_year_extended(
      lib = lib_saba, year = 2050, n_ts = 2, n_hur = 1
    )
    sampled <- perturb_event(sampled, delta_sst = dsst, cc_params = default_cc_params())
    events_list_sens[[rep]] <- sampled
  }
  all_events <- bind_rows(events_list_sens)
  sens_results[[paste0("dSST_", dsst)]] <- tibble(
    delta_sst = dsst,
    mean_V_peak = mean(all_events$V_peak, na.rm = TRUE),
    p90_V_peak = quantile(all_events$V_peak, 0.9, na.rm = TRUE),
    mean_dur = mean(all_events$dur_days, na.rm = TRUE),
    mean_precip_scale = mean(all_events$precip_scaling, na.rm = TRUE)
  )
}

sens_df <- bind_rows(sens_results)

cat("\nL3 Sensitivity (500 reps, n_ts=2 + n_hur=1 per year):\n")
print(knitr::kable(sens_df |> mutate(across(where(is.numeric), ~round(., 2))),
                   format = "pipe"))

# Monotonicity check
mono_V <- all(diff(sens_df$mean_V_peak) > 0)
mono_d <- all(diff(sens_df$mean_dur) > 0)
mono_p <- all(diff(sens_df$mean_precip_scale) > 0)

cat(sprintf("\nMonotonicity check: V_peak=%s, dur=%s, precip=%s\n",
            ifelse(mono_V, "PASS", "FAIL"),
            ifelse(mono_d, "PASS", "FAIL"),
            ifelse(mono_p, "PASS", "FAIL")))

# Identity check at delta_SST = 0
id_check <- sens_df$mean_precip_scale[sens_df$delta_sst == 0] == 1.0
cat(sprintf("Identity check (dSST=0): precip_scaling=1.0? %s\n",
            ifelse(id_check, "PASS", "FAIL")))


message("\n>>> Pipeline assessment complete.")
