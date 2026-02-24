
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

library(tidyr)
library(stringr)
library(lubridate)
library(geosphere)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

source("R/hazard_core.R")
source("R/hazard_climate.R")
source("R/hazard_downscale.R")
source("R/hazard_ibtracs.R")
source("R/hazard_utils.R")
source("R/hazard_run.R")

set.seed(123)


# =============================================================================
# Global configuration
# =============================================================================

cfg <- make_hazard_cfg(
  data_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
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
  "Miami",       25.7617,  -80.1918
)

per_target_cfg <- list(
  Saba    = list(thr_port = 40, thr_infra = 55),
  Statia  = list(thr_port = 38, thr_infra = 52),
  St_Martin = list(thr_port = 45, thr_infra = 60)
)


# =============================================================================
# 1) STATIONARY BASELINE
# =============================================================================

message("\n", strrep("=", 72))
message("  1) STATIONARY BASELINE")
message(strrep("=", 72))

out_stat <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = NULL
)

cat("\nStationary lambdas:\n")
print(out_stat$lambda_all)


# =============================================================================
# 2) LEVEL 1 ONLY: rate scaling, no intensity shift
# =============================================================================

message("\n", strrep("=", 72))
message("  2) LEVEL 1 ONLY (SSP5-8.5 rate scaling, gamma=0)")
message(strrep("=", 72))

sst_L1_only <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  scenario = "ssp585",
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,         # estimate from data
    beta_prior = 0.6,
    gamma_intensity = 0,     # explicitly disable intensity shift
    gamma_prior = 0
  )
)

out_L1 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_L1_only
)

cat("\nbeta_SST:", out_L1$beta_sst, "\n")
cat("gamma:", out_L1$gamma_intensity, "\n")


# =============================================================================
# 3) LEVEL 1+2: rate scaling + intensity shift under SSP scenarios
# =============================================================================

message("\n", strrep("=", 72))
message("  3) LEVEL 1+2 (SSP2-4.5)")
message(strrep("=", 72))

sst_245 <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",
  baseline_years = 1991L:2020L,
  scenario = "ssp245",
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,
    beta_prior = 0.6,
    gamma_intensity = NULL,  # estimate from data
    gamma_prior = 0.065      # 6.5% HUR fraction increase per C
  )
)

out_245 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_245
)

message("\n", strrep("=", 72))
message("  3b) LEVEL 1+2 (SSP5-8.5)")
message(strrep("=", 72))

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

out_585 <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = sst_585
)

cat("\nEstimated climate parameters:\n")
cat("  beta_SST:", out_585$beta_sst, "\n")
cat("  gamma:", out_585$gamma_intensity, "\n")
cat("  p_HUR_base:", out_585$p_hur_base, "\n")


# =============================================================================
# 4) COMPARE ACROSS SCENARIOS
# =============================================================================

message("\n", strrep("=", 72))
message("  4) SCENARIO COMPARISON")
message(strrep("=", 72))

sim_compare <- bind_rows(
  out_stat$sim_all |> mutate(scenario = "Stationary"),
  out_L1$sim_all   |> mutate(scenario = "L1 only (SSP5-8.5)"),
  out_245$sim_all  |> mutate(scenario = "L1+L2 SSP2-4.5"),
  out_585$sim_all  |> mutate(scenario = "L1+L2 SSP5-8.5")
)

# --- Activity summary ---
activity_summary <- sim_compare |>
  group_by(scenario, island) |>
  summarise(
    mean_total = mean(n_total),
    mean_TS = mean(n_TS),
    mean_HUR = mean(n_HUR),
    mean_p_hur = mean(n_HUR / pmax(1, n_total)),
    mean_A = mean(A),
    .groups = "drop"
  ) |>
  arrange(island, scenario)

cat("\n--- Mean Annual Activity by Scenario ---\n")
print(knitr::kable(activity_summary |> mutate(across(where(is.numeric), ~round(., 3))),
                   format = "pipe"))

# --- Severity split comparison (L2 effect) ---
split_summary <- sim_compare |>
  group_by(scenario, island) |>
  summarise(
    p_hur_sim = mean(n_HUR / pmax(1, n_total)),
    mean_p_hur_param = mean(p_hur),
    .groups = "drop"
  ) |>
  arrange(island, scenario)

cat("\n--- Hurricane Fraction by Scenario (L2 Effect) ---\n")
print(knitr::kable(split_summary |> mutate(across(where(is.numeric), ~round(., 4))),
                   format = "pipe"))

# --- Percentage changes vs stationary ---
change_table <- sim_compare |>
  group_by(scenario, island) |>
  summarise(
    mean_total = mean(n_total),
    mean_HUR = mean(n_HUR),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = scenario,
              values_from = c(mean_total, mean_HUR)) |>
  mutate(
    pct_total_L1L2_585 = 100 * (`mean_total_L1+L2 SSP5-8.5` - mean_total_Stationary) / mean_total_Stationary,
    pct_HUR_L1L2_585 = 100 * (`mean_HUR_L1+L2 SSP5-8.5` - mean_HUR_Stationary) / mean_HUR_Stationary,
    pct_HUR_L1only_585 = 100 * (`mean_HUR_L1 only (SSP5-8.5)` - mean_HUR_Stationary) / mean_HUR_Stationary
  )

cat("\n--- Percentage Change vs Stationary (SSP5-8.5) ---\n")
cat("  (L1+L2 captures both more activity AND higher HUR fraction)\n")
print(knitr::kable(
  change_table |>
    select(island, pct_total_L1L2_585, pct_HUR_L1only_585, pct_HUR_L1L2_585) |>
    mutate(across(where(is.numeric), ~round(., 1))),
  format = "pipe"
))


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

if (!is.null(out_585$sst_info$gamma_info$fit_data)) {

  dir.create("output/climate", recursive = TRUE, showWarnings = FALSE)

  # --- SST-activity plot (L1) ---
  if (!is.null(out_585$sst_info$beta_info$fit_data)) {
    fit_df <- out_585$sst_info$beta_info$fit_data

    p1 <- ggplot(fit_df, aes(x = sst_anomaly, y = N)) +
      geom_point(color = "steelblue", size = 2.5, alpha = 0.7) +
      geom_smooth(method = "glm", method.args = list(family = "poisson"),
                  formula = y ~ x, color = "red", linewidth = 0.8, se = TRUE) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = "MDR SST Anomaly (C)",
        y = "Annual TC count",
        title = "L1: SST-Activity Relationship",
        subtitle = sprintf("beta_SST = %.3f | +1C -> %+.0f%% activity",
                            out_585$beta_sst,
                            100 * (exp(out_585$beta_sst) - 1))
      ) +
      theme_light(base_size = 11)

    ggsave("output/climate/L1_sst_activity.png", p1, width = 8, height = 5, dpi = 150)
    message("Saved: output/climate/L1_sst_activity.png")
  }

  # --- HUR fraction vs SST plot (L2) ---
  gamma_df <- out_585$sst_info$gamma_info$fit_data

  p2 <- ggplot(gamma_df, aes(x = sst_anomaly, y = p_hur)) +
    geom_point(color = "darkred", size = 2.5, alpha = 0.7) +
    geom_smooth(method = "glm",
                method.args = list(family = binomial(link = "logit")),
                formula = y ~ x, color = "red", linewidth = 0.8, se = TRUE) +
    geom_hline(yintercept = out_585$p_hur_base, linetype = "dashed", color = "grey50") +
    labs(
      x = "MDR SST Anomaly (C)",
      y = "Hurricane fraction (n_HUR / n_total)",
      title = "L2: Intensity Shift - HUR Fraction vs SST",
      subtitle = sprintf("gamma = %.4f | p_HUR_base = %.3f | +1C -> p_HUR %+.1f%%",
                          out_585$gamma_intensity,
                          out_585$p_hur_base,
                          100 * out_585$gamma_intensity)
    ) +
    theme_light(base_size = 11)

  ggsave("output/climate/L2_intensity_shift.png", p2, width = 8, height = 5, dpi = 150)
  message("Saved: output/climate/L2_intensity_shift.png")


  # --- Scenario trajectories ---
  scenarios <- bind_rows(
    generate_sst_scenario(76, mode = "stationary", start_year = 2025),
    generate_sst_scenario(76, mode = "ssp245", start_year = 2025),
    generate_sst_scenario(76, mode = "ssp585", start_year = 2025)
  )

  gamma <- out_585$gamma_intensity
  p_base <- out_585$p_hur_base
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
