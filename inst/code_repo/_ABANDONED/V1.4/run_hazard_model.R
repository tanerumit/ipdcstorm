# run_hazard_model_enhanced.R
# 
# ENHANCEMENTS:
#   1. Refined Saffir-Simpson severity classification (7 classes: none, TS, CAT1-5)
#   2. Integrated Wind Exposure Index for cumulative damage assessment
#   3. Flexible severity scheme selection (Saffir-Simpson, simple, or major)
#   4. Peak gust estimates and multiple exposure metrics

rm(list = ls())

library(dplyr)
library(tibble)
library(ggplot2)

source("R/V1.4/functions_hazard.R")

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
# Per-island config: thresholds for infrastructure disruption
# -----------------------------------------------------------------------------

per_target_cfg <- list(
  Saba = list(
    thr_port  = 40,  # port operational disruption threshold
    thr_infra = 55   # critical infrastructure disruption threshold
  ),
  Statia = list(
    thr_port  = 38,
    thr_infra = 52
  ),
  St_Martin = list(
    thr_port  = 45,
    thr_infra = 60
  )
)

# -----------------------------------------------------------------------------
# Run model with Saffir-Simpson classification
# -----------------------------------------------------------------------------

cat("\n=== Running Enhanced Model with Saffir-Simpson Classification ===\n\n")

out_saffir <- run_hazard_model_enhanced(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  severity_scheme = "saffir_simpson"
  # This will use: TS, CAT1, CAT2, CAT3, CAT4, CAT5
)

# -----------------------------------------------------------------------------
# Examine enhanced outputs
# -----------------------------------------------------------------------------

cat("Enhanced event columns (first 3 rows):\n")
print(head(out_saffir$events_all %>% select(island, year, name, 
                                            V_site_max_kt, V_gust_max_kt,
                                            wind_exposure_index, hours_ge34, 
                                            severity), 3))

cat("\n\nSeverity distribution across all islands:\n")
print(out_saffir$events_all %>% 
        filter(severity != "none") %>%
        count(severity) %>%
        arrange(desc(n)))

cat("\n\nAnnual occurrence rates (lambda) by severity class:\n")
print(out_saffir$lambda_all %>% 
        arrange(island, desc(lambda)))

# -----------------------------------------------------------------------------
# DEMONSTRATION: Compare integrated wind exposure vs. simple peak wind
# -----------------------------------------------------------------------------

cat("\n\n=== Wind Exposure Index vs. Peak Wind Analysis ===\n")

exposure_analysis <- out_saffir$events_all %>%
  filter(severity %in% c("TS", "CAT1", "CAT2", "CAT3", "CAT4", "CAT5")) %>%
  group_by(island, severity) %>%
  summarise(
    n_events = n(),
    mean_peak_wind = mean(V_site_max_kt, na.rm = TRUE),
    mean_exposure_index = mean(wind_exposure_index, na.rm = TRUE),
    mean_duration_hrs = mean(hours_ge34, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(island, desc(mean_exposure_index))

print(exposure_analysis)

# -----------------------------------------------------------------------------
# DEMONSTRATION: Alternative severity schemes
# -----------------------------------------------------------------------------

cat("\n\n=== Running with simplified 3-class scheme ===\n")

out_simple <- run_hazard_model_enhanced(
  cfg = cfg,
  targets = targets %>% filter(name == "Saba"),  # Just Saba for demo
  severity_scheme = "simple"
  # This will use: TS, HUR64plus
)

cat("\nSimple scheme lambdas (Saba only):\n")
print(out_simple$lambda_all)

cat("\n\n=== Running with 4-class major hurricane scheme ===\n")

out_major <- run_hazard_model_enhanced(
  cfg = cfg,
  targets = targets %>% filter(name == "Saba"),
  severity_scheme = "major"
  # This will use: TS, HUR (Cat 1-2), MAJOR (Cat 3+)
)

cat("\nMajor scheme lambdas (Saba only):\n")
print(out_major$lambda_all)

# -----------------------------------------------------------------------------
# DEMONSTRATION: Top 10 most damaging events (by integrated exposure)
# -----------------------------------------------------------------------------

cat("\n\n=== Top 10 Most Damaging Events by Wind Exposure Index ===\n")

top_events <- out_saffir$events_all %>%
  filter(!is.na(wind_exposure_index)) %>%
  arrange(desc(wind_exposure_index)) %>%
  select(island, year, name, severity, V_site_max_kt, V_gust_max_kt,
         wind_exposure_index, hours_ge34, hours_ge64) %>%
  head(10)

print(top_events)

# -----------------------------------------------------------------------------
# DEMONSTRATION: Event comparison - same storm, different metrics
# -----------------------------------------------------------------------------

cat("\n\n=== Event Comparison: Peak Wind vs. Integrated Exposure ===\n")
cat("Events where exposure index ranking differs from peak wind ranking\n\n")

# Find events where high exposure doesn't match high peak wind (long-duration moderate events)
comparison <- out_saffir$events_all %>%
  filter(severity %in% c("TS", "CAT1", "CAT2")) %>%
  mutate(
    wind_rank = rank(-V_site_max_kt, ties.method = "average"),
    exposure_rank = rank(-wind_exposure_index, ties.method = "average"),
    rank_diff = abs(wind_rank - exposure_rank)
  ) %>%
  filter(rank_diff > 10) %>%  # Significant ranking differences
  arrange(desc(rank_diff)) %>%
  select(island, year, name, V_site_max_kt, hours_ge34, wind_exposure_index,
         wind_rank, exposure_rank, rank_diff) %>%
  head(5)

if (nrow(comparison) > 0) {
  print(comparison)
  cat("\nInterpretation: These events had moderate peak winds but long exposure duration,\n")
  cat("resulting in high cumulative damage potential despite lower peak intensity.\n")
} else {
  cat("No significant ranking differences found in this sample.\n")
}

# -----------------------------------------------------------------------------
# Output files
# -----------------------------------------------------------------------------

# cat("\n\n=== Writing output files ===\n")
# 
# # Main events table with all enhanced metrics
# readr::write_csv(out_saffir$events_all, "/mnt/user-data/outputs/hazard_events_enhanced.csv")
# cat("✓ Written: hazard_events_enhanced.csv\n")
# 
# # Annual occurrence rates by severity
# readr::write_csv(out_saffir$lambda_all, "/mnt/user-data/outputs/lambda_rates_saffir_simpson.csv")
# cat("✓ Written: lambda_rates_saffir_simpson.csv\n")
# 
# # Simulation results
# readr::write_csv(out_saffir$sim_all, "/mnt/user-data/outputs/simulation_saffir_simpson.csv")
# cat("✓ Written: simulation_saffir_simpson.csv\n")
# 
# # Exposure analysis summary
# readr::write_csv(exposure_analysis, "/mnt/user-data/outputs/exposure_analysis_summary.csv")
# cat("✓ Written: exposure_analysis_summary.csv\n")

# -----------------------------------------------------------------------------
# VISUALIZATION
# -----------------------------------------------------------------------------

library(ggplot2)

# Plot 1: Lambda rates by severity class
p1 <- ggplot(out_saffir$lambda_all, aes(x = severity, y = lambda, fill = island)) +
geom_col(position = "dodge") +
labs(
  title = "Annual Occurrence Rates by Severity Class",
  subtitle = "Saffir-Simpson Classification (1970-present)",
  x = "Severity Class",
  y = "Lambda (events/year)",
  fill = "Island"
) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave("/lambda_rates_plot.png", p1, width = 10, height = 6, dpi = 300)
cat("✓ Created: lambda_rates_plot.png\n")

# Plot 2: Wind Exposure Index vs Peak Wind (scatter)
p2 <- out_saffir$events_all %>%
filter(severity != "none", !is.na(wind_exposure_index)) %>%
ggplot(aes(x = V_site_max_kt, y = wind_exposure_index, color = severity)) +
geom_point(alpha = 0.6, size = 2) +
scale_color_manual(
  values = c("TS" = "#74add1", "CAT1" = "#fee090", "CAT2" = "#fdae61",
             "CAT3" = "#f46d43", "CAT4" = "#d73027", "CAT5" = "#a50026")
) +
labs(
  title = "Integrated Wind Exposure Index vs. Peak Wind Speed",
  subtitle = "Events vary in damage potential based on both intensity and duration",
  x = "Peak Wind Speed (kt)",
  y = "Wind Exposure Index (kt³·hrs / 1000)",
  color = "Severity"
) +
theme_minimal() +
facet_wrap(~island, scales = "free")
