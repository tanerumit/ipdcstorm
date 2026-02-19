
### RUN VALIDATION SUITE
### Sources the hazard model, runs it, then executes all three validation tiers.
### Produces console report + diagnostic plots.

rm(list = ls())

library(tidyr)
library(stringr)
library(lubridate)
library(geosphere)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

# Source model code
source("R/hazard_core.R")
source("R/hazard_downscale.R")
source("R/hazard_ibtracs.R")
source("R/hazard_utils.R")
source("R/hazard_run.R")
source("R/hazard_validation.R")

set.seed(123)

# =============================================================================
# Configuration (same as run_hazard_model.R)
# =============================================================================

cfg <- list(
  ibtracs_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
  min_year     = 1970,
  gate_km      = 800,
  thr_ts       = 34,
  thr_50       = 50,
  thr_hur      = 64,
  cap_r34_nm   = 600,
  cap_r50_nm   = 400,
  cap_r64_nm   = 250,
  n_years_sim  = 1000
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
  Saba      = list(thr_port = 40, thr_infra = 55),
  Statia    = list(thr_port = 38, thr_infra = 52),
  St_Martin = list(thr_port = 45, thr_infra = 60)
)

# =============================================================================
# Run the hazard model
# =============================================================================

message("\n>>> Running hazard model...")
out <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  severities = c("TS", "HUR64plus")
)

# =============================================================================
# Run the full validation suite
# =============================================================================

message("\n>>> Running validation suite...")
val <- run_validation_suite(
  out = out,
  holdout_years = 10,
  n_sim = 5000,
  return_periods = c(5, 10, 25, 50),
  seed = 42
)

# =============================================================================
# Diagnostic Plots
# =============================================================================

dir.create("output/validation", recursive = TRUE, showWarnings = FALSE)

ggtheme <- theme_light(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))


# --- Plot 1: Return-level comparison (hindcast) ---

if (!is.null(val$hindcast) && nrow(val$hindcast$comparison) > 0) {

  p_rl <- ggplot(val$hindcast$comparison,
                 aes(x = factor(return_period), y = sim_rl)) +
    geom_errorbar(aes(ymin = sim_lo_90, ymax = sim_hi_90),
                  width = 0.3, color = "steelblue", linewidth = 0.8) +
    geom_point(aes(y = sim_median), size = 3, color = "steelblue", shape = 16) +
    geom_point(aes(y = obs_full_rl), size = 3, color = "red", shape = 17) +
    geom_hline(yintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    facet_wrap(~ island, scales = "free_y", ncol = 3) +
    labs(
      x = "Return period (years)",
      y = "Return level — peak site wind (kt)",
      title = "Hindcast Validation: Simulated vs Observed Return Levels",
      subtitle = "Blue = simulated median [90% CI]; Red triangle = observed (full record); dashed = 64 kt HUR threshold"
    ) +
    ggtheme

  ggsave("output/validation/hindcast_return_levels.png", p_rl,
         width = 12, height = 7, dpi = 150)
  message("Saved: output/validation/hindcast_return_levels.png")

  # --- Plot 1b: Simulated vs observed annual max distributions ---

  for (isl in names(val$hindcast$per_island)) {
    hc_isl <- val$hindcast$per_island[[isl]]
    if (is.null(hc_isl)) next

    obs_df <- hc_isl$obs_annual_max
    sim_df <- tibble::tibble(V_max_kt = hc_isl$sim_annual_max)

    p_dist <- ggplot() +
      geom_histogram(data = sim_df, aes(x = V_max_kt, y = after_stat(density)),
                     fill = "steelblue", alpha = 0.4, bins = 40) +
      geom_density(data = sim_df, aes(x = V_max_kt), color = "steelblue",
                   linewidth = 0.8) +
      geom_rug(data = obs_df, aes(x = V_max_kt, color = period),
               linewidth = 0.8, alpha = 0.8) +
      scale_color_manual(values = c(train = "grey40", test = "red")) +
      geom_vline(xintercept = 64, linetype = "dashed", color = "grey50") +
      labs(
        x = "Annual maximum site wind (kt)",
        y = "Density",
        title = paste("Annual Max Distribution:", isl),
        subtitle = "Blue = simulated; rug ticks = observed (grey=train, red=test)",
        color = "Observed"
      ) +
      ggtheme

    ggsave(sprintf("output/validation/hindcast_dist_%s.png", tolower(isl)), p_dist,
           width = 8, height = 5, dpi = 150)
  }
}


# --- Plot 2: Rate comparison ---

if (!is.null(val$rate_check) && nrow(val$rate_check) > 0) {

  rc <- val$rate_check |>
    filter(!is.na(lambda_ref))

  p_rate <- ggplot(rc, aes(x = lambda_ref, y = lambda_model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 2, intercept = 0, linetype = "dotted", color = "grey70") +
    geom_abline(slope = 0.5, intercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(aes(color = flag, shape = severity), size = 3.5) +
    geom_text(aes(label = island), nudge_y = 0.08, size = 3, check_overlap = TRUE) +
    scale_color_manual(values = c(
      "OK" = "forestgreen",
      "elevated (expected: model gate > ref radius)" = "orange",
      "HIGH: model >> reference (check gate_km or min_year)" = "red",
      "LOW: model << reference (check severity filter or data)" = "purple"
    )) +
    coord_equal(xlim = c(0, max(c(rc$lambda_model, rc$lambda_ref), na.rm = TRUE) + 0.3),
                ylim = c(0, max(c(rc$lambda_model, rc$lambda_ref), na.rm = TRUE) + 0.3)) +
    labs(
      x = "Reference λ (published climatology)",
      y = "Model λ (fitted)",
      title = "Rate Sanity Check: Model vs Published Annual Rates",
      subtitle = "Dashed = 1:1; dotted = 0.5x and 2x bounds",
      color = "Flag", shape = "Severity"
    ) +
    ggtheme

  ggsave("output/validation/rate_comparison.png", p_rate,
         width = 8, height = 6, dpi = 150)
  message("Saved: output/validation/rate_comparison.png")
}


# --- Plot 3: Wind field spot-checks ---

if (!is.null(val$wind_field) && any(is.finite(val$wind_field$model_V_site_kt))) {

  wf <- val$wind_field |>
    filter(is.finite(model_V_site_kt))

  p_wf <- ggplot(wf, aes(x = obs_1min_equiv_kt, y = model_V_site_kt)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 20, linetype = "dotted", color = "grey70") +
    geom_abline(slope = 1, intercept = -20, linetype = "dotted", color = "grey70") +
    geom_point(aes(color = island, shape = island), size = 3.5) +
    geom_text(aes(label = storm_name), nudge_y = 4, size = 3, check_overlap = TRUE) +
    coord_equal(xlim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20),
                ylim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20)) +
    labs(
      x = "Observed wind (1-min equiv, kt)",
      y = "Model V_site_kt (kt)",
      title = "Wind Field Validation: Model vs Station Observations",
      subtitle = "Dashed = 1:1; dotted = ±20 kt bounds",
      color = "Island", shape = "Island"
    ) +
    ggtheme

  ggsave("output/validation/wind_field_scatter.png", p_wf,
         width = 8, height = 6, dpi = 150)
  message("Saved: output/validation/wind_field_scatter.png")

  # Detailed per-storm distance profile for Irma at St. Martin
  if (!is.null(out$trackpoints$St_Martin)) {
    irma_tp <- out$trackpoints$St_Martin |>
      filter(SID == "2017242N16333") |>
      filter(is.finite(V_site_kt), is.finite(dist_km))

    if (nrow(irma_tp) > 0) {
      p_irma <- ggplot(irma_tp, aes(x = dist_km, y = V_site_kt)) +
        geom_line(color = "steelblue", linewidth = 0.8) +
        geom_point(color = "steelblue", size = 1.5) +
        geom_hline(yintercept = c(34, 64), linetype = "dashed",
                   color = c("orange", "red")) +
        geom_hline(yintercept = 155, linetype = "solid", color = "red", alpha = 0.5) +
        annotate("text", x = max(irma_tp$dist_km) * 0.7, y = 155,
                 label = "NHC best-track Vmax = 155 kt", size = 3, color = "red") +
        labs(
          x = "Distance from St. Martin (km)",
          y = "Model site wind (kt)",
          title = "Hurricane Irma (2017): Wind Profile at St. Martin",
          subtitle = "V_site_kt vs distance at each 6-hourly track point"
        ) +
        ggtheme

      ggsave("output/validation/irma_stmartin_profile.png", p_irma,
             width = 8, height = 5, dpi = 150)
      message("Saved: output/validation/irma_stmartin_profile.png")
    }
  }
}


# =============================================================================
# Print final results tables
# =============================================================================

message("\n\n>>> FINAL VALIDATION TABLES")

if (!is.null(val$hindcast$comparison)) {
  message("\n--- Hindcast Return Levels ---")
  print(knitr::kable(val$hindcast$comparison |>
                       select(island, return_period, obs_full_rl, sim_median,
                              sim_lo_90, sim_hi_90, bias_pct, obs_in_90ci) |>
                       mutate(across(where(is.numeric), ~round(., 1))),
                     format = "pipe"))
}

if (!is.null(val$rate_check)) {
  message("\n--- Rate Comparison ---")
  print(knitr::kable(val$rate_check |>
                       select(island, severity, lambda_model, lambda_ref, ratio, flag) |>
                       mutate(across(where(is.numeric), ~round(., 3))),
                     format = "pipe"))
}

if (!is.null(val$wind_field)) {
  message("\n--- Wind Field Spot-Checks ---")
  print(knitr::kable(val$wind_field |>
                       select(storm_name, island, obs_1min_equiv_kt,
                              model_V_site_kt, bias_kt, bias_pct, min_dist_km) |>
                       mutate(across(where(is.numeric), ~round(., 1))),
                     format = "pipe"))
}

message("\n--- Summary ---")
print(knitr::kable(val$summary, format = "pipe"))

message("\n>>> Validation complete.")
