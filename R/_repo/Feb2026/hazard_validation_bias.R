### BIAS DECOMPOSITION DIAGNOSTIC
### Run AFTER run_validation.R (requires `val` and `out` objects in memory)
###
### Decomposes the systematic overestimation into:
###   (a) Frequency bias   — too many events per year?
###   (b) Intensity bias   — events too strong?
###   (c) Zero-year bias   — too few quiet years?

library(dplyr)
library(tibble)
library(ggplot2)

# =============================================================================
# 1) Build comparison table across all islands
# =============================================================================

diag_rows <- list()

for (isl in names(val$hindcast$per_island)) {
  hc <- val$hindcast$per_island[[isl]]
  if (is.null(hc)) next
  
  obs_am  <- hc$obs_annual_max$V_max_kt          # full observed annual max vector
  sim_am  <- hc$sim_annual_max                     # 5000 simulated annual maxima
  
  # --- (a) Fraction of zero-event years ---
  obs_p0 <- mean(obs_am <= 0)
  sim_p0 <- mean(sim_am <= 0)
  
  # --- (b) Mean annual max CONDITIONAL on occurrence (nonzero years) ---
  obs_nz <- obs_am[obs_am > 0]
  sim_nz <- sim_am[sim_am > 0]
  
  obs_mean_nz <- if (length(obs_nz) > 0) mean(obs_nz) else NA
  sim_mean_nz <- if (length(sim_nz) > 0) mean(sim_nz) else NA
  
  # --- (c) Overall mean annual max (including zeros) ---
  obs_mean_all <- mean(obs_am)
  sim_mean_all <- mean(sim_am)
  
  # --- (d) Tail metrics: 90th and 95th percentile of annual max ---
  obs_p90 <- quantile(obs_am, 0.90)
  sim_p90 <- quantile(sim_am, 0.90)
  obs_p95 <- quantile(obs_am, 0.95)
  sim_p95 <- quantile(sim_am, 0.95)
  
  # --- (e) Frequency: mean total events/year from lambda table ---
  lt <- hc$train_params$lambda_table
  # lambda_table has columns: severity, lambda (mean rate)
  model_lambda_total <- sum(lt$lambda, na.rm = TRUE)
  
  # Observed: count event-years (V > 0) / total years as simple rate proxy
  n_obs_years <- length(obs_am)
  obs_event_rate <- sum(obs_am > 0) / n_obs_years
  sim_event_rate <- sum(sim_am > 0) / length(sim_am)
  
  # --- (f) Standard deviation (spread of annual max) ---
  obs_sd <- sd(obs_am)
  sim_sd <- sd(sim_am)
  
  diag_rows[[isl]] <- tibble(
    island           = isl,
    n_obs_years      = n_obs_years,
    # Zero-year fraction
    obs_p_zero       = round(obs_p0, 3),
    sim_p_zero       = round(sim_p0, 3),
    p_zero_diff      = round(sim_p0 - obs_p0, 3),
    # Event rate (fraction of years with any event)
    obs_event_rate   = round(obs_event_rate, 3),
    sim_event_rate   = round(sim_event_rate, 3),
    # Model fitted total lambda
    model_lambda     = round(model_lambda_total, 3),
    # Intensity: mean conditional on occurrence
    obs_mean_nz_kt   = round(obs_mean_nz, 1),
    sim_mean_nz_kt   = round(sim_mean_nz, 1),
    intensity_bias_kt = round(sim_mean_nz - obs_mean_nz, 1),
    intensity_bias_pct = round(100 * (sim_mean_nz - obs_mean_nz) / obs_mean_nz, 1),
    # Overall mean
    obs_mean_kt      = round(obs_mean_all, 1),
    sim_mean_kt      = round(sim_mean_all, 1),
    overall_bias_kt  = round(sim_mean_all - obs_mean_all, 1),
    overall_bias_pct = round(100 * (sim_mean_all - obs_mean_all) / obs_mean_all, 1),
    # Spread
    obs_sd_kt        = round(obs_sd, 1),
    sim_sd_kt        = round(sim_sd, 1),
    # Tail
    obs_p90_kt       = round(obs_p90, 1),
    sim_p90_kt       = round(sim_p90, 1),
    obs_p95_kt       = round(obs_p95, 1),
    sim_p95_kt       = round(sim_p95, 1),
    tail_bias_p95_pct = round(100 * (sim_p95 - obs_p95) / obs_p95, 1),
    # GEV shape (positive = heavier tail)
    gev_xi           = round(hc$diagnostics$gev_xi, 3),
    k_hat            = round(hc$train_params$k_hat, 2)
  )
}

diag_tbl <- bind_rows(diag_rows)

# =============================================================================
# 2) Print summary
# =============================================================================

message("\n", paste(rep("=", 72), collapse = ""))
message("  BIAS DECOMPOSITION DIAGNOSTIC")
message(paste(rep("=", 72), collapse = ""))

message("\n--- (a) Zero-year fraction: obs vs sim ---")
message("  (sim_p0 < obs_p0 means model produces fewer quiet years → inflates annual max)")
print(knitr::kable(
  diag_tbl |> select(island, n_obs_years, obs_p_zero, sim_p_zero, p_zero_diff),
  format = "pipe"
))

message("\n--- (b) Mean intensity conditional on occurrence ---")
message("  (positive bias = model events are too strong)")
print(knitr::kable(
  diag_tbl |> select(island, obs_mean_nz_kt, sim_mean_nz_kt,
                     intensity_bias_kt, intensity_bias_pct),
  format = "pipe"
))

message("\n--- (c) Overall mean annual max ---")
message("  (combines frequency and intensity effects)")
print(knitr::kable(
  diag_tbl |> select(island, obs_mean_kt, sim_mean_kt,
                     overall_bias_kt, overall_bias_pct),
  format = "pipe"
))

message("\n--- (d) Tail behaviour (P90, P95) ---")
message("  (where return-level overestimation lives)")
print(knitr::kable(
  diag_tbl |> select(island, obs_p90_kt, sim_p90_kt, obs_p95_kt, sim_p95_kt,
                     tail_bias_p95_pct, gev_xi),
  format = "pipe"
))

message("\n--- (e) Variability & overdispersion ---")
print(knitr::kable(
  diag_tbl |> select(island, obs_sd_kt, sim_sd_kt, model_lambda, k_hat),
  format = "pipe"
))

# =============================================================================
# 3) Dominant bias source summary
# =============================================================================

message("\n--- DOMINANT BIAS SOURCE PER ISLAND ---")
for (i in seq_len(nrow(diag_tbl))) {
  r <- diag_tbl[i, ]
  
  freq_signal   <- abs(r$p_zero_diff) > 0.05
  intens_signal <- abs(r$intensity_bias_pct) > 5
  tail_signal   <- abs(r$tail_bias_p95_pct) > 10
  
  sources <- character(0)
  if (freq_signal)
    sources <- c(sources, sprintf("frequency (p0 diff = %+.3f)", r$p_zero_diff))
  if (intens_signal)
    sources <- c(sources, sprintf("intensity (%+.1f%% | %+.1f kt)", 
                                  r$intensity_bias_pct, r$intensity_bias_kt))
  if (tail_signal)
    sources <- c(sources, sprintf("upper tail (P95 %+.1f%%)", r$tail_bias_p95_pct))
  
  if (length(sources) == 0) sources <- "within tolerance"
  
  message(sprintf("  %-12s: %s", r$island, paste(sources, collapse = " + ")))
}

# =============================================================================
# 4) Diagnostic plots
# =============================================================================

dir.create("output/validation", recursive = TRUE, showWarnings = FALSE)

ggtheme <- theme_light(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# --- QQ-style: simulated vs observed quantiles of annual max ---
qq_list <- list()
for (isl in names(val$hindcast$per_island)) {
  hc <- val$hindcast$per_island[[isl]]
  if (is.null(hc)) next
  
  obs_q <- quantile(hc$obs_annual_max$V_max_kt, probs = seq(0.05, 0.95, 0.05))
  sim_q <- quantile(hc$sim_annual_max, probs = seq(0.05, 0.95, 0.05))
  
  qq_list[[isl]] <- tibble(
    island = isl,
    prob = seq(0.05, 0.95, 0.05),
    obs_q = obs_q,
    sim_q = sim_q
  )
}

qq_df <- bind_rows(qq_list)

p_qq <- ggplot(qq_df, aes(x = obs_q, y = sim_q)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ island, scales = "free", ncol = 3) +
  labs(
    x = "Observed quantile (kt)",
    y = "Simulated quantile (kt)",
    title = "QQ Plot: Simulated vs Observed Annual Maximum Wind",
    subtitle = "Points above dashed line = model overestimation at that quantile"
  ) +
  ggtheme

ggsave("output/validation/bias_qq_annual_max.png", p_qq,
       width = 12, height = 7, dpi = 150)
message("\nSaved: output/validation/bias_qq_annual_max.png")


# --- Bar chart: bias decomposition ---
decomp_df <- diag_tbl |>
  select(island, p_zero_diff, intensity_bias_pct, tail_bias_p95_pct) |>
  tidyr::pivot_longer(-island, names_to = "source", values_to = "value") |>
  mutate(source = case_when(
    source == "p_zero_diff"        ~ "Frequency\n(Δ p_zero)",
    source == "intensity_bias_pct" ~ "Intensity\n(% bias | V>0)",
    source == "tail_bias_p95_pct"  ~ "Upper tail\n(% bias P95)"
  ))

p_decomp <- ggplot(decomp_df, aes(x = source, y = value, fill = island)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(
    x = NULL, y = "Bias metric",
    title = "Bias Decomposition by Source",
    subtitle = "Frequency: Δp_zero (negative = fewer quiet years); Intensity & Tail: % overestimation",
    fill = "Island"
  ) +
  ggtheme

ggsave("output/validation/bias_decomposition.png", p_decomp,
       width = 10, height = 5, dpi = 150)
message("Saved: output/validation/bias_decomposition.png")


# --- CDF comparison per island ---
cdf_list <- list()
for (isl in names(val$hindcast$per_island)) {
  hc <- val$hindcast$per_island[[isl]]
  if (is.null(hc)) next
  
  obs_sorted <- sort(hc$obs_annual_max$V_max_kt)
  n_obs <- length(obs_sorted)
  
  sim_sorted <- sort(hc$sim_annual_max)
  n_sim <- length(sim_sorted)
  
  cdf_list[[isl]] <- bind_rows(
    tibble(island = isl, source = "Observed",
           V_kt = obs_sorted,
           cdf = seq_len(n_obs) / n_obs),
    tibble(island = isl, source = "Simulated",
           V_kt = sim_sorted,
           cdf = seq_len(n_sim) / n_sim)
  )
}

cdf_df <- bind_rows(cdf_list)

p_cdf <- ggplot(cdf_df, aes(x = V_kt, y = cdf, color = source)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ island, scales = "free_x", ncol = 3) +
  scale_color_manual(values = c(Observed = "red", Simulated = "steelblue")) +
  labs(
    x = "Annual maximum wind (kt)",
    y = "CDF",
    title = "CDF Comparison: Simulated vs Observed Annual Maxima",
    subtitle = "Simulated curve right of observed = overestimation",
    color = NULL
  ) +
  ggtheme

ggsave("output/validation/bias_cdf_comparison.png", p_cdf,
       width = 12, height = 7, dpi = 150)
message("Saved: output/validation/bias_cdf_comparison.png")


message("\n>>> Bias decomposition complete.")
message(">>> Add 'source(\"diagnose_bias.R\")' at the end of run_validation.R")