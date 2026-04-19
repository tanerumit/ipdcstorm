# =============================================================================
# Model validation pipeline
#
# Purpose:
#   Assess the skill of the stationary baseline hazard model through a
#   three-tier validation framework:
#
#   Tier 1 (Hindcast) — Hold out the last HOLDOUT_YEARS of the historical
#     record.  Re-fit on the training period, simulate N_SIM synthetic years,
#     and compare simulated return levels against observed return levels at
#     the specified return periods.  A well-calibrated model places simulated
#     return levels within the observed confidence interval at most sites.
#
#   Tier 2 (Rate check) — Compare calibrated lambda values (annual TS and HUR
#     rates per location) against published HURDAT2/IBTrACS climatologies.
#     Flags rates that deviate substantially from the reference range.
#
#   Tier 3 (Wind field) — Spot-check the Holland wind-profile model against
#     historical station observations for individual storm passages.  Reports
#     site-level wind bias.
#
# Key config switches (all in the Config block below):
#   SEED           — RNG seed; must match the baseline run for reproducibility
#   HOLDOUT_YEARS  — number of years withheld from training for hindcast
#   N_SIM          — synthetic years for hindcast comparison (NULL = inherit)
#   RETURN_PERIODS — return periods (yr) to evaluate in Tier 1
#   CONF_LEVEL     — width of the confidence interval around observed RL
#   SAVE_OUTPUTS   — write plots and tables to OUT_DIR when TRUE
#
# Outputs:
#   out_base      — run_hazard_model() output used as validation input
#   val_cfg       — validation_cfg object (validation parameters)
#   val_out       — run_validation_suite() output (summary, comparison,
#                   rate_check, wind_field, plots)
# =============================================================================

library(dplyr)
library(ipdcstorm)

# =============================================================================
# Config
# =============================================================================

SEED           <- 42L
N_SIM          <- 5000L
HIST_START_YEAR  <- 1970L
SEARCH_RADIUS_KM <- 800

# Holdout: number of years withheld from rate calibration for the hindcast.
# Common choices:
#   10 — standard; withholds the most recent decade (e.g. 2015–2024 if data
#        runs to 2024); large enough to include at least one or two major
#        events without over-shrinking the training period.
#   15 — stricter test; appropriate when the record extends to ≥ 1985.
#    5 — light holdout; useful for short records or high-return-period focus
#        where a longer training period is needed.
HOLDOUT_YEARS <- 10L

# Return periods (years) for hindcast comparison.
# Choose values where the observed record has meaningful empirical support
# (roughly ≤ record length / 3 for stability of empirical quantiles).
RETURN_PERIODS <- c(5, 10, 25, 50, 100)

# Confidence level for the GEV-based return-level confidence interval.
# The simulated RL is considered consistent if it falls within this CI.
# 0.95 is the standard validation threshold in this package.
CONF_LEVEL <- 0.95

# Output directory and file-saving flags.
# Set SAVE_OUTPUTS = TRUE to write PNG plots and CSV tables to OUT_DIR.
SAVE_OUTPUTS <- TRUE
OUT_DIR      <- "output/validation"

# Target locations (Dutch Caribbean islands)
targets <- tibble::tribble(
  ~name,         ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740
)

# Data: packaged demo subset — replace with full IBTrACS CSV for production runs
DATA_PATH <- system.file("extdata", "ibtracs_1970.csv", package = "ipdcstorm")


# =============================================================================
# Part 1: Run stationary baseline model
#
# Validation operates on the same model output as a standard analysis run.
# If you have already run 01-baseline-analysis.R in the same session, you can
# skip this block and pass `out_base` directly to run_validation_suite() below.
# =============================================================================

message("Part 1: Running stationary baseline model for validation")
message("  Seed            : ", SEED)
message("  Simulation years: ", N_SIM)
message("  Historical start: ", HIST_START_YEAR)
message("  Search radius   : ", SEARCH_RADIUS_KM, " km")

hazard_cfg <- make_hazard_cfg(
  data_path             = DATA_PATH,
  search_radius_km      = SEARCH_RADIUS_KM,
  historical_start_year = HIST_START_YEAR,
  simulation_years      = N_SIM,
  climate               = make_climate_cfg(scenario = "stationary")
)

out_base <- run_hazard_model(
  cfg     = hazard_cfg,
  targets = targets,
  seed    = SEED,
  verbose = FALSE
)

# =============================================================================
# Part 2: Configure and run validation function
#
# Runs all three tiers in sequence and (optionally) and write plots and tables
# to OUT_DIR.The returned object `val_out` contains structured results for each tier
# accessible as named list elements.
# =============================================================================

message("\nPart 2: Running three-tier validation suite")
message("  Holdout years  : ", HOLDOUT_YEARS)
message("  Return periods : ", paste(RETURN_PERIODS, collapse = ", "), " yr")
message("  Conf. level    : ", CONF_LEVEL * 100, "%")
message("  Save outputs   : ", SAVE_OUTPUTS)

val_cfg <- make_validation_cfg(
  holdout_years  = HOLDOUT_YEARS,
  n_sim          = N_SIM,
  return_periods = RETURN_PERIODS,
  conf_level     = CONF_LEVEL,
  seed           = SEED,
  out_dir        = OUT_DIR,
  save_plots     = SAVE_OUTPUTS,
  save_tables    = SAVE_OUTPUTS
)

val_out <- run_validation_suite(
  out = out_base,
  cfg = val_cfg
)

# =============================================================================
# Part 3: Inspect Tier 1 — Hindcast return levels
#
# val_out$comparison is a tibble with one row per (location, return_period)
# combination.  Key columns:
#   obs_full_rl — observed return level from the full record (kt)
#   sim_rl      — simulated return level from the N_SIM ensemble (kt)
#   obs_ci_lo   — lower bound of the CONF_LEVEL confidence interval (kt)
#   obs_ci_hi   — upper bound of the CONF_LEVEL confidence interval (kt)
#   obs_in_ci   — TRUE if sim_rl falls within [obs_ci_lo, obs_ci_hi]
# =============================================================================

message("\nPart 3: Tier 1 — Hindcast return-level comparison")

if (!is.null(val_out$comparison)) {
  hc_tbl <- val_out$comparison

  n_checks  <- nrow(hc_tbl)
  n_in_ci   <- sum(hc_tbl$obs_in_ci, na.rm = TRUE)
  pct_in_ci <- round(100 * n_in_ci / n_checks, 1)

  message(sprintf("  Pass rate: %d / %d checks within %d%% CI  (%s%%)",
                  n_in_ci, n_checks,
                  round(CONF_LEVEL * 100),
                  pct_in_ci))

  cat("\nHindcast return-level comparison:\n")
  print(
    hc_tbl |>
      mutate(across(where(is.numeric), ~ round(.x, 1))) |>
      arrange(location, return_period)
  )
} else {
  message("  Hindcast results not available (insufficient historical events?).")
}

# =============================================================================
# Part 4: Inspect Tier 2 — Rate sanity check
#
# val_out$rate_check compares calibrated site lambdas against reference
# climatology bounds.  The `status` column is "OK" or "FLAG".
# Flagged rates warrant inspection of the search radius or historical record.
# =============================================================================

message("\nPart 4: Tier 2 — Storm rate sanity check")

if (!is.null(val_out$rate_check) && nrow(val_out$rate_check) > 0) {
  cat("\nRate check results:\n")
  print(
    val_out$rate_check |>
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )

  n_flagged <- sum(val_out$rate_check$flag == "FLAG", na.rm = TRUE)
  if (n_flagged > 0) {
    message(sprintf("  WARNING: %d rate(s) flagged — review search_radius_km or historical_start_year.",
                    n_flagged))
  } else {
    message("  All rates within reference bounds.")
  }
} else {
  message("  Rate check results not available.")
}

# =============================================================================
# Part 5: Inspect Tier 3 — Wind-field spot-check
#
# val_out$wind_field compares Holland model site winds against historical
# station observations for individual storms.  Reports signed bias (kt) and
# RMSE per location.  A bias close to 0 kt indicates good wind-profile
# calibration; persistent high bias suggests the radius-wind relationship
# may need tuning for this basin sub-region.
# =============================================================================

message("\nPart 5: Tier 3 — Wind-field spot-check")

if (!is.null(val_out$wind_field) && nrow(val_out$wind_field) > 0) {
  cat("\nWind-field bias estimates:\n")
  print(
    val_out$wind_field |>
      mutate(across(where(is.numeric), ~ round(.x, 2)))
  )
} else {
  message("  Wind-field results not available (station observations required).")
}

# =============================================================================
# Part 6: Overall validation summary
#
# val_out$summary is a compact tibble with aggregate diagnostics across all
# tiers.  Use this for quick pass/fail assessment before running scenario
# analyses — a model that fails Tier 1 at more than one site should be
# recalibrated (adjust search_radius_km or historical_start_year) before
# drawing conclusions from climate projections.
# =============================================================================

message("\nPart 6: Overall validation summary")

if (!is.null(val_out$summary) && nrow(val_out$summary) > 0) {
  cat("\nValidation summary:\n")
  print(val_out$summary)
} else {
  message("  Summary not available.")
}

# =============================================================================
# Summary
# =============================================================================

message("\n--- Model validation complete ---")
message("Seed              : ", SEED)
message("Simulation years  : ", N_SIM)
message("Holdout years     : ", HOLDOUT_YEARS)
message("Return periods    : ", paste(RETURN_PERIODS, collapse = ", "), " yr")
message("Conf. level       : ", CONF_LEVEL * 100, "%")
message("Save outputs      : ", SAVE_OUTPUTS)
if (SAVE_OUTPUTS) {
  message("Output directory  : ", OUT_DIR)
  message("  Plots (ggplot)  : ", paste(names(val_out$plots), collapse = ", "))
}
message("")
message("Objects in workspace:")
message("  out_base  — list: run_hazard_model() output (events, fit, rates, sim, config)")
message("  val_cfg   — validation_cfg: validation parameters")
message("  val_out   — list: three-tier validation results")
message("              val_out$summary    — compact pass/fail summary across all tiers")
message("              val_out$comparison — Tier 1 hindcast return-level comparison")
message("              val_out$rate_check — Tier 2 storm rate sanity check")
message("              val_out$wind_field — Tier 3 wind-field bias estimates")
message("              val_out$plots      — named list of ggplot objects")
