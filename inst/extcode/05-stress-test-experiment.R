# =============================================================================
# Stress-test experiment: natural variability selection + climate change impact
#
# Purpose:
#   Part 1 — Select 5 representative high-impact years from 2,000 synthetic
#             baseline years. Candidates must include at least one event at or
#             above the IRMA 2017 site-level impact (peak wind >= 80 kt). Years
#             are rated on a two-metric composite score (r_peak, r_damage) with
#             configurable weights, then sampled at Q50/Q80/Q90/Q95/Q99.
#
#   Part 2 — Replay those same 5 year indices under KNMI'23 climate scenarios
#             (knmi_Ld, knmi_Hd) using the same seed (common random numbers).
#             Two experiment modes are available via FREEZE_FREQUENCY:
#               FALSE — full response: both frequency and intensity change.
#               TRUE  — intensity-only: frequency frozen at baseline, only
#                       event intensities are perturbed by warming.
#
# Key config switches (all in the Config block below):
#   TARGET_YEAR      — 2050 or 2100
#   FREEZE_FREQUENCY — FALSE (full) or TRUE (intensity-only)
#   w_peak / w_damage — relative weights for the composite stress rating
#
# Outputs:
#   selected_ids    — integer[5]: replicable year indices (seed = SEED)
#   baseline_stress — named list[location] of daily tibbles for 5 baseline years
#   cc_stress       — named list[scenario][location] of daily tibbles for 5 years
#
# Data path:
#   Uses the packaged demo IBTrACS subset by default.
#   For production, replace `data_path` with the full IBTrACS NA CSV path.
# =============================================================================

library(dplyr)
library(ipdcstorm)

# =============================================================================
# Config
# =============================================================================

SEED           <- 42L
N_SIM          <- 2000L
IRMA_WIND_KT   <- 80     # site-level peak wind threshold for IRMA-like impact (kt)
QUANTILE_PROBS <- c(0.50, 0.80, 0.90, 0.95, 0.99)
KNMI_SCENARIOS <- c("knmi_Ld", "knmi_Hd")

# Future target year: set to 2050 (near-term) or 2100 (end-of-century).
# Controls the SST delta used for all KNMI climate scenarios in Part 2.
# KNMI'23 provides distinct delta_sst values for each horizon:
#   knmi_Ld / knmi_Hd at 2050 → moderate forcing
#   knmi_Ld / knmi_Hd at 2100 → full end-of-century forcing (knmi_Hd reaches ~2.6 °C)
TARGET_YEAR <- 2100   # switch to 2100 for end-of-century analysis

# Controlled experiment mode (TRUE / FALSE).
#
# FALSE (default): full climate response — KNMI scenarios change both storm
#   frequency (via climate-adjusted Poisson/NB rates in out_cc$sim) and
#   individual event intensity (via perturb_event()). Years with different
#   storm counts from the baseline are expected and represent the genuine
#   frequency signal of warming.
#
# TRUE: intensity-only response — storm frequency is frozen at the baseline
#   by substituting out_base$sim into the KNMI out object before daily
#   generation. Only perturb_event() is active, so every scenario year
#   contains the same events as the corresponding baseline year, each
#   shifted upward in intensity by the scenario's delta_sst. Use this mode
#   to isolate the pure intensity signal and produce directly comparable
#   traces across baseline and scenarios.
FREEZE_FREQUENCY <- FALSE

# Target locations (Dutch Caribbean + Puerto Rico + Miami)
targets <- tibble::tribble(
  ~name,         ~lat,     ~lon,
  "St_Martin",   18.0708, -63.0501,
  "Saba",        17.6350, -63.2300,
  "Statia",      17.4890, -62.9740,
  "Puerto_Rico", 18.4655, -66.1057,
  "Miami",       25.7617, -80.1918
)

# Data: packaged demo subset — replace with full IBTrACS CSV for production runs
data_path <- system.file("extdata", "ibtracs_demo.csv", package = "ipdcstorm")
# data_path <- "/path/to/ibtracs.NA.list.v04r01.csv"

# =============================================================================
# Part 1: Stationary baseline — 2,000 synthetic years
# =============================================================================

message("Part 1: Stationary baseline (", N_SIM, " years, seed = ", SEED, ")")

cfg_base <- make_hazard_cfg(
  data_path        = data_path,
  simulation_years = N_SIM,
  climate          = make_climate_cfg(scenario = "stationary")
)

out_base <- run_hazard_model(cfg_base, targets = targets, seed = SEED, verbose = FALSE)

# --- Generate daily hazard for all 2,000 years ---
message("  Generating ", N_SIM, " years of daily hazard ...")

daily_base_list <- generate_daily_hazard_impact_spatial(
  out         = out_base,
  location    = targets$name,
  sim_years   = seq_len(N_SIM),
  year0       = 2000L,
  gust_factor = 1.0,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "baseline",
  seed        = SEED
)



# Combine into a single tibble (location column already present in each element)
daily_base <- dplyr::bind_rows(daily_base_list)

# --- Compute annual metrics aggregated across all target locations ---

annual_loc <- daily_base |>
  group_by(sim_year, location) |>
  summarise(
    peak_wind_kt = max(wind_kt, na.rm = TRUE),
    n_hur_days   = sum(wind_kt >= 64, na.rm = TRUE),
    n_storm_days = sum(wind_kt >= 34, na.rm = TRUE),
    total_damage = sum(damage_rate, na.rm = TRUE),
    .groups = "drop"
  )

annual_all <- annual_loc |>
  group_by(sim_year) |>
  summarise(
    peak_wind_max    = max(peak_wind_kt),
    hur_days_total   = sum(n_hur_days),
    storm_days_total = sum(n_storm_days),
    damage_total     = sum(total_damage),
    .groups = "drop"
  )

# Attach annual storm count from simulation table for compounding metric.
# Use the first location — counts are shared across a single sim draw.
sim_counts <- out_base$sim |>
  filter(location == targets$name[1]) |>
  select(sim_year, n_total)

year_metrics <- annual_all |>
  left_join(sim_counts, by = "sim_year")

# --- Filter: require at least one IRMA-level event in the year ---
irma_candidates <- year_metrics |>
  filter(peak_wind_max >= IRMA_WIND_KT)

message("  IRMA-like candidates (peak wind >= ", IRMA_WIND_KT, " kt): ",
        nrow(irma_candidates), " of ", N_SIM, " years")

# --- Rank candidates on a composite stress rating ---
#
# Two metrics, each rank-normalised to [0, 1] within the candidate pool:
#
#   r_peak   : max(peak_wind_kt) across all sites — worst instantaneous
#              intensity; governs single-event structural damage threshold
#              exceedance (loss scales ~wind^3).
#
#   r_damage : sum(damage_rate) across all days and sites — integrated
#              portfolio impact; embeds intensity, duration, and multi-event
#              accumulation in one physically meaningful number.
#
# Weights (w_peak + w_damage must sum to 1):
#   Equal weighting by default; increase w_peak to favour severity-driven
#   selection, increase w_damage to favour accumulation-driven selection.

w_peak   <- 0.5
w_damage <- 0.5

n_cand <- nrow(irma_candidates)

irma_rated <- irma_candidates |>
  mutate(
    r_peak   = rank(peak_wind_max, ties.method = "average") / n_cand,
    r_damage = rank(damage_total,  ties.method = "average") / n_cand,
    rating   = w_peak * r_peak + w_damage * r_damage
  ) |>
  arrange(rating)

# --- Select 5 years at target quantile levels ---
sel_pos    <- pmax(1L, ceiling(QUANTILE_PROBS * n_cand))
selected_ids <- irma_rated$sim_year[sel_pos]

message("  Selected stress-test year IDs:")
for (i in seq_along(QUANTILE_PROBS)) {
  yr  <- selected_ids[i]
  row <- irma_rated[irma_rated$sim_year == yr, ]
  message(sprintf(
    "    Q%02d: year %4d | peak=%3.0f kt  hur_days=%2d  n_storms=%d  damage=%.4f  rating=%.3f",
    round(QUANTILE_PROBS[i] * 100), yr,
    row$peak_wind_max, row$hur_days_total,
    row$n_total, row$damage_total, row$rating
  ))
}

# --- Extract baseline daily series for the 5 selected years ---
# Filter from the already-generated data to avoid redundant computation.
baseline_stress <- lapply(daily_base_list, function(df) {
  df[df$sim_year %in% selected_ids, ]
})

# --- Identify focal SID for each selected year (basin-level peak HUR event) ---
#
# The event_id in the daily output is formatted as "{SID}_y{cal_year}_{counter}".
# Stripping the suffix recovers the original IBTrACS SID.
#
# Basin-level: pool all target locations, pick the SID driving the highest
# single-day wind reading across all sites in that year. This SID will be
# pinned into every scenario run so the focal event always appears, with its
# intensity and duration perturbed by the scenario's delta_sst.

message("\nExtracting basin-level focal SIDs for selected years ...")

focal_sids <- stats::setNames(vector("list", length(selected_ids)),
                               as.character(selected_ids))

for (yr in selected_ids) {
  yr_hur <- dplyr::bind_rows(lapply(baseline_stress, function(df) {
    df[df$sim_year == yr & !is.na(df$event_class) & df$event_class == "HUR", ]
  }))
  if (nrow(yr_hur) == 0L) {
    warning("No HUR-class days found in baseline for year ", yr, ". Focal pin skipped.")
    focal_sids[[as.character(yr)]] <- NA_character_
    next
  }
  peak_row <- yr_hur[which.max(yr_hur$wind_kt), ]
  focal_sids[[as.character(yr)]] <- sub("_y[0-9]+_[0-9]+$", "", peak_row$event_id[1L])
}

message("  Focal SIDs:")
for (i in seq_along(selected_ids)) {
  yr  <- selected_ids[i]
  sid <- focal_sids[[as.character(yr)]]
  message(sprintf("    Q%02d (year %4d): %s",
                  round(QUANTILE_PROBS[i] * 100), yr,
                  if (is.na(sid)) "<none>" else sid))
}

# =============================================================================
# Part 2: KNMI'23 climate scenarios — same seed, same year indices
#
# Experiment mode: focal event pinned across all scenarios.
#
#   The dominant HUR event from each baseline stress-test year is guaranteed
#   to appear in both KNMI scenario runs (same SID → same calendar timing).
#   perturb_event() shifts its intensity and duration by the scenario delta_sst.
#   All other storm slots are drawn freely from the climate-adjusted pool, so
#   compound-event frequency varies naturally between scenarios.
#
# FREEZE_FREQUENCY is retained for reference but is no longer the primary
# mechanism — set it TRUE only if you want to also lock non-focal storm counts.
# =============================================================================

message("\nPart 2: KNMI'23 scenarios at target year ", TARGET_YEAR)
message("  Scenarios       : ", paste(KNMI_SCENARIOS, collapse = ", "))
message("  Freeze frequency: ", FREEZE_FREQUENCY)
message("  Year indices    : ", paste(selected_ids, collapse = ", "))

cc_stress <- vector("list", length(KNMI_SCENARIOS))
names(cc_stress) <- KNMI_SCENARIOS

for (scen in KNMI_SCENARIOS) {

  message("  [", scen, "] Running model ...")

  cfg_cc <- make_hazard_cfg(
    data_path        = data_path,
    simulation_years = N_SIM,
    climate          = make_climate_cfg(scenario = scen,
                                        target_year = TARGET_YEAR,
                                        perturb = TRUE)
  )

  # Same seed → common random numbers; differences isolate climate forcing only.
  out_cc <- run_hazard_model(cfg_cc, targets = targets, seed = SEED, verbose = FALSE)

  # Controlled mode: substitute baseline sim table so storm counts per year
  # are identical to the baseline. Only perturb_event() (intensity) differs.
  if (FREEZE_FREQUENCY) {
    out_cc$sim    <- out_base$sim
    out_cc$config <- out_cc$config   # leave climate config intact for perturb
  }

  # Generate daily series for ALL N_SIM years, then filter to selected_ids.
  #
  # CRN requirement: both baseline and scenario runs must iterate through the
  # same sim_years vector (1:N_SIM) so that the RNG is consumed at the same
  # rate for every year index. Passing only selected_ids here would put the
  # scenario run's first iteration at the initial RNG state, while the
  # baseline reached those same year indices after consuming hundreds of prior
  # draws — resulting in completely different event sampling that has nothing
  # to do with climate forcing.
  cc_daily_all <- generate_daily_hazard_impact_spatial(
    out         = out_cc,
    location    = targets$name,
    sim_years   = seq_len(N_SIM),
    year0       = 2000L,
    gust_factor = 1.0,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = scen,
    seed        = SEED,
    pinned_sids = focal_sids
  )

  # Retain only the 5 selected stress-test years for downstream use
  cc_stress[[scen]] <- lapply(cc_daily_all, function(df) {
    df[df$sim_year %in% selected_ids, ]
  })
}

# =============================================================================
# Summary
# =============================================================================

message("\n--- Stress-test experiment complete ---")
message("Seed              : ", SEED)
message("Simulation years  : ", N_SIM)
message("Target year       : ", TARGET_YEAR)
message("IRMA threshold    : ", IRMA_WIND_KT, " kt")
message("Candidate years   : ", n_cand, " / ", N_SIM)
message("Rating metrics    : r_peak (w=", w_peak, ") + r_damage (w=", w_damage, ")")
message("Selected year IDs : ", paste(selected_ids, collapse = ", "),
        "  [", paste0("Q", round(QUANTILE_PROBS * 100), collapse = " / "), "]")
message("CC scenarios      : ", paste(KNMI_SCENARIOS, collapse = ", "))
message("Freeze frequency  : ", FREEZE_FREQUENCY)
message("")
message("Objects in workspace:")
message("  selected_ids      — integer[5]: year indices for all scenario comparisons")
message("  focal_sids        — list[sim_year]: basin-level focal SID pinned across scenarios")
message("  baseline_stress   — list[location]: daily tibbles for 5 baseline years")
message("  cc_stress         — list[scenario][location]: daily tibbles for 5 CC years")
