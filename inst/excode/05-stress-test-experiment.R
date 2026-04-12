# =============================================================================
# Stress-test experiment: natural variability selection + climate change impact
#
# Purpose:
#   Part 1 — Select 5 representative high-impact years from 2,000 synthetic
#             baseline years. Candidates must include at least one event at or
#             above the IRMA 2017 site-level impact (peak wind >= 80 kt). Years
#             are rated on severity, hurricane duration, compounding storm
#             activity, and total damage, then sampled at Q50/Q80/Q90/Q95/Q99.
#
#   Part 2 — Replay those same 5 year indices under four KNMI'23 climate
#             scenarios using the same seed (common random numbers), so that
#             any difference in outcomes is attributable to climate forcing
#             rather than sampling variability.
#
# Outputs:
#   selected_ids      — integer[5]: replicable year indices (seed = SEED)
#   baseline_stress   — named list[location] of daily tibbles for 5 baseline years
#   cc_stress         — named list[scenario][location] of daily tibbles for 5 years
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
TARGET_YEAR    <- 2050
IRMA_WIND_KT   <- 80     # site-level peak wind threshold for IRMA-like impact (kt)
QUANTILE_PROBS <- c(0.50, 0.80, 0.90, 0.95, 0.99)
KNMI_SCENARIOS <- c("knmi_Ld", "knmi_Hd")

# Target locations (Dutch Caribbean)
targets <- tibble::tribble(
  ~name,       ~lat,     ~lon,
  "St_Martin", 18.0708, -63.0501,
  "Saba",      17.6350, -63.2300,
  "Statia",    17.4890, -62.9740
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
# Four metrics (all higher = more stressful):
#   r_peak     : site-level peak wind (severity)
#   r_duration : total hurricane-intensity days across sites (persistence)
#   r_compound : number of storms in the year (compounding exposure)
#   r_damage   : cumulative damage rate across sites and year
# Each metric is rank-normalised within the candidate pool [0, 1].
# Equal weights; final rating averaged across four components.

n_cand <- nrow(irma_candidates)

irma_rated <- irma_candidates |>
  mutate(
    r_peak     = rank(peak_wind_max,  ties.method = "average") / n_cand,
    r_duration = rank(hur_days_total, ties.method = "average") / n_cand,
    r_compound = rank(n_total,        ties.method = "average") / n_cand,
    r_damage   = rank(damage_total,   ties.method = "average") / n_cand,
    rating     = (r_peak + r_duration + r_compound + r_damage) / 4
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

# =============================================================================
# Part 2: KNMI'23 climate scenarios — same seed, same year indices
# =============================================================================

message("\nPart 2: KNMI'23 scenarios at target year ", TARGET_YEAR)
message("  Scenarios: ", paste(KNMI_SCENARIOS, collapse = ", "))
message("  Replaying year indices: ", paste(selected_ids, collapse = ", "))

cc_stress <- vector("list", length(KNMI_SCENARIOS))
names(cc_stress) <- KNMI_SCENARIOS

for (scen in KNMI_SCENARIOS) {

  message("  [", scen, "] Running model ...")

  cfg_cc <- make_hazard_cfg(
    data_path        = data_path,
    simulation_years = N_SIM,
    climate          = make_climate_cfg(scenario = scen, target_year = TARGET_YEAR)
  )

  # Same seed → common random numbers; differences isolate climate forcing only.
  out_cc <- run_hazard_model(cfg_cc, targets = targets, seed = SEED, verbose = FALSE)

  # Generate daily series for ALL N_SIM years, then filter to selected_ids.
  #
  # CRN requirement: both baseline and scenario runs must iterate through the
  # same sim_years vector (1:N_SIM) so that the RNG is consumed at the same
  # rate for every year index. Passing only selected_ids here would put the
  # scenario run's first iteration at the initial RNG state, while the
  # baseline reached those same year indices after consuming hundreds of prior
  # draws — resulting in completely different event sampling that has nothing
  # to do with climate forcing.
  message("  [", scen, "] Generating ", N_SIM, " years of daily hazard ...")
  cc_daily_all <- generate_daily_hazard_impact_spatial(
    out         = out_cc,
    location    = targets$name,
    sim_years   = seq_len(N_SIM),
    year0       = 2000L,
    gust_factor = 1.0,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = scen,
    seed        = SEED
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
message("IRMA threshold    : ", IRMA_WIND_KT, " kt")
message("Candidate years   : ", n_cand, " / ", N_SIM)
message("Selected year IDs : ", paste(selected_ids, collapse = ", "),
        "  [", paste0("Q", round(QUANTILE_PROBS * 100), collapse = " / "), "]")
message("CC scenarios      : ", paste(KNMI_SCENARIOS, collapse = ", "))
message("")
message("Objects in workspace:")
message("  selected_ids      — integer[5]: year indices for all scenario comparisons")
message("  baseline_stress   — list[location]: daily tibbles for 5 baseline years")
message("  cc_stress         — list[scenario][location]: daily tibbles for 5 CC years")
