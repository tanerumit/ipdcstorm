# =============================================================================
# Stress-test experiment: focal-anchored compound stress years under climate change
#
# Purpose
# -------
#   Goal: assess climate vulnerability to major hurricanes and explore the
#   impact of single and compounding events at Dutch Caribbean sites.
#
#   Design: every synthetic year is anchored with one historical major-
#   hurricane SID sampled (with replacement) from a fixed pool. Year-to-year
#   variability comes from (a) which pool SID was pinned, and (b) the
#   stochastic non-focal storm draws in that year. A guaranteed major event
#   plus random background activity yields realistic compound stress years
#   at a modest ensemble size, without requiring a 10^4-year random ensemble
#   to find rare Irma-class events by chance.
#
#   Part 1  -  Baseline. Run an anchored ensemble of N_BASE years, rank on
#              compound damage across Saba + Statia over a 60-day window
#              starting at focal-event onset, and select 5 quantile years
#              (Q50/Q80/Q90/Q95/Q99).
#
#   Part 2  -  KNMI'23 scenarios. Replay the 5 selected years under warming.
#              Each year's pinned focal SID is pinned again so the anchored
#              event carries scenario-scaled intensity, radius and duration
#              via perturb_event().
#
# Key config switches
# -------------------
#   PINNED_SID_POOL   - SIDs of historical major hurricanes to anchor into
#                       baseline years. Every SID must exist in the IBTrACS
#                       input or the pin will silently fail to land.
#   N_BASE            - anchored baseline ensemble size
#   TARGET_YEAR       - 2050 or 2100, controls KNMI delta_sst
#
# Outputs
# -------
#   selected_ids      - integer[5]: replicable baseline year indices (seed = SEED)
#   focal_sids        - list[sim_year]: pinned focal SID per selected year
#   baseline_stress   - named list[location] of daily tibbles for the 5 selected years
#   cc_stress         - named list[scenario][location] of replayed daily tibbles
# =============================================================================

library(dplyr)
library(ipdcstorm)

# =============================================================================
# Config
# =============================================================================

SEED             <- 123
N_BASE           <- 2500L
QUANTILE_PROBS   <- c(0.50, 0.75, 0.90, 0.95, 0.99, 0.995)
KNMI_SCENARIOS   <- c("knmi_Ld", "knmi_Hd")
OUTPUT_DIR       <- file.path("output", "raw")

# Compound-stress window (days after focal-event start). 90 days captures
# late-season pins (Lenny in November) against pre-focal Aug/Sep activity
# plus post-focal aftermath through January.
COMPOUND_WINDOW_DAYS <- 90L

# --- Stochastic pin jitter (internal variability of pinned focal events) ----
# The deterministic IBTrACS track is one realisation out of a distribution of
# similar storms. Per-year jitter perturbs the pinned event's timing,
# peak intensity and RMW, so each anchored year is a distinct realisation
# instead of a verbatim replay. Jitter applies BEFORE climate-scenario
# perturbation, and the SAME per-year jitter is used in baseline and KNMI
# scenarios so each year keeps its focal date across scenarios.
#
# Timing modes:
#   "uniform_season" (default) - resample focal DOY uniformly within the
#     Atlantic hurricane season [PIN_DOY_MIN, PIN_DOY_MAX]. The historical
#     DOY of each SID is discarded, so a given year's focal may be Aug, Sep,
#     Oct or Nov regardless of when the SID actually occurred.
#   "gaussian"       - small jitter around the historical DOY (std dev
#     PIN_DOY_SD). Preserves each SID's climatological date.
#
# Intensity/RMW jitter is always Gaussian with small SD (~5%).
PIN_TIMING_MODE  <- "uniform_season"
PIN_DOY_MIN      <- 152L     # June 1  (DOY, non-leap)
PIN_DOY_MAX      <- 304L     # October 31 (leave 90-day window before year end)
PIN_DOY_SD       <- 7        # used only when PIN_TIMING_MODE == "gaussian"
PIN_V_SCALE_SD   <- 0.05
PIN_R_SCALE_SD   <- 0.05

# --- State-dependent damage (HAZUS-HM-style progressive vulnerability) ------
# Damaged structures are more vulnerable to subsequent wind exposure. Each
# day's damage rate is amplified by cumulative season-to-date damage:
#   d_new(t) = d_raw(t) * (1 + DAMAGE_ALPHA * min(V(t-1), DAMAGE_V_CAP))
#
# Defaults:
#   DAMAGE_ALPHA = 3    (moderate amplification; lit range 2-5)
#   DAMAGE_V_CAP = 0.5  (cap at 50% cumulative damage state)
DAMAGE_ALPHA     <- 3
DAMAGE_V_CAP     <- 0.5

# Pool of historical major hurricanes used as focal anchors. Each baseline
# year is assigned one pool SID by random sampling with replacement; that
# SID is pinned into the year's HUR storm slots.
#
# Pool selection criteria (Cat-4+ at lifetime peak AND closest approach
# <= 300 km from the Saba/Statia mid-point), chosen for diversity in track
# geometry and intensity. All SIDs verified present in
# inst/extdata/ibtracs_1970.csv as of the most recent build.
#
#   IRMA    2017  155 kt   57 km   Cat-5 direct hit
#   MARIA   2017  150 kt  115 km   Cat-5 passing south of the islands
#   HUGO    1989  140 kt   85 km   classic Cat-4 reference
#   LUIS    1995  130 kt   83 km   devastating St Martin landfall
#   LENNY   1999  135 kt   49 km   unusual east-bound November Cat-4
#   GEORGES 1998  135 kt   46 km   closest Cat-4 approach in the record
PINNED_SID_POOL  <- c(
  "2017242N16333",   # Irma     (2017)
  "2017260N12310",   # Maria    (2017)
  "1989254N13340",   # Hugo     (1989)
  "1995241N11333",   # Luis     (1995)
  "1999318N17278",   # Lenny    (1999)
  "1998259N10335"    # Georges  (1998)
)

# Human-readable label for each pool SID (for plots and tables).
PINNED_SID_NAMES <- c(
  "2017242N16333" = "Irma 2017",
  "2017260N12310" = "Maria 2017",
  "1989254N13340" = "Hugo 1989",
  "1995241N11333" = "Luis 1995",
  "1999318N17278" = "Lenny 1999",
  "1998259N10335" = "Georges 1998"
)

# Future target year for KNMI'23 scenarios. 2100 gives full end-of-century
# forcing (knmi_Hd reaches ~2.6 deg C). Switch to 2050 for near-term.
TARGET_YEAR      <- 2100

# KNMI scenarios always apply the full climate response: storm frequency
# changes via climate-adjusted Poisson/NB rates in out_cc$sim AND individual
# event intensity/radius/duration change via perturb_event. The focal SID is
# still pinned in every replay year, so each year guarantees a scenario-
# scaled major event while non-focal storms reflect the scenario's frequency
# signal.

# Target locations (Dutch Caribbean + Puerto Rico + Miami)
targets <- tibble::tribble(
  ~name,         ~lat,     ~lon,
  "St_Martin",   18.0708, -63.0501,
  "Saba",        17.6350, -63.2300,
  "Statia",      17.4890, -62.9740,
  "Puerto_Rico", 18.4655, -66.1057,
  "Miami",       25.7617, -80.1918
)

# Use the repo-tracked IBTrACS file directly so edits to the source CSV are
# picked up without reinstalling the package. Falls back to the installed
# copy if the script is sourced from outside the repo.
data_path <- file.path("inst", "extdata", "ibtracs_1970.csv")
if (!file.exists(data_path)) {
  data_path <- system.file("extdata", "ibtracs_1970.csv", package = "ipdcstorm")
}

# =============================================================================
# Background wind configuration
# =============================================================================
#
# Literature basis:
#   - Puerto Rico / U.S. Virgin Islands Wind Energy Resource Atlas
#     (Elliott et al., 1990; DOE/PNL) reports Rayleigh-like trade-wind
#     distributions and annual mean near-surface winds around 4.4 m/s at
#     San Juan, with stronger exposed coastal and ridge sites.
#   - Caribbean low-level-jet climatology shows the eastern Caribbean and
#     Lesser Antilles sit in a persistently strong easterly trade-wind belt,
#     supporting slightly stronger background winds there than at San Juan.
#   - South Florida station climatology is weaker than the eastern Caribbean,
#     so Miami is assigned a lower mean background wind.
#
# These are representative stress-test values, not local station fits. For
# production work, refit Weibull marginals from local non-TC observations.
# The Weibull marginals themselves are encoded in
# ipdcstorm:::.make_stress_test_background_wind_cfg().

bg_cfg <- ipdcstorm:::.make_stress_test_background_wind_cfg()

# =============================================================================
# Part 1: Anchored baseline
# =============================================================================

message("Part 1: Anchored baseline (", N_BASE, " years, seed = ", SEED, ")")
message("  Focal pool: ", paste(PINNED_SID_POOL, collapse = ", "))

baseline_pin_map    <- ipdcstorm:::.build_pin_map(seq_len(N_BASE), PINNED_SID_POOL, SEED)
baseline_pin_jitter <- ipdcstorm:::.build_pin_jitter(
  sim_years  = seq_len(N_BASE),
  mode       = PIN_TIMING_MODE,
  doy_sd     = PIN_DOY_SD,
  doy_min    = PIN_DOY_MIN,
  doy_max    = PIN_DOY_MAX,
  v_scale_sd = PIN_V_SCALE_SD,
  r_scale_sd = PIN_R_SCALE_SD,
  seed       = SEED + 1L
)

cfg_base <- make_hazard_cfg(
  data_path        = data_path,
  simulation_years = N_BASE,
  climate          = make_climate_cfg(scenario = "stationary"),
  background_wind  = bg_cfg
)

out_base <- run_hazard_model(cfg_base, targets = targets, seed = SEED, verbose = FALSE)

message("  Generating ", N_BASE, " years of daily hazard with pinned focal events ...")

daily_base_list <- generate_daily_hazard_impact_spatial(
  out         = out_base,
  location    = targets$name,
  sim_years   = seq_len(N_BASE),
  year0       = 2000L,
  gust_factor = 1.0,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "baseline",
  seed        = SEED,
  pinned_sids = baseline_pin_map,
  pin_jitter  = baseline_pin_jitter
)

# State-dependent damage amplification: d(t) = d_raw(t) * (1 + alpha * V(t-1))
message("  Applying state-dependent damage amplification ",
        "(alpha = ", DAMAGE_ALPHA, ", cap = ", DAMAGE_V_CAP, ") ...")
daily_base_list <- ipdcstorm:::.apply_state_dependent_damage(
  daily_list = daily_base_list,
  alpha      = DAMAGE_ALPHA,
  cap        = DAMAGE_V_CAP
)

# --- Verify every year's pinned SID actually landed ---
missed_baseline <- ipdcstorm:::.verify_pins_landed(daily_base_list, baseline_pin_map)
if (nrow(missed_baseline) > 0L) {
  warning(
    nrow(missed_baseline), " of ", N_BASE,
    " baseline years did not contain the pinned focal SID. ",
    "This typically occurs when the HUR count draw is 0 for those years. ",
    "Inspect `missed_baseline` in the workspace for affected sim_years.",
    call. = FALSE
  )
  assign("missed_baseline", missed_baseline, envir = globalenv())
} else {
  message("  Pin check: all ", N_BASE, " baseline years carry their focal SID.")
}

baseline_paths <- save_daily_hazard_csvs(
  daily    = daily_base_list,
  scenario = "baseline",
  out_dir  = OUTPUT_DIR
)
message("  Wrote baseline CSVs to ", OUTPUT_DIR, " (", length(baseline_paths), " files)")

# --- Compute compound stress metrics from Saba + Statia only ---
compound_metrics <- ipdcstorm:::compute_compound_stress_year_metrics(
  daily       = daily_base_list,
  sim_years   = seq_len(N_BASE),
  location    = c("Saba", "Statia"),
  window_days = COMPOUND_WINDOW_DAYS
)

sim_counts <- out_base$sim |>
  filter(location == targets$name[1]) |>
  select(sim_year, n_total)

# Annual total damage at Saba + Statia — tie-breaker for the rating.
# Pinned focal events are deterministic, so years that share a pool SID and
# have no follow-on storms within the 60-day window all tie at the same
# compound_cum_damage. The annual sum captures any non-focal damage
# elsewhere in the year (e.g., an August storm in a November-Lenny year),
# which separates otherwise-degenerate ties.
annual_damage <- dplyr::bind_rows(
  daily_base_list[["Saba"]],
  daily_base_list[["Statia"]]
) |>
  dplyr::group_by(.data$sim_year) |>
  dplyr::summarise(
    annual_total_damage = sum(.data$damage_rate, na.rm = TRUE),
    .groups = "drop"
  )

year_metrics <- compound_metrics |>
  left_join(sim_counts, by = "sim_year") |>
  left_join(annual_damage, by = "sim_year")

# --- Rank anchored years: compound damage primary, annual damage secondary ---
#
# Every year is anchored by a pool SID, so the focal peak varies only across
# the (small) pool. The meaningful stress gradient is compound cumulative
# damage over the 60-day window at Saba + Statia, which embeds focal-event
# losses plus any compounding follow-on damage in that window.
#
# Tie-breaker: annual_total_damage at Saba + Statia. Within ties on the
# compound metric (typical for late-season pins where the 60-day window
# falls in a quiet period), this distinguishes years by total in-year stress.
#
# Focal peak is retained as a diagnostic, not a ranking axis.
n_cand <- nrow(year_metrics)

if (n_cand == 0L) {
  stop("Compound metrics table is empty; check daily output and location filter.",
       call. = FALSE)
}

compound_rated <- year_metrics |>
  mutate(
    r_damage = rank(compound_cum_damage, ties.method = "average") / n_cand
  ) |>
  arrange(compound_cum_damage, annual_total_damage) |>
  mutate(rating = seq_len(dplyr::n()) / n_cand)

sel_pos      <- pmax(1L, ceiling(QUANTILE_PROBS * n_cand))
selected_ids <- compound_rated$sim_year[sel_pos]
n_selected   <- length(selected_ids)

# Format a quantile probability as "Q50", "Q99", "Q99.5", etc. — handles both
# integer percentiles (0.50, 0.99) and fractional ones (0.995).
fmt_q_label <- function(p) {
  pct <- p * 100
  if (isTRUE(all.equal(pct, round(pct)))) {
    sprintf("Q%d", as.integer(round(pct)))
  } else {
    sprintf("Q%g", pct)
  }
}

message("  Selected stress-test year IDs:")
for (i in seq_along(QUANTILE_PROBS)) {
  yr  <- selected_ids[i]
  row <- compound_rated[compound_rated$sim_year == yr, ]
  message(sprintf(
    paste0(
      "    %-5s: year %4d | focal=%s  peak=%3.0f kt  follow_on=%d  ",
      "compound_dmg=%.4f  annual_dmg=%.4f  n_storms=%d"
    ),
    fmt_q_label(QUANTILE_PROBS[i]), yr,
    baseline_pin_map[[as.character(yr)]],
    row$focal_peak_wind_kt, row$compound_n_aftermath_events,
    row$compound_cum_damage, row$annual_total_damage, row$n_total
  ))
}

# --- Extract daily series for selected years ---
baseline_stress <- lapply(daily_base_list, function(df) {
  df[df$sim_year %in% selected_ids, ]
})

# --- Focal SIDs for replay ---
# Taken directly from the pin map rather than re-derived from the daily
# output: the pin IS the focal by design, even if a non-focal storm in some
# year happened to produce a slightly higher instantaneous site-level wind.
focal_sids <- baseline_pin_map[as.character(selected_ids)]

message("\nFocal SIDs pinned for Part 2 replay:")
for (i in seq_along(selected_ids)) {
  yr  <- selected_ids[i]
  sid <- focal_sids[[as.character(yr)]]
  message(sprintf("    %-5s (year %4d): %s",
                  fmt_q_label(QUANTILE_PROBS[i]), yr, sid))
}

# =============================================================================
# Part 2: KNMI'23 scenarios - anchored replay of selected years
#
# Each replay run draws n_selected years under the scenario's climate-adjusted
# rate, with the same focal SID pinned into each slot. Replay slots are then
# remapped onto the baseline selected_ids for direct paired comparison.
# =============================================================================

message("\nPart 2: KNMI'23 scenarios at target year ", TARGET_YEAR)
message("  Scenarios     : ", paste(KNMI_SCENARIOS, collapse = ", "))
message("  Baseline years: ", paste(selected_ids, collapse = ", "))
message("  Replay slots  : ", n_selected)

cc_stress <- vector("list", length(KNMI_SCENARIOS))
names(cc_stress) <- KNMI_SCENARIOS

# Pin map for replay: slot i (1..n_selected) anchored by the same SID as
# baseline selected_ids[i]. After daily generation, slots are remapped back
# to selected_ids for comparison.
replay_pins <- stats::setNames(
  focal_sids[as.character(selected_ids)],
  as.character(seq_along(selected_ids))
)

# Replay jitter: for each selected year, pass the SAME jitter spec that was
# applied in the baseline run so the pinned event has identical timing/
# intensity/RMW perturbation under baseline and scenarios. Only the climate-
# scenario delta_sst differs across baseline and KNMI runs.
replay_jitter <- stats::setNames(
  baseline_pin_jitter[as.character(selected_ids)],
  as.character(seq_along(selected_ids))
)

for (scen in KNMI_SCENARIOS) {

  message("  [", scen, "] Running model ...")

  cfg_cc <- make_hazard_cfg(
    data_path        = data_path,
    simulation_years = n_selected,
    climate          = make_climate_cfg(scenario = scen,
                                        target_year = TARGET_YEAR,
                                        perturb = TRUE),
    background_wind  = bg_cfg
  )

  out_cc <- run_hazard_model(cfg_cc, targets = targets, seed = SEED, verbose = FALSE)

  cc_daily_all <- generate_daily_hazard_impact_spatial(
    out         = out_cc,
    location    = targets$name,
    sim_years   = seq_len(n_selected),
    year0       = 2000L,
    gust_factor = 1.0,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = scen,
    seed        = SEED,
    pinned_sids = replay_pins,
    pin_jitter  = replay_jitter
  )

  cc_daily_all <- ipdcstorm:::.remap_replay_years(
    daily_list   = cc_daily_all,
    baseline_ids = selected_ids
  )

  # Apply the same state-dependent damage amplification as the baseline.
  cc_daily_all <- ipdcstorm:::.apply_state_dependent_damage(
    daily_list = cc_daily_all,
    alpha      = DAMAGE_ALPHA,
    cap        = DAMAGE_V_CAP
  )

  # Verify replay pins landed in the remapped output
  replay_check_map <- stats::setNames(
    focal_sids[as.character(selected_ids)],
    as.character(selected_ids)
  )
  missed_cc <- ipdcstorm:::.verify_pins_landed(cc_daily_all, replay_check_map)
  if (nrow(missed_cc) > 0L) {
    warning("[", scen, "] ", nrow(missed_cc), " of ", n_selected,
            " replay years did not contain the pinned focal SID.",
            call. = FALSE)
  }

  cc_paths <- save_daily_hazard_csvs(
    daily    = cc_daily_all,
    scenario = scen,
    out_dir  = OUTPUT_DIR
  )
  message("  [", scen, "] Wrote ", length(cc_paths), " CSVs to ", OUTPUT_DIR)

  cc_stress[[scen]] <- cc_daily_all
}

# =============================================================================
# Summary
# =============================================================================

# =============================================================================
# Diagnostic: scenario comparison for selected stress-test years
# =============================================================================
#
# Builds a tidy per-(year, location, scenario) summary table using the same
# compound window as the selection step, and produces overlay plots
# of daily wind so the climate signal on the pinned focal event is visible.

message("\nDiagnostic: building scenario comparison table ...")

daily_by_scen <- list(
  baseline = baseline_stress,
  knmi_Ld  = cc_stress[["knmi_Ld"]],
  knmi_Hd  = cc_stress[["knmi_Hd"]]
)

scen_rows <- list()
for (scen in names(daily_by_scen)) {
  daily_list <- daily_by_scen[[scen]]
  for (loc in targets$name) {
    m <- ipdcstorm:::compute_compound_stress_year_metrics(
      daily       = stats::setNames(list(daily_list[[loc]]), loc),
      sim_years   = selected_ids,
      location    = loc,
      window_days = COMPOUND_WINDOW_DAYS
    )
    if (nrow(m) == 0L) next
    m$scenario <- scen
    m$location <- loc
    m$pinned_sid <- unlist(focal_sids[as.character(m$sim_year)])
    scen_rows[[length(scen_rows) + 1L]] <- m
  }
}

scenario_comparison <- dplyr::bind_rows(scen_rows) |>
  dplyr::select(
    sim_year, location, scenario, pinned_sid,
    focal_event_id, focal_peak_wind_kt,
    compound_n_aftermath_events, compound_cum_damage,
    focal_start_doy, focal_end_doy, compound_window_end_doy
  ) |>
  dplyr::arrange(sim_year, location,
                 factor(scenario, levels = c("baseline", "knmi_Ld", "knmi_Hd")))

comp_path <- file.path(OUTPUT_DIR, "scenario_comparison.csv")
readr::write_csv(scenario_comparison, comp_path)
message("  Wrote scenario comparison to ", comp_path)

# Compact wide tables for at-a-glance climate-signal inspection
message("\nFocal-event peak wind (kt) by year x location x scenario:")
peak_wide <- scenario_comparison |>
  dplyr::select(sim_year, location, scenario, focal_peak_wind_kt) |>
  tidyr::pivot_wider(names_from = scenario, values_from = focal_peak_wind_kt) |>
  dplyr::arrange(sim_year, location)
print(peak_wide, n = nrow(peak_wide))

message(sprintf(
  "\nCompound cumulative damage (%dd) by year x location x scenario:",
  COMPOUND_WINDOW_DAYS
))
dmg_wide <- scenario_comparison |>
  dplyr::select(sim_year, location, scenario, compound_cum_damage) |>
  tidyr::pivot_wider(names_from = scenario, values_from = compound_cum_damage) |>
  dplyr::arrange(sim_year, location)
print(dmg_wide, n = nrow(dmg_wide))

# Per-location overlay plots: daily wind across the 3 scenarios, faceted by
# selected year. Dates are normalised to day-of-year so scenarios overlay
# cleanly within each year (cc_stress uses different calendar years).
# Axis ranges are fixed across locations so plots are directly comparable.
message("\nGenerating per-location daily-wind overlay plots ...")
plot_dir <- file.path(OUTPUT_DIR, "..", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

year_labels <- stats::setNames(
  paste0("Q", round(QUANTILE_PROBS * 100, 2), " (year ", selected_ids, ")"),
  as.character(selected_ids)
)

scen_colours <- c(baseline = "black",
                  knmi_Ld  = "#1f77b4",
                  knmi_Hd  = "#d62728")

# Saffir-Simpson thresholds for site-level wind context.
TS_THRESHOLD_KT  <- 34L    # tropical-storm force
HUR_THRESHOLD_KT <- 64L    # Cat-1 hurricane force

# Fixed axes across locations so panels are directly comparable.
#
# Linear wind axis ceiling chosen to accommodate perturbed Irma-class peaks
# under knmi_Hd (~170 kt). Cumulative damage ceiling sized for Lenny-anchored
# years amplified by HAZUS-style state dependency under knmi_Hd (~0.10-0.12).
WIND_Y_MAX        <- 180L
WIND_Y_STEP       <- 30L
CUM_DAMAGE_Y_MAX  <- 0.15
CUM_DAMAGE_Y_STEP <- 0.025

# Toggle log10 y-axis for the daily wind plot. On log scale the threshold
# lines (TS, HUR) remain horizontal and the upper tail (major hurricanes)
# is compressed; the lower tail (background trade-wind days) is expanded.
# Useful for examining background-wind structure side-by-side with peak
# events; less useful for reading off peak-wind values.
WIND_LOG_Y       <- FALSE
WIND_LOG_Y_MIN   <- 1     # kt - log floor; days with wind < 1 kt are clipped
WIND_LOG_Y_MAX   <- 200   # kt - log ceiling above the strongest perturbed peak

# Month-start day-of-year ticks (non-leap reference).
doy_month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
doy_month_labels <- month.abb

for (loc in targets$name) {
  long_loc <- dplyr::bind_rows(lapply(names(daily_by_scen), function(scen) {
    df <- daily_by_scen[[scen]][[loc]]
    df$scenario <- scen
    df
  }))
  long_loc <- long_loc[long_loc$sim_year %in% selected_ids, , drop = FALSE]
  long_loc$doy <- as.integer(long_loc$doy)
  long_loc$year_label <- factor(
    year_labels[as.character(long_loc$sim_year)],
    levels = year_labels
  )

  if (WIND_LOG_Y) {
    # Log10 y axis: clamp wind to [WIND_LOG_Y_MIN, WIND_LOG_Y_MAX] and use
    # decade-aware breaks. The TS/HUR threshold lines remain horizontal and
    # split the visual into background-wind (below 34 kt) and storm regimes.
    long_loc_wind <- long_loc
    long_loc_wind$wind_kt <- pmin(WIND_LOG_Y_MAX,
                                  pmax(WIND_LOG_Y_MIN, long_loc_wind$wind_kt))
    wind_y_breaks <- c(1, 5, 10, 20, TS_THRESHOLD_KT,
                       HUR_THRESHOLD_KT, 100, 150, 200)
    wind_y_scale <- scale_y_log10(
      limits = c(WIND_LOG_Y_MIN, WIND_LOG_Y_MAX),
      breaks = wind_y_breaks,
      labels = wind_y_breaks,
      expand = c(0, 0)
    )
    wind_y_label <- "wind (kt, log scale)"
    wind_subtitle <- paste0(
      "baseline vs KNMI'23 @ ", TARGET_YEAR,
      " (focal SID pinned in every year). ",
      "Dashed lines: TS 34 kt (orange), HUR 64 kt (red). ",
      "Wind clamped to [", WIND_LOG_Y_MIN, ", ", WIND_LOG_Y_MAX, "] kt."
    )
  } else {
    long_loc_wind <- long_loc
    wind_y_scale <- scale_y_continuous(
      limits = c(0, WIND_Y_MAX),
      breaks = seq(0, WIND_Y_MAX, by = WIND_Y_STEP),
      expand = c(0, 0)
    )
    wind_y_label <- "wind (kt)"
    wind_subtitle <- paste0(
      "baseline vs KNMI'23 @ ", TARGET_YEAR,
      "Dashed lines: TS 34 kt (orange), HUR 64 kt (red)"
    )
  }

  season_band <- data.frame(xmin = 152, xmax = 335, ymin = -Inf, ymax = Inf)

  p_wind <- ggplot(long_loc_wind,
                            aes(doy, wind_kt)) +
    geom_hline(yintercept = TS_THRESHOLD_KT,
                        linetype = "dashed", colour = "#ff7f00",
                        linewidth = 0.7, alpha = 0.8) +
    geom_hline(yintercept = HUR_THRESHOLD_KT,
                        linetype = "dashed", colour = "#e31a1c",
                        linewidth = 0.7, alpha = 0.8) +
    geom_line(aes(colour = scenario), linewidth = 0.7, alpha = 0.85) +
    geom_rect(
      data = season_band,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "green", alpha = 0.1
    ) +
    facet_wrap(~ year_label, ncol = 1) +
    scale_colour_manual(values = scen_colours) +
    scale_x_continuous(
      breaks = doy_month_breaks,
      labels = doy_month_labels,
      limits = c(1, 366),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0,160,40),
      limits = c(0, 160),
      expand = c(0, 0)) +
    labs(
      title = paste0("Daily wind: ", loc, " - selected synthetic years"),
      subtitle = wind_subtitle,
      x = NULL, y = wind_y_label, colour = "scenario"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  png_wind <- file.path(plot_dir, paste0("scenario_compare__", loc, "__wind.png"))
  ggsave(png_wind, p_wind, width = 10, height = 12, dpi = 120)

  # Cumulative damage overlay - integrated impact through the season.
  # Within each facet (year) the curve is monotonically non-decreasing and
  # rises in steps when a storm strikes. The vertical separation between
  # baseline / knmi_Ld / knmi_Hd at year-end shows the climate-driven
  # increment in seasonal loss.
  p_damage <- ggplot(long_loc,
                              aes(doy, cum_damage, colour = scenario)) +
    geom_line(linewidth = 0.6, alpha = 0.9) +
    facet_wrap(~ year_label, ncol = 1) +
    scale_colour_manual(values = scen_colours) +
    scale_x_continuous(
      breaks = doy_month_breaks,
      labels = doy_month_labels,
      limits = c(1, 366),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      limits = c(0, CUM_DAMAGE_Y_MAX),
      breaks = seq(0, CUM_DAMAGE_Y_MAX, by = CUM_DAMAGE_Y_STEP),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Cumulative damage: ", loc, " - selected stress years"),
      subtitle = paste0("baseline vs KNMI'23 @ ", TARGET_YEAR,
                        " (focal SID pinned in every year). ",
                        "State-dependent amplification: alpha = ",
                        DAMAGE_ALPHA, ", cap = ", DAMAGE_V_CAP, "."),
      x = NULL,
      y = "cumulative damage (fraction of total loss)",
      colour = "scenario"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  png_damage <- file.path(plot_dir,
                          paste0("scenario_compare__", loc, "__cum_damage.png"))
  ggsave(png_damage, p_damage, width = 10, height = 12, dpi = 120)
}
message("  Wrote per-location wind + damage plots to ", plot_dir)

# --- Baseline stress distribution plot (all 1000 years, Qs highlighted) ------
message("\nGenerating baseline compound-stress distribution plot ...")

rank_dist <- compound_rated |>
  dplyr::mutate(
    rank_pct   = 100 * rating,
    pinned_sid = unname(unlist(baseline_pin_map[as.character(sim_year)])),
    pinned     = factor(PINNED_SID_NAMES[pinned_sid],
                        levels = PINNED_SID_NAMES[PINNED_SID_POOL])
  )

quantile_rows <- rank_dist |>
  dplyr::filter(sim_year %in% selected_ids) |>
  dplyr::mutate(
    quantile_label = paste0(
      "Q", round(QUANTILE_PROBS[match(sim_year, selected_ids)] * 100,2)
    ),
    point_label = paste0(quantile_label, " (yr ", sim_year, ")")
  ) |>
  dplyr::arrange(rating)

pool_palette <- stats::setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  PINNED_SID_NAMES[PINNED_SID_POOL]
)

p_dist <- ggplot(
  rank_dist,
  aes(rank_pct, compound_cum_damage)
) +
  geom_vline(
    data = quantile_rows, aes(xintercept = rank_pct),
    linetype = "dashed", colour = "grey40", linewidth = 0.3, alpha = 0.6
  ) +
  geom_point(size = 0.9, alpha = 0.55) +
  geom_point(
    data = quantile_rows,
    aes(rank_pct, compound_cum_damage),
    shape = 21, fill = "white", colour = "black",
    size = 4.2, stroke = 1.2
  ) +
  geom_text(
    data = quantile_rows,
    aes(
      rank_pct, compound_cum_damage,
      label = point_label
    ),
    vjust = -1.2, hjust = 0.5, size = 4, fontface = "bold"
  ) +
  scale_colour_manual(values = pool_palette, name = "Pinned focal") +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),
    limits = c(0, 103),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(rank_dist$compound_cum_damage, na.rm = TRUE) * 1.15),
    expand = c(0, 0),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    x = "Rating percentile",
    y = "Compound cumulative damage (60d, Saba + Statia)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  guides(colour = guide_legend(
    nrow = 1, override.aes = list(size = 3, alpha = 1)
  ))

png_dist <- file.path(plot_dir, "baseline_stress_distribution.png")
ggsave(png_dist, p_dist, width = 12, height = 7, dpi = 120)
message("  Wrote ", png_dist)

message("\n--- Stress-test experiment complete ---")
message("Seed               : ", SEED)
message("Focal pool         : ", paste(PINNED_SID_POOL, collapse = ", "))
message("Baseline years     : ", N_BASE)
message("Replay years       : ", n_selected)
message("Target year        : ", TARGET_YEAR)
message("Rating metric      : compound_cum_damage (",
        COMPOUND_WINDOW_DAYS, "-day window at Saba + Statia)")
message("Selected year IDs  : ", paste(selected_ids, collapse = ", "),
        "  [", paste0("Q", round(QUANTILE_PROBS * 100, 2), collapse = " / "), "]")
message("CC scenarios       : ", paste(KNMI_SCENARIOS, collapse = ", "))
message("Raw CSV output     : ", OUTPUT_DIR)
message("Experiment type    : anchored compound stress years under climate forcing")
message("")
message("Objects in workspace:")
message("  selected_ids        -  integer[5]: baseline year IDs for all scenario comparisons")
message("  focal_sids          -  list[sim_year]: pinned focal SID per selected year")
message("  baseline_stress     -  list[location]: daily tibbles for 5 baseline years")
message("  cc_stress           -  list[scenario][location]: replayed daily tibbles keyed to baseline year IDs")
message("  scenario_comparison -  tidy tibble: per-(year, location, scenario) compound metrics")
