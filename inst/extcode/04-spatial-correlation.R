# =============================================================================
# Spatial correlation of synthetic daily wind series
#
# Purpose:
#   Assess how correlated the synthetic daily wind series are across the
#   target location portfolio, focusing on two complementary aspects:
#
#   1) Co-occurrence rates — how often do pairs of locations exceed a wind
#      threshold on the same calendar day, and how does that change from
#      tropical-storm to hurricane intensity?
#
#   2) Event-sharing rates — what fraction of simulation years contain the
#      same storm event at more than one location, and which location pairs
#      are most frequently co-affected by the same event?
#
# Scientific note:
#   The daily series are heavily zero-inflated (most days have no TC
#   activity).  Standard Pearson correlation is dominated by shared calm days
#   and is not a useful measure of spatial dependence for this type of data.
#   Co-occurrence and event-sharing metrics bypass the zero-inflation problem
#   by operating directly on exceedance or event-identity information.
#
# Outputs (printed and/or returned as objects for downstream use):
#   - co_occurrence_tbl   : pairwise daily co-occurrence rates at TS / hurricane
#   - cond_exceed_tbl     : conditional exceedance P(B >= T | A >= T)
#   - annual_co_occ_tbl   : fraction of sim_years with joint annual exceedance
#   - event_share_tbl     : pairwise event-sharing rates across all events
#   - shared_events_tbl   : year-level log of shared events per location pair
# =============================================================================

# %%
library(ipdcstorm)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# 1) Global settings
# =============================================================================

ibtracs_file_path <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
seed              <- 2026L
simulation_years  <- 2000L
year0             <- 2025L

# Wind thresholds used throughout
thresh_ts  <- 34L   # tropical-storm force (kt)
thresh_hur <- 64L   # hurricane force (kt)

# %%

# =============================================================================
# 2) Target locations
# =============================================================================

targets <- tibble::tribble(
  ~name         , ~lat    , ~lon     ,
  "St_Martin"   , 18.0708 , -63.0501 ,
  "Saba"        , 17.6350 , -63.2300 ,
  "Statia"      , 17.4890 , -62.9740 ,
  "Puerto_Rico" , 18.2208 , -66.5901 ,
  "Miami"       , 25.7617 , -80.1918
)

# %%

# =============================================================================
# 3) Run hazard model and generate daily series
# =============================================================================
# (Skip if hazard_out and daily_out are already in the workspace from a
#  previous script run — e.g. after sourcing 03-stress-test-forcings.R)

if (!exists("hazard_out")) {
  hazard_cfg <- make_hazard_cfg(
    data_path             = ibtracs_file_path,
    search_radius_km      = 600,
    historical_start_year = 1970L,
    simulation_years      = simulation_years,
    climate = make_climate_cfg(
      scenario    = "stationary",
      delta_sst   = 0,
      target_year = year0
    )
  )
  hazard_out <- run_hazard_model(
    cfg     = hazard_cfg,
    targets = targets,
    seed    = seed,
    verbose = TRUE
  )
}

if (!exists("daily_out")) {
  daily_out <- generate_daily_hazard_impact(
    out         = hazard_out,
    location    = targets$name,
    sim_years   = seq_len(simulation_years),
    year0       = year0,
    gust_factor = 1.3,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = "stationary"
  )
}

# Flatten to a single tibble for pairwise operations
daily_flat <- dplyr::bind_rows(daily_out)

locs   <- targets$name
n_locs <- length(locs)

# %%

# =============================================================================
# 4) Daily co-occurrence rates
# =============================================================================
#
# For each ordered pair (loc_a, loc_b) and threshold T, we compute:
#
#   joint_rate      = P(A >= T  and  B >= T)    [fraction of all days]
#   cond_rate_a_b   = P(B >= T  |    A >= T)    [given A exceeds, does B?]
#   cond_rate_b_a   = P(A >= T  |    B >= T)    [given B exceeds, does A?]
#
# joint_rate tells you how common simultaneous hits are in absolute terms.
# cond_rate captures the spatial association: conditional on one site being
# affected, how likely is the other site to also be affected?
#
# We compute this for both TS (34 kt) and hurricane (64 kt) thresholds.

# --- Pivot to wide format: one row per (sim_year, date), one col per location
daily_wide <- daily_flat |>
  dplyr::select("location", "sim_year", "date", "wind_kt") |>
  tidyr::pivot_wider(names_from = "location", values_from = "wind_kt",
                     values_fill = 0)

n_days <- nrow(daily_wide)

# Helper: compute pairwise co-occurrence stats for a given threshold
compute_co_occurrence <- function(wide, locs, threshold, n_days) {
  pairs <- expand.grid(loc_a = locs, loc_b = locs, stringsAsFactors = FALSE) |>
    dplyr::filter(loc_a != loc_b)

  results <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$loc_a[i]
    b <- pairs$loc_b[i]

    exceed_a  <- wide[[a]] >= threshold
    exceed_b  <- wide[[b]] >= threshold
    joint     <- exceed_a & exceed_b

    n_a    <- sum(exceed_a)
    n_b    <- sum(exceed_b)
    n_both <- sum(joint)

    results[[i]] <- tibble::tibble(
      loc_a      = a,
      loc_b      = b,
      threshold  = as.integer(threshold),
      n_days_a   = n_a,
      n_days_b   = n_b,
      n_days_both = n_both,
      joint_rate  = n_both / n_days,
      # P(B >= T | A >= T)
      cond_rate   = if (n_a > 0L) n_both / n_a else NA_real_
    )
  }
  dplyr::bind_rows(results)
}

co_occ_ts  <- compute_co_occurrence(daily_wide, locs, thresh_ts,  n_days)
co_occ_hur <- compute_co_occurrence(daily_wide, locs, thresh_hur, n_days)
co_occurrence_tbl <- dplyr::bind_rows(co_occ_ts, co_occ_hur)

cat("\n=== Daily co-occurrence rates (joint) ===\n")
co_occurrence_tbl |>
  dplyr::mutate(
    label      = paste0(loc_a, " & ", loc_b),
    threshold  = paste0(threshold, " kt"),
    joint_pct  = round(joint_rate * 100, 3),
    cond_pct   = round(cond_rate  * 100, 2)
  ) |>
  dplyr::select(loc_a, loc_b, threshold, joint_pct, cond_pct) |>
  dplyr::arrange(threshold, dplyr::desc(joint_pct)) |>
  print(n = 40)

# %%

# =============================================================================
# 5) Annual co-occurrence rates
# =============================================================================
#
# Collapsing to annual maxima removes the daily zero-inflation entirely.
# For each sim_year and location we compute the annual peak wind; then ask:
# what fraction of simulation years saw BOTH locations exceed the threshold?
#
# This is the metric most directly comparable to historical event records,
# where we know (e.g.) "Irma hit both St Martin and Saba in 2017".

annual_max <- daily_flat |>
  dplyr::group_by(.data$location, .data$sim_year) |>
  dplyr::summarise(annual_peak = max(.data$wind_kt, na.rm = TRUE),
                   .groups = "drop")

annual_wide <- annual_max |>
  tidyr::pivot_wider(names_from = "location", values_from = "annual_peak",
                     values_fill = 0)

n_years <- nrow(annual_wide)

compute_annual_co_occ <- function(ann_wide, locs, threshold, n_years) {
  pairs <- expand.grid(loc_a = locs, loc_b = locs, stringsAsFactors = FALSE) |>
    dplyr::filter(loc_a != loc_b)

  results <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$loc_a[i]
    b <- pairs$loc_b[i]

    exceed_a  <- ann_wide[[a]] >= threshold
    exceed_b  <- ann_wide[[b]] >= threshold
    n_both    <- sum(exceed_a & exceed_b)
    n_a       <- sum(exceed_a)

    results[[i]] <- tibble::tibble(
      loc_a             = a,
      loc_b             = b,
      threshold         = as.integer(threshold),
      annual_joint_rate = n_both / n_years,
      # P(both exceed in same year | loc_a exceeds)
      annual_cond_rate  = if (n_a > 0L) n_both / n_a else NA_real_
    )
  }
  dplyr::bind_rows(results)
}

annual_co_occ_tbl <- dplyr::bind_rows(
  compute_annual_co_occ(annual_wide, locs, thresh_ts,  n_years),
  compute_annual_co_occ(annual_wide, locs, thresh_hur, n_years)
)

cat("\n=== Annual co-occurrence rates ===\n")
annual_co_occ_tbl |>
  dplyr::mutate(
    threshold       = paste0(threshold, " kt"),
    annual_joint_pct = round(annual_joint_rate * 100, 2),
    annual_cond_pct  = round(annual_cond_rate  * 100, 1)
  ) |>
  dplyr::select(loc_a, loc_b, threshold, annual_joint_pct, annual_cond_pct) |>
  dplyr::arrange(threshold, dplyr::desc(annual_joint_pct)) |>
  print(n = 40)

# %%

# =============================================================================
# 6) Event-sharing rates
# =============================================================================
#
# Two locations "share an event" in a given sim_year if the same event_id
# appears in both daily series for that year.  This is a structural check on
# the spatial footprint of individual storms: does a single storm in the
# simulation propagate wind to multiple sites simultaneously?
#
# We report three metrics per location pair:
#
#   share_rate         = fraction of sim_years with >= 1 shared event
#   mean_shared_events = average number of distinct shared events per year
#                        (years with no shared event contribute 0)
#   p_multi_shared     = fraction of sim_years with >= 2 shared events
#                        (consecutive hits in the same season)
#
# We also produce shared_events_tbl: the full year-level log of which
# event_ids were shared, useful for downstream inspection.

# Extract active days only (event_id is not NA)
active_days <- daily_flat |>
  dplyr::filter(!is.na(.data$event_id)) |>
  dplyr::select("location", "sim_year", "event_id") |>
  dplyr::distinct()

# For each location pair, find shared event_ids per sim_year
pairs_df <- expand.grid(loc_a = locs, loc_b = locs, stringsAsFactors = FALSE) |>
  dplyr::filter(loc_a < loc_b)   # unordered pairs only (a < b lexicographically)

shared_events_list <- vector("list", nrow(pairs_df))

for (i in seq_len(nrow(pairs_df))) {
  a <- pairs_df$loc_a[i]
  b <- pairs_df$loc_b[i]

  events_a <- active_days |> dplyr::filter(.data$location == a)
  events_b <- active_days |> dplyr::filter(.data$location == b)

  # Inner join on (sim_year, event_id) — rows that exist at both locations
  shared <- dplyr::inner_join(
    events_a |> dplyr::select("sim_year", "event_id"),
    events_b |> dplyr::select("sim_year", "event_id"),
    by = c("sim_year", "event_id")
  ) |>
    dplyr::mutate(loc_a = a, loc_b = b) |>
    dplyr::select("loc_a", "loc_b", "sim_year", "event_id")

  shared_events_list[[i]] <- shared
}

shared_events_tbl <- dplyr::bind_rows(shared_events_list)

# Summarise to pairwise annual counts
shared_annual <- shared_events_tbl |>
  dplyr::count(.data$loc_a, .data$loc_b, .data$sim_year,
               name = "n_shared_events")

# Roll up to pair-level stats across all sim_years
event_share_tbl <- shared_annual |>
  # Add zero-rows for years with no shared event
  tidyr::complete(
    tidyr::nesting(loc_a, loc_b),
    sim_year = seq_len(simulation_years),
    fill = list(n_shared_events = 0L)
  ) |>
  dplyr::group_by(.data$loc_a, .data$loc_b) |>
  dplyr::summarise(
    share_rate         = mean(.data$n_shared_events > 0),
    mean_shared_events = mean(.data$n_shared_events),
    p_multi_shared     = mean(.data$n_shared_events > 1),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data$share_rate))

cat("\n=== Event-sharing rates (all sim_years) ===\n")
event_share_tbl |>
  dplyr::mutate(
    share_pct        = round(share_rate         * 100, 1),
    mean_shared_evts = round(mean_shared_events,        2),
    multi_pct        = round(p_multi_shared      * 100, 1)
  ) |>
  dplyr::select(loc_a, loc_b, share_pct, mean_shared_evts, multi_pct) |>
  print(n = 20)

# %%

# =============================================================================
# 7) Most frequently shared individual storms
# =============================================================================
#
# Identifies which specific historical storm events (by event_id prefix =
# original IBTrACS SID) are most commonly shared across multiple locations.
# A high-frequency shared storm is a persistent basin-wide threat.

# Extract the IBTrACS SID from event_id ("<SID>_y<year>_<counter>")
shared_events_annotated <- shared_events_tbl |>
  dplyr::mutate(
    source_sid = sub("_y[0-9]+_.*$", "", .data$event_id)
  )

top_shared_storms <- shared_events_annotated |>
  dplyr::count(.data$loc_a, .data$loc_b, .data$source_sid,
               name = "n_years_shared") |>
  dplyr::group_by(.data$loc_a, .data$loc_b) |>
  dplyr::slice_max(.data$n_years_shared, n = 5, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$loc_a, .data$loc_b, dplyr::desc(.data$n_years_shared))

cat("\n=== Top 5 most-shared source storms per location pair ===\n")
print(top_shared_storms, n = 50)

# %%

# =============================================================================
# 8) Visualisation — co-occurrence heatmaps
# =============================================================================

# Helper: pivot pairwise stat to a symmetric matrix for heatmap display
to_sym_matrix <- function(pair_tbl, loc_a, loc_b, value_col, all_locs) {
  # Both directions: (a, b) and (b, a)
  fwd <- pair_tbl |>
    dplyr::select(row = {{ loc_a }}, col = {{ loc_b }},
                  value = {{ value_col }})
  rev <- pair_tbl |>
    dplyr::select(row = {{ loc_b }}, col = {{ loc_a }},
                  value = {{ value_col }})
  dplyr::bind_rows(fwd, rev) |>
    dplyr::filter(.data$row != .data$col) |>
    dplyr::bind_rows(
      tibble::tibble(row = all_locs, col = all_locs, value = NA_real_)
    ) |>
    dplyr::mutate(
      row = factor(.data$row, levels = all_locs),
      col = factor(.data$col, levels = all_locs)
    )
}

# --- 8a: Daily co-occurrence at hurricane threshold (cond_rate) -------------

hur_cond <- co_occurrence_tbl |>
  dplyr::filter(.data$threshold == thresh_hur)

mat_hur_cond <- to_sym_matrix(hur_cond, loc_a, loc_b, cond_rate, locs)

p_daily_cond <- ggplot2::ggplot(
  mat_hur_cond,
  ggplot2::aes(x = col, y = row, fill = value)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(!is.na(value), scales::percent(value, 1), "")),
    size = 3.2
  ) +
  ggplot2::scale_fill_distiller(
    palette  = "YlOrRd",
    direction = 1,
    na.value = "grey92",
    limits   = c(0, 1),
    labels   = scales::percent_format(1),
    name     = "P(B ≥ 64 kt\n| A ≥ 64 kt)"
  ) +
  ggplot2::labs(
    title    = "Conditional daily co-occurrence — hurricane force (≥ 64 kt)",
    subtitle = sprintf("P(column site exceeds | row site exceeds)  ·  %d simulation years",
                       simulation_years),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

print(p_daily_cond)

# --- 8b: Annual co-occurrence at hurricane threshold ------------------------

ann_hur <- annual_co_occ_tbl |>
  dplyr::filter(.data$threshold == thresh_hur)

mat_ann_hur <- to_sym_matrix(ann_hur, loc_a, loc_b, annual_joint_rate, locs)

p_annual_joint <- ggplot2::ggplot(
  mat_ann_hur,
  ggplot2::aes(x = col, y = row, fill = value)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(!is.na(value), scales::percent(value, 1), "")),
    size = 3.2
  ) +
  ggplot2::scale_fill_distiller(
    palette   = "Blues",
    direction = 1,
    na.value  = "grey92",
    limits    = c(0, NA),
    labels    = scales::percent_format(1),
    name      = "Joint annual\nexceedance rate"
  ) +
  ggplot2::labs(
    title    = "Annual joint exceedance rate — hurricane force (≥ 64 kt)",
    subtitle = sprintf("Fraction of years in which both sites exceed threshold  ·  %d simulation years",
                       simulation_years),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

print(p_annual_joint)

# --- 8c: Event-sharing rate heatmap ----------------------------------------

mat_share <- to_sym_matrix(
  event_share_tbl, loc_a, loc_b, share_rate, locs
)

p_share <- ggplot2::ggplot(
  mat_share,
  ggplot2::aes(x = col, y = row, fill = value)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(!is.na(value), scales::percent(value, 1), "")),
    size = 3.2
  ) +
  ggplot2::scale_fill_distiller(
    palette   = "Greens",
    direction = 1,
    na.value  = "grey92",
    limits    = c(0, NA),
    labels    = scales::percent_format(1),
    name      = "Event-sharing\nrate"
  ) +
  ggplot2::labs(
    title    = "Event-sharing rate — fraction of years with ≥ 1 shared storm",
    subtitle = sprintf("Same event_id appearing at both sites  ·  %d simulation years",
                       simulation_years),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

print(p_share)

# %%

# =============================================================================
# 9) Summary comparison table
# =============================================================================
#
# Merges daily co-occurrence, annual co-occurrence, and event-sharing into a
# single wide table for quick inspection or export.

# For unordered pairs only (loc_a < loc_b)
daily_hur_cond_asym <- co_occurrence_tbl |>
  dplyr::filter(.data$threshold == thresh_hur, .data$loc_a < .data$loc_b) |>
  dplyr::select(loc_a, loc_b,
                daily_joint_rate_hur = joint_rate,
                daily_cond_rate_hur  = cond_rate)

annual_hur_asym <- annual_co_occ_tbl |>
  dplyr::filter(.data$threshold == thresh_hur, .data$loc_a < .data$loc_b) |>
  dplyr::select(loc_a, loc_b,
                annual_joint_rate_hur = annual_joint_rate,
                annual_cond_rate_hur  = annual_cond_rate)

summary_tbl <- event_share_tbl |>
  dplyr::left_join(daily_hur_cond_asym,  by = c("loc_a", "loc_b")) |>
  dplyr::left_join(annual_hur_asym,      by = c("loc_a", "loc_b")) |>
  dplyr::mutate(
    pair = paste0(loc_a, " / ", loc_b),
    across(where(is.numeric), ~ round(.x, 4))
  ) |>
  dplyr::select(
    pair,
    share_rate, mean_shared_events, p_multi_shared,
    daily_joint_rate_hur, daily_cond_rate_hur,
    annual_joint_rate_hur, annual_cond_rate_hur
  ) |>
  dplyr::arrange(dplyr::desc(share_rate))

cat("\n=== Summary: all pairwise spatial correlation metrics ===\n")
print(summary_tbl, n = 20)
