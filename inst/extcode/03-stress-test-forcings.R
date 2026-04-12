# =============================================================================
# Stress-test forcing pipeline
#
# Purpose:
#   Extract targeted stress-test forcing years from the stochastic ensemble
#   for use in downstream impact modelling (e.g. economic loss, infrastructure
#   resilience, insurance pricing).
#
# What this script does:
#   1) Configure and run the baseline stochastic hazard model
#   2) Generate a full synthetic daily hazard-impact series
#   3) Track-based selection  — identify years where a specific historical
#      storm (e.g. Irma) was resampled from the event library, preserving its
#      original approach geometry and wind profile
#   4) Impact-based selection — identify years where any storm produced at
#      least the same site-level impact as that reference storm, regardless of
#      which storm caused it
#   5) Compare the two selection sets and summarise stress-test candidates
#
# Scientific note:
#   Track-based years are by definition rare (they contain exactly one storm
#   with Irma's geometry) and useful for deterministic stress tests.
#   Impact-based years are the broader exceedance set and are more appropriate
#   for frequency-of-exceedance or return-period analysis.
# =============================================================================

# %%
library(ipdcstorm)
library(dplyr)
library(ggplot2)

# =============================================================================
# 1) Global settings
# =============================================================================

ibtracs_file_path     <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
baseline_output_dir   <- "output/baseline"
validation_output_dir <- "output/validation"

seed             <- 2026L
simulation_years <- 2000L
year0            <- 2025L    # calendar year assigned to sim_year = 1

# Reference storm for stress-test selection (Hurricane Irma, 2017).
# Uses the IBTrACS native SID (YYYYDDDLLLBBB format), NOT the ATCF ID
# ("AL112017"). Use lookup_storm_id(hazard_out, year = 2017, min_wind_kt = 50)
# to verify or find SIDs for other storms.
irma_sid <- "2017242N16333"

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
# 3) Configure and run the hazard model
# =============================================================================

# Stationary baseline: no climate delta applied; all stochastic variability
# comes from Poisson sampling of historical rates and event-library resampling.
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

# seed is stored in hazard_out$run_metadata$seed and automatically inherited
# by generate_daily_hazard_impact() when called without an explicit seed.
hazard_out <- run_hazard_model(
  cfg     = hazard_cfg,
  targets = targets,
  seed    = seed,
  verbose = TRUE
)

# %%

# =============================================================================
# 4) Generate daily hazard-impact series
# =============================================================================

# Spatially coherent variant (Option 1 — shared event pool).
#
# Storms are drawn once at basin level per simulated year; each drawn storm is
# assigned to every location whose event library contains it. This enforces
# co-occurrence: if Irma (2017242N16333) is drawn in year 47, it appears at
# Saba, St. Martin, and any other location with Irma in its library, each with
# its own site-level wind profile.
#
# Drop-in replacement for generate_daily_hazard_impact(); output schema is
# identical. Seed is inherited from hazard_out$run_metadata$seed automatically.
daily_out <- generate_daily_hazard_impact_spatial(
  out         = hazard_out,
  location    = targets$name,
  sim_years   = seq_len(simulation_years),
  year0       = year0,
  gust_factor = 1.3,
  damage      = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario    = "stationary"
  # seed = NULL  ->  inherits hazard_out$run_metadata$seed automatically
)

# To switch back to independent (per-location) sampling, replace with:
# daily_out <- generate_daily_hazard_impact(
#   out = hazard_out, location = targets$name, sim_years = seq_len(simulation_years),
#   year0 = year0, gust_factor = 1.3, damage = list(method = "intensity"),
#   pulse_shape = "cosine", scenario = "stationary"
# )

# %%

# =============================================================================
# 5) Track-based stress selection
# =============================================================================
#
# Question: "Give me every synthetic year in which Hurricane Irma's exact
# track geometry was resampled — same approach angle, same temporal wind
# profile, same storm size."
#
# The event-library sampler encodes the source storm in each event_id as
# "<SID>_y<year>_<counter>", so this query matches by SID prefix.
# Results represent a subset of the ensemble with Irma's physical fingerprint.

track_years <- query_storm_track_years(
  daily    = daily_out,
  storm_id = irma_sid    # "AL112017" — Hurricane Irma 2017
)

# How many times was Irma's track resampled at each location?
track_freq <- track_years |>
  count(location, name = "n_irma_track_years") |>
  mutate(pct_of_sim = round(100 * n_irma_track_years / simulation_years, 2))

cat("\n--- Track-based selection: Irma resampling frequency ---\n")
print(track_freq)

# Inspect the actual sim_years for one focal location
irma_track_years_saba <- track_years |>
  filter(location == "Saba") |>
  pull(sim_year)

cat(sprintf(
  "\nSaba: Irma track appeared in %d of %d simulated years (%.2f%%)\n",
  length(irma_track_years_saba),
  simulation_years,
  100 * length(irma_track_years_saba) / simulation_years
))

# Extract the full daily series for those years at Saba — these are the
# stress-test forcing time series that retain Irma's exact wind geometry.
daily_irma_track_saba <- daily_out$Saba |>
  filter(sim_year %in% irma_track_years_saba)

# %%

# =============================================================================
# 6) Impact-based stress selection
# =============================================================================
#
# Question: "Give me every synthetic year in which a site experienced at least
# the same impact as the historical Irma passage — regardless of which storm
# caused it."
#
# Three complementary metrics are supported:
#   peak_wind_kt    maximum daily sustained wind in the year (kt)
#   cum_damage      total cumulative damage fraction over the year
#   max_damage_rate worst single-day damage rate
#
# Threshold derivation (hybrid):
#   ref_threshold   -> Irma's site wind from hazard_out$events (peak_wind_kt)
#                      or median metric across track-based years (damage metrics)
#   percentile gate -> empirical quantile of annual metric distribution
#   effective       -> max(ref_threshold, percentile_threshold)
#   min_threshold   -> optional absolute floor applied after the above
#
# Tune `percentile_gate` and `wind_floor_kt` in the parameters block below.

# --- Selection parameters (adjust to control stress-year frequency) ----------

# percentile_gate:  selection gate on the annual-metric distribution.
#   0.95 -> top 5% of years (~100 out of 2000); 0.90 -> top 10% (~200 / 2000).
#   Combined with the storm-based lower bound via max(), so both must be cleared.
percentile_gate <- 0.95

# wind_floor_kt: absolute lower bound on the peak-wind threshold.
#   NULL -> no floor (percentile gate + Irma's site wind determine threshold).
#   Set to e.g. 64 (Cat-1) or 83 (Cat-2) to enforce a physical minimum.
wind_floor_kt <- NULL

# --- 6a) Peak-wind exceedance ---

# Effective threshold per location = max(Irma's historical site wind,
#   percentile_gate quantile of annual peak wind).  wind_floor_kt is applied
#   afterwards as an absolute floor.
impact_years_wind <- query_impact_years(
  daily         = daily_out,
  storm_id      = irma_sid,
  out           = hazard_out,
  metric        = "peak_wind_kt",
  percentile    = percentile_gate,
  min_threshold = wind_floor_kt
)

wind_freq <- impact_years_wind |>
  count(location, name = "n_irma_impact_years") |>
  mutate(pct_of_sim = round(100 * n_irma_impact_years / simulation_years, 2))

cat("\n--- Impact-based selection: peak-wind exceedance frequency ---\n")
print(wind_freq)

# --- 6b) Cumulative damage exceedance ---

# Same percentile gate applied to cum_damage distribution.
impact_years_damage <- query_impact_years(
  daily      = daily_out,
  storm_id   = irma_sid,
  out        = hazard_out,
  metric     = "cum_damage",
  percentile = percentile_gate
)

damage_freq <- impact_years_damage |>
  count(location, name = "n_damage_years") |>
  mutate(pct_of_sim = round(100 * n_damage_years / simulation_years, 2))

cat("\n--- Impact-based selection: cumulative-damage exceedance frequency ---\n")
print(damage_freq)

# %%

# =============================================================================
# 7) Compare track-based vs impact-based sets (Saba focus)
# =============================================================================
#
# Track-based set:  years with Irma's geometry (may include weak-impact years
#                   if the storm passed at the edge of the search radius)
# Impact-based set: years with Irma-level wind from any storm
# Intersection:     the "worst of both worlds" — Irma geometry AND Irma impact

impact_years_wind_saba <- impact_years_wind |>
  filter(location == "Saba") |>
  pull(sim_year)

n_track        <- length(irma_track_years_saba)
n_impact       <- length(impact_years_wind_saba)
n_both         <- length(intersect(irma_track_years_saba, impact_years_wind_saba))
n_impact_only  <- length(setdiff(impact_years_wind_saba, irma_track_years_saba))

cat("\n--- Saba stress-test set comparison ---\n")
cat(sprintf(
  "  Track-based years    : %d  (%.1f%% of ensemble)\n",
  n_track, 100 * n_track / simulation_years
))
cat(sprintf(
  "  Impact-based years   : %d  (%.1f%% of ensemble)\n",
  n_impact, 100 * n_impact / simulation_years
))
cat(sprintf(
  "  Overlap (both sets)  : %d  — Irma geometry AND Irma-level wind\n",
  n_both
))
cat(sprintf(
  "  Impact-only years    : %d  — other storms that matched or exceeded Irma\n",
  n_impact_only
))

# Peak-wind distribution within each selection set
peak_by_set <- bind_rows(
  daily_out$Saba |>
    filter(sim_year %in% irma_track_years_saba) |>
    group_by(sim_year) |>
    summarise(peak_wind_kt = max(wind_kt, na.rm = TRUE), .groups = "drop") |>
    mutate(set = "Track-based (Irma geometry)"),
  daily_out$Saba |>
    filter(sim_year %in% impact_years_wind_saba) |>
    group_by(sim_year) |>
    summarise(peak_wind_kt = max(wind_kt, na.rm = TRUE), .groups = "drop") |>
    mutate(set = "Impact-based (\u2265 Irma wind)")
)

wind_summary <- peak_by_set |>
  group_by(set) |>
  summarise(
    n         = n(),
    p5_kt     = quantile(peak_wind_kt, 0.05, na.rm = TRUE),
    median_kt = median(peak_wind_kt,         na.rm = TRUE),
    p95_kt    = quantile(peak_wind_kt, 0.95, na.rm = TRUE),
    max_kt    = max(peak_wind_kt,            na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 1)))

cat("\n--- Peak-wind distribution by selection set (Saba) ---\n")
print(wind_summary)

# %%

# =============================================================================
# 8) Optional: visualise stress-test candidates
# =============================================================================

# Distribution of annual peak winds — full ensemble vs each stress-test set.
# Useful for communicating the severity tail to stakeholders.

annual_peak_all <- daily_out$Saba |>
  group_by(sim_year) |>
  summarise(peak_wind_kt = max(wind_kt, na.rm = TRUE), .groups = "drop") |>
  mutate(set = "Full ensemble")

p_stress <- bind_rows(annual_peak_all, peak_by_set) |>
  mutate(set = factor(set, levels = c(
    "Full ensemble",
    "Track-based (Irma geometry)",
    "Impact-based (\u2265 Irma wind)"
  ))) |>
  ggplot(aes(x = peak_wind_kt, fill = set, colour = set)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  scale_fill_manual(
    values = c(
      "Full ensemble"                     = "#4575b4",
      "Track-based (Irma geometry)"       = "#d73027",
      "Impact-based (\u2265 Irma wind)"   = "#f46d43"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Full ensemble"                     = "#4575b4",
      "Track-based (Irma geometry)"       = "#d73027",
      "Impact-based (\u2265 Irma wind)"   = "#f46d43"
    )
  ) +
  labs(
    x      = "Annual peak sustained wind (kt)",
    y      = "Density",
    title  = "Saba \u2014 Annual peak wind: full ensemble vs stress-test selections",
    fill   = NULL,
    colour = NULL
  ) +
  theme_light(base_size = 11) +
  theme(legend.position = "bottom")

print(p_stress)

# %%

# =============================================================================
# 9) Representative subset selection — portfolio level
# =============================================================================
#
# The impact-based candidate set may contain hundreds of years. For downstream
# deterministic modelling (e.g. running a loss model or infrastructure
# simulation) we typically need a much smaller set — say 10 to 20 years — that
# covers the full range of severity across the entire location portfolio.
#
# Workflow (three steps):
#   1) compute_stress_year_metrics()   — 7-metric summary per (location, sim_year)
#   2) aggregate_stress_metrics()      — aggregate across locations then score
#                                        → one row per sim_year (no location col)
#   3) select_stress_years()           — pick k years from the single pool

# --- Step 1: compute per-year metrics for the impact-based candidate set -----

# Pass the query-output tibble directly as sim_years; per-location filtering
# is applied automatically inside compute_stress_year_metrics().
stress_metrics <- compute_stress_year_metrics(
  daily     = daily_out,
  sim_years = impact_years_wind   # tibble with (location, sim_year) columns
)

cat("\n--- Stress-year metrics (first 6 rows) ---\n")
print(head(stress_metrics))

# --- Step 2: aggregate across locations and compute composite score ----------

# Default: aggregate by mean across locations, all seven metrics, uniform weights.
# Output has one row per sim_year and no location column (portfolio-level table).
stress_scored <- aggregate_stress_metrics(
  stress_metrics,
  location_agg = "mean"     # portfolio average across all five target locations
)

# Custom weighting example: up-weight intensity and cumulative damage,
# include compound-event count, ignore duration metrics.
# location_agg = "mean" is the default; shown explicitly for clarity.
stress_scored_custom <- aggregate_stress_metrics(
  stress_metrics,
  metrics_used = c("peak_wind_kt", "cum_damage", "max_damage_rate", "n_events"),
  weights      = c(
    peak_wind_kt    = 3,   # intensity matters most
    cum_damage      = 2,   # total damage second
    max_damage_rate = 1,   # acute shock
    n_events        = 1    # compound-event penalty
  ),
  location_agg = "mean"
)

cat("\n--- Composite score summary (custom weights, portfolio) ---\n")
stress_scored_custom |>
  summarise(
    n          = n(),
    score_min  = round(min(composite_score),    3),
    score_med  = round(median(composite_score), 3),
    score_max  = round(max(composite_score),    3)
  ) |>
  print()

# --- Step 3: select k representative years from the single portfolio pool ----

# Stratified (default): one year per severity bin — ensures mild, moderate,
# and extreme years are all represented.
selected_stratified <- select_stress_years(
  stress_scored_custom,
  k      = 10,
  method = "stratified"
)

# Diverse: spread across multi-dimensional metric space; useful when metrics
# are partially uncorrelated (e.g. long-duration low-wind vs short-duration
# high-wind events).
selected_diverse <- select_stress_years(
  stress_scored_custom,
  k      = 10,
  method = "diverse"
)

# Top-k: highest composite score only; use when the goal is a purely
# high-severity set (e.g. for PML estimation).
selected_top <- select_stress_years(
  stress_scored_custom,
  k      = 10,
  method = "top"
)

# The output has one row per selected sim_year — no location grouping.
cat("\n--- Selected stress-test years: portfolio (stratified, k = 10) ---\n")
selected_stratified |>
  select(sim_year, selection_rank, peak_wind_kt, cum_damage, n_events, composite_score) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  print()

# Extract the daily series for the final stress-test portfolio (all locations)
selected_years <- selected_stratified$sim_year

stress_portfolio <- lapply(daily_out, function(loc_tbl) {
  dplyr::filter(loc_tbl, sim_year %in% selected_years)
})

# Single-location convenience extract (e.g. for loss modelling at Saba)
stress_portfolio_saba <- stress_portfolio$Saba
