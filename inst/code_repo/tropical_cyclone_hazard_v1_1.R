
##############################################################
# Tropical Cyclone Hazard Assessment Model
# Using IBTrACS Data with Directional Wind Radii
##############################################################
#
# WORKFLOW OVERVIEW:
# 
# This script implements a stochastic hazard model for tropical cyclones
# at user-defined locations using the IBTrACS North Atlantic best-track 
# dataset. The key improvement over basic approaches is the use of 
# directional wind radii rather than assuming symmetric wind fields.
#
# For each target location (e.g., Miami, St. Martin), the workflow:
#
#   1. Extracts all storm track points (6-hourly observations) that pass 
#      within a specified buffer distance (default 200 km) of the location.
#
#   2. For each observation, calculates the bearing from the storm center 
#      to the target location and determines which meteorological quadrant 
#      (NE/SE/SW/NW) the location occupies relative to the storm.
#
#   3. Applies the appropriate directional wind radius (R34, R50, R64) for 
#      that specific quadrant, rather than assuming the maximum radius 
#      applies in all directions. This accounts for the asymmetric nature 
#      of tropical cyclone wind fields.
#
#   4. Aggregates individual track points to storm-level "events" with 
#      hazard metrics including minimum distance, maximum winds, exposure 
#      duration, and wind threshold exceedances.
#
#   5. Classifies each event into categorical severity classes based on 
#      actual wind exposure: none, Tropical Storm (34-63 kt), Category 1-2 
#      Hurricane (64+ kt), or Major Hurricane (Category 3+).
#
#   6. Counts how many events of each severity occur per year over a chosen 
#      reliable period (default: 1970-present, the satellite era).
#
#   7. Fits a Poisson model for each severity class to estimate the mean 
#      annual occurrence rate (lambda, λ) and the probability of at least 
#      one event occurring in any given year.
#
#   8. Uses Monte Carlo simulation to project the frequency and probability 
#      of events over a synthetic 1000-year period.
#
# The result is a statistically robust model of tropical cyclone occurrence 
# by severity at each location, suitable for risk assessment, planning, and 
# long-term hazard characterization. The directional radii approach provides 
# more accurate exposure estimates than methods that assume symmetric wind 
# fields, reducing false positives and improving rate estimates.
#
################################################################################

##############################################################
# Tropical Cyclone Hazard Model (Wind-threshold exceedance)
# At-island wind threshold exceedance (R34 / R64)
##############################################################

library(readr)
library(dplyr)
library(lubridate)
library(geosphere)
library(tidyr)
library(stringr)

set.seed(123)

################################################################################
# HELPERS
################################################################################

to_num_quiet <- function(x) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
  x <- str_replace_all(x, regex("\\b(degrees?_north|degrees?_south|degrees?_east|degrees?_west)\\b", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("\\b(deg|degrees?)\\b", ignore_case = TRUE), "")
  x <- str_trim(x)
  out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
  out[out %in% c(-999, -99, 999)] <- NA_real_
  out
}

max_na <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else max(x) }
min_na <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else min(x) }
mean_na <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else mean(x) }

# Great-circle distance storm center -> target (km)
dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}

# Bearing storm center -> target (degrees; geosphere may return -180..180)
calculate_bearing <- function(lat, lon, t_lat, t_lon) {
  geosphere::bearing(cbind(lon, lat), cbind(t_lon, t_lat))
}

# Quadrant from bearing (0..360)
get_quadrant <- function(bearing) {
  b <- (bearing + 360) %% 360
  dplyr::case_when(
    b >= 0   & b < 90  ~ "NE",
    b >= 90  & b < 180 ~ "SE",
    b >= 180 & b < 270 ~ "SW",
    b >= 270 & b < 360 ~ "NW",
    TRUE ~ NA_character_
  )
}

# Select directional radius value (nautical miles)
get_directional_radius <- function(quadrant, r_ne, r_se, r_sw, r_nw) {
  dplyr::case_when(
    quadrant == "NE" ~ r_ne,
    quadrant == "SE" ~ r_se,
    quadrant == "SW" ~ r_sw,
    quadrant == "NW" ~ r_nw,
    TRUE ~ NA_real_
  )
}

# Large computational gate; NOT part of impact definition
gate_km <- 800

# Conservative gating radius (km): include storms that could plausibly affect site
# Uses max R34 directional radius + margin. If radii missing, falls back to fallback_km.
compute_gate_km <- function(df, margin_km = 50, fallback_km = 800) {
  r34_max_nm <- suppressWarnings(max(df$USA_R34_NE, df$USA_R34_SE, df$USA_R34_SW, df$USA_R34_NW, na.rm = TRUE))
  if (!is.finite(r34_max_nm)) return(fallback_km)
  r34_max_nm * 1.852 + margin_km
}


to_radii_nm <- function(x, cap_nm) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
  
  out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
  
  # Common missing/sentinel patterns (extend aggressively)
  out[out %in% c(-999, -99, 0, 99, 999, 9999, 99999)] <- NA_real_
  
  # Anything implausibly large is missing/bad
  out[is.finite(out) & out > cap_nm] <- NA_real_
  
  out
}


make_storm_events <- function(df) {
  
  df <- df %>% distinct(SID, ISO_TIME, .keep_all = TRUE)
  
  df %>%
    group_by(SID, SEASON) %>%
    summarise(
      year = as.integer(first(SEASON)),
      name = dplyr::first(na.omit(NAME)),
      
      min_dist_km = min_na(dist_km),
      
      # Continuous at-site wind (kt)
      V_site_max_kt = max_na(V_site_kt),
      V_site_p95_kt = if (all(is.na(V_site_kt))) NA_real_ else
        as.numeric(stats::quantile(V_site_kt, probs = 0.95, na.rm = TRUE, names = FALSE)),
      
      # Exposure duration from continuous wind (hours)
      hours_ge34 = sum(V_site_kt >= 34, na.rm = TRUE) * 6,
      hours_ge50 = sum(V_site_kt >= 50, na.rm = TRUE) * 6,
      hours_ge64 = sum(V_site_kt >= 64, na.rm = TRUE) * 6,
      
      # Impact flags from continuous wind
      exposed_to_34kt = any(V_site_kt >= 34, na.rm = TRUE),
      exposed_to_50kt = any(V_site_kt >= 50, na.rm = TRUE),
      exposed_to_64kt = any(V_site_kt >= 64, na.rm = TRUE),
      
      # QA
      any_Vsite_known = any(!is.na(V_site_kt)),
      frac_Vsite_unknown = mean(is.na(V_site_kt)),
      
      # Metadata
      storm_speed_mean = mean_na(STORM_SPEED),
      max_center_wind  = max_na(USA_WIND),
      
      .groups = "drop"
    ) %>%
    filter(!is.na(year))
}


estimate_site_wind <- function(Vmax, r_km, R34_km, R50_km, R64_km, r0_mult = 1.5) {
  # Returns estimated sustained wind at site in knots (numeric scalar).
  # Robust to missing radii: falls back to the best available anchor.
  
  if (!is.finite(Vmax) || !is.finite(r_km) || r_km < 0) return(NA_real_)
  
  # Collect available anchors (radius -> wind)
  # Only keep finite, positive radii
  anchors_r <- c(R64_km, R50_km, R34_km)
  anchors_v <- c(64,     50,     34)
  
  ok <- is.finite(anchors_r) & anchors_r > 0
  anchors_r <- anchors_r[ok]
  anchors_v <- anchors_v[ok]
  
  # If no radii anchors, cannot estimate
  if (length(anchors_r) == 0L) return(NA_real_)
  
  # Ensure anchors are ordered from inner to outer radius
  o <- order(anchors_r)
  anchors_r <- anchors_r[o]
  anchors_v <- anchors_v[o]
  
  # Outer decay to zero at R0 = r0_mult * outermost anchor (typically R34)
  R_outer <- max(anchors_r)
  R0 <- r0_mult * R_outer
  
  # If beyond R0, wind ~ 0
  if (r_km >= R0) return(0)
  
  # If inside innermost anchor, interpolate between (0, Vmax) and (R_in, V_in)
  R_in <- anchors_r[1]
  V_in <- anchors_v[1]
  if (r_km <= R_in) {
    # linear decay from center to innermost anchor
    # cap at Vmax (cannot exceed)
    v <- Vmax + (V_in - Vmax) * (r_km / R_in)
    return(max(0, min(Vmax, v)))
  }
  
  # Between anchors: piecewise linear
  if (length(anchors_r) >= 2L) {
    for (j in 1:(length(anchors_r) - 1L)) {
      r1 <- anchors_r[j];   v1 <- anchors_v[j]
      r2 <- anchors_r[j+1]; v2 <- anchors_v[j+1]
      if (r_km > r1 && r_km <= r2) {
        v <- v1 + (v2 - v1) * ((r_km - r1) / (r2 - r1))
        return(max(0, min(Vmax, v)))
      }
    }
  }
  
  # Outside outermost anchor but within R0: decay from (R_outer, V_outer) to (R0, 0)
  V_outer <- anchors_v[length(anchors_v)]
  v <- V_outer + (0 - V_outer) * ((r_km - R_outer) / (R0 - R_outer))
  max(0, min(Vmax, v))
}

################################################################################
# 1. Read IBTrACS
################################################################################

ibtracs <- read_csv("./data/ibtracs.NA.list.v04r01.csv")

################################################################################
# 2. Target locations
################################################################################

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

################################################################################
# 3. Select and clean required IBTrACS variables
################################################################################

ib_sub <- ibtracs %>%
  select(
    SID, SEASON, BASIN, SUBBASIN, NAME, ISO_TIME,
    LAT, LON,
    USA_WIND, WMO_WIND,
    starts_with("USA_R34"), starts_with("USA_R50"), starts_with("USA_R64"),
    STORM_SPEED
  ) %>%
  mutate(
    LAT = to_num_quiet(LAT),
    LON = to_num_quiet(LON),
    LON = if_else(!is.na(LON) & LON > 180, LON - 360, LON),
    
    USA_WIND = to_num_quiet(USA_WIND),
    WMO_WIND = to_num_quiet(WMO_WIND),
    STORM_SPEED = to_num_quiet(STORM_SPEED),
    
    ISO_TIME = parse_date_time(
      ISO_TIME,
      orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S"),
      tz = "UTC"
    ),
    
    across(starts_with("USA_R34"), ~to_radii_nm(., cap_nm = 600)),
    across(starts_with("USA_R50"), ~to_radii_nm(., cap_nm = 400)),
    across(starts_with("USA_R64"), ~to_radii_nm(., cap_nm = 250))
  ) %>%
  filter(!is.na(LAT), !is.na(LON))


################################################################################
# 4. Track-point exposure computation PER LOCATION
################################################################################

results <- list()

for (i in seq_len(nrow(targets))) {
  loc <- targets[i, ]
  
  dat_loc <- ib_sub %>%
    mutate(
      dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon),
      bearing_to_target = calculate_bearing(LAT, LON, loc$lat, loc$lon),
      quadrant = get_quadrant(bearing_to_target),
      
      R34_nm = get_directional_radius(quadrant, USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
      R50_nm = get_directional_radius(quadrant, USA_R50_NE, USA_R50_SE, USA_R50_SW, USA_R50_NW),
      R64_nm = get_directional_radius(quadrant, USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW),
      
      R34_km = R34_nm * 1.852,
      R50_km = R50_nm * 1.852,
      R64_km = R64_nm * 1.852,
      # enforce monotone radii where possible
      R64_km = if_else(is.finite(R64_km) & is.finite(R50_km), pmin(R64_km, R50_km), R64_km),
      R50_km = if_else(is.finite(R50_km) & is.finite(R34_km), pmin(R50_km, R34_km), R50_km),

      Vmax_kt = USA_WIND,
      
      V_site_kt = mapply(
        estimate_site_wind,
        Vmax = Vmax_kt,
        r_km = dist_km,
        R34_km = R34_km,
        R50_km = R50_km,
        R64_km = R64_km,
        MoreArgs = list(r0_mult = 1.5)
      )
    ) %>%
    filter(dist_km <= gate_km)
  
  
  results[[loc$name]] <- dat_loc
}

################################################################################
# 5. Storm-event aggregation (ONE ROW PER SID PER ISLAND)
################################################################################

hazard_events <- lapply(results, make_storm_events)

################################################################################
# 6. Event severity (CONSISTENT: exposure-only)
################################################################################

hazard_events <- lapply(hazard_events, function(df) {
  df %>%
    mutate(
      severity = case_when(
        is.na(V_site_max_kt) ~ "unknown",
        V_site_max_kt < 34   ~ "none",
        V_site_max_kt < 64   ~ "TS",
        TRUE                 ~ "HUR64plus"
      )
    )
})

################################################################################
# 7. Location selection + STORM-EVENT annual counts (FULL FIX)
################################################################################

loc_name <- "Saba"  # "Statia", etc.

haz_loc <- hazard_events[[loc_name]] %>%
  filter(year >= 1970) %>%
  filter(any_Vsite_known) %>%      # optional but recommended
  filter(severity != "unknown")

# This MUST be 1 now.
haz_loc_counts <- haz_loc %>%
  filter(severity %in% c("TS", "HUR64plus")) %>%
  distinct(year, severity, SID) %>%
  count(year, severity, name = "n_events") %>%
  tidyr::complete(
    year = full_seq(range(year), 1),
    severity = c("TS", "HUR64plus"),
    fill = list(n_events = 0)
  )

# ---- Build a total-count series per year (two severities) ----
annual_total <- haz_loc_counts %>%
  group_by(year) %>%
  summarise(N = sum(n_events), .groups = "drop")

mu <- mean(annual_total$N)
va <- stats::var(annual_total$N)

# If va <= mu, no evidence of overdispersion; fall back to large k (≈ Poisson)

#small k_hat (~1–5) → strong year-to-year clustering
#large k_hat (>>50) → close to independent Poisson

k_hat <- if (is.finite(va) && va > mu && mu > 0) mu^2 / (va - mu) else 1e6
k_hat

# ---- Diagnostic: confirm one row per storm-event ----
# If this ever shows n > 1, your pipeline is not at event level.
diag_rows_per_sid <- haz_loc %>%
  count(year, SID) %>%
  summarise(max_rows_per_sid = max(n), .groups = "drop")
print(diag_rows_per_sid)
 
chk1 <- results[["Saba"]] %>%
  summarise(
    max_ratio = max(V_site_kt / Vmax_kt, na.rm = TRUE),
    max_diff  = max(V_site_kt - Vmax_kt, na.rm = TRUE)
  )

print(chk1) # max_ratio <= 1, max_diff <= 0 (allow tiny numerical jitter ~1e-8)

chk2 <- results[["Saba"]] %>%
  filter(is.finite(R34_km), is.finite(R50_km), is.finite(R64_km)) %>%
  summarise(
    frac_bad_order = mean(!(R64_km <= R50_km & R50_km <= R34_km)),
    n_checked = n()
  )

print(chk2) 
# If frac_bad_order is large (>5–10%), you need to enforce monotonicity per point (simple fix below).

chk3 <- results[["Saba"]] %>%
  summarise(
    frac_ge34 = mean(V_site_kt >= 34, na.rm = TRUE),
    frac_ge64 = mean(V_site_kt >= 64, na.rm = TRUE)
  )
# “What fraction of 6-hour points experienced TS / hurricane winds at the site?”
print(chk3)

################################################################################
# 8. Estimate Poisson λ (baseline)
################################################################################

lambda_table <- haz_loc_counts %>%
  group_by(severity) %>%
  summarise(
    lambda = mean(n_events),
    n_years = n(),
    p_at_least_one = 1 - exp(-lambda),
    p_zero         = exp(-lambda),
    .groups = "drop"
  )

print(lambda_table)

################################################################################
# 9. Monte Carlo simulation
################################################################################

#Quantity	Your λ	Typical observed range	Verdict
#TS (34+ kt)	0.88 / yr	~0.5 – 1.5 / yr	✔ Reasonable
#HUR64+	0.059 / yr	~0.05 – 0.3 / yr	✔ Reasonable

n_years_sim <- 1000

# Baseline lambdas
lam_TS  <- lambda_table$lambda[lambda_table$severity == "TS"]
lam_HUR <- lambda_table$lambda[lambda_table$severity == "HUR64plus"]

# Shared annual activity factor A_y ~ Gamma(k, k)  (mean=1, var=1/k)
A <- rgamma(n_years_sim, shape = k_hat, rate = k_hat)

sim_twolevel <- tibble::tibble(
  sim_year = 1:n_years_sim,
  A = A,
  n_TS = rpois(n_years_sim, lam_TS  * A),
  n_HUR64plus = rpois(n_years_sim, lam_HUR * A)
)

# Any-event probability
p_any_event <- sim_twolevel %>%
  summarise(p_at_least_one = mean((n_TS + n_HUR64plus) >= 1))

p_any_event

#compute implied correlation (sanity check)
cor(sim_twolevel$n_TS, sim_twolevel$n_HUR64plus)


annual_total %>% summarise(mean=mean(N), var=var(N))
#If var ≈ mean, the two-level model will behave almost like Poisson (fine).
#### Double check/validation 


out <- haz_loc %>%
  summarise(
    years = n_distinct(year),
    n_TS = sum(severity == "TS"),
    n_HUR = sum(severity == "HUR64plus")
  )

out$n_TS / out$years #≈ 0.8–1.0
out$n_HUR / out$years #≈ 0.05–0.1


