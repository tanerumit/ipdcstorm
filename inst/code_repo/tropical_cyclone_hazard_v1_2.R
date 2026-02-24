################################################################################
# Tropical Cyclone Hazard Assessment Model (IBTrACS + Directional Radii)
# - Impact definition: at-island sustained wind (V_site_kt) and thresholds
# - Event definition: one row per storm (SID, SEASON) per island
# - Rate model: 2-level Poisson with shared annual activity (Gamma-Poisson)
################################################################################

library(readr)
library(dplyr)
library(lubridate)
library(geosphere)
library(tidyr)
library(stringr)


set.seed(123)

################################################################################
# 0. CONFIGURATION
################################################################################

cfg <- list(
  # Data
  ibtracs_path = "./data/ibtracs.NA.list.v04r01.csv",
  min_year = 1970,
  
  # Computational gate (NOT impact definition)
  gate_km = 800,
  
  # Wind thresholds (kt)
  thr_ts = 34,
  thr_hur = 64,
  thr_50 = 50,
  
  # Radii parsing caps (nautical miles)
  cap_r34_nm = 600,
  cap_r50_nm = 400,
  cap_r64_nm = 250,
  
  # Wind field model selection
  use_asymmetry = TRUE,          # NEW: Add forward motion asymmetry
  use_holland = TRUE,
  
  # Piecewise decay tail
  r0_mult = 1.5
)

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

################################################################################
# 1. HELPERS (I/O cleaning, geometry, radii selection)
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

max_na  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else max(x) }
min_na  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else min(x) }
mean_na <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else mean(x) }

to_radii_nm <- function(x, cap_nm) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
  out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
  out[out %in% c(-999, -99, 0, 99, 999, 9999, 99999)] <- NA_real_
  out[is.finite(out) & out > cap_nm] <- NA_real_
  out
}

dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}

calculate_bearing <- function(lat, lon, t_lat, t_lon) {
  geosphere::bearing(cbind(lon, lat), cbind(t_lon, t_lat))
}

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

get_directional_radius <- function(quadrant, r_ne, r_se, r_sw, r_nw) {
  dplyr::case_when(
    quadrant == "NE" ~ r_ne,
    quadrant == "SE" ~ r_se,
    quadrant == "SW" ~ r_sw,
    quadrant == "NW" ~ r_nw,
    TRUE ~ NA_real_
  )
}

enforce_monotone_radii <- function(R34_km, R50_km, R64_km) {
  # enforce R64 <= R50 <= R34 where finite
  R64_km <- if_else(is.finite(R64_km) & is.finite(R50_km), pmin(R64_km, R50_km), R64_km)
  R50_km <- if_else(is.finite(R50_km) & is.finite(R34_km), pmin(R50_km, R34_km), R50_km)
  list(R34_km = R34_km, R50_km = R50_km, R64_km = R64_km)
}

################################################################################
# 2. SITE WIND ESTIMATOR (piecewise decay anchored to R64/R50/R34)
################################################################################

estimate_site_wind <- function(Vmax, r_km, R34_km, R50_km, R64_km, r0_mult = 1.5) {
  if (!is.finite(Vmax) || !is.finite(r_km) || r_km < 0) return(NA_real_)
  
  anchors_r <- c(R64_km, R50_km, R34_km)
  anchors_v <- c(64,     50,     34)
  
  ok <- is.finite(anchors_r) & anchors_r > 0
  anchors_r <- anchors_r[ok]
  anchors_v <- anchors_v[ok]
  if (length(anchors_r) == 0L) return(NA_real_)
  
  o <- order(anchors_r)
  anchors_r <- anchors_r[o]
  anchors_v <- anchors_v[o]
  
  R_outer <- max(anchors_r)
  R0 <- r0_mult * R_outer
  if (r_km >= R0) return(0)
  
  R_in <- anchors_r[1]
  V_in <- anchors_v[1]
  if (r_km <= R_in) {
    v <- Vmax + (V_in - Vmax) * (r_km / R_in)
    return(max(0, min(Vmax, v)))
  }
  
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
  
  V_outer <- anchors_v[length(anchors_v)]
  v <- V_outer + (0 - V_outer) * ((r_km - R_outer) / (R0 - R_outer))
  max(0, min(Vmax, v))
}

estimate_site_wind_holland <- function(
    Vmax_kt,           
    r_km,              
    R34_km,            
    R50_km = NA,       
    R64_km = NA,       
    RMW_km = NA,       
    Pn = 1013,         
    Pc = NA            
) {
  
  # Input validation
  if (!is.finite(Vmax_kt) || Vmax_kt <= 0) return(NA_real_)
  if (!is.finite(r_km) || r_km < 0) return(NA_real_)
  
  # ---- Step 1: Estimate Radius of Maximum Winds (RMW) ----
  if (!is.finite(RMW_km)) {
    if (is.finite(R64_km) && R64_km > 0) {
      RMW_km <- 0.35 * R64_km  
    } else if (is.finite(R50_km) && R50_km > 0) {
      RMW_km <- 0.40 * R50_km  
    } else if (is.finite(R34_km) && R34_km > 0) {
      RMW_km <- 0.50 * R34_km  
    } else {
      # Knaff et al. (2015) climatology
      RMW_nm <- 35 - 0.11 * Vmax_kt + 0.00043 * Vmax_kt^2
      RMW_km <- pmax(10, pmin(150, RMW_nm * 1.852))  
    }
  }
  
  RMW_km <- pmax(5, pmin(200, RMW_km))
  
  # ---- Step 2: Estimate Holland B parameter ----
  if (is.finite(Pc) && is.finite(Pn) && Pc < Pn) {
    deltaP <- Pn - Pc
    B <- -4.4e-5 * deltaP^2 + 0.01 * deltaP + 0.03 * (deltaP - 25)
    B <- pmax(1.0, pmin(2.5, B))  
  } else {
    if (Vmax_kt >= 100) {
      B <- 1.2 + (120 - Vmax_kt) * 0.005  
    } else if (Vmax_kt >= 64) {
      B <- 1.4 + (100 - Vmax_kt) * 0.008  
    } else {
      B <- 1.8 + (64 - Vmax_kt) * 0.01    
    }
    B <- pmax(0.9, pmin(2.5, B))
  }
  
  # ---- Step 3: Apply Holland Model ----
  if (r_km < 0.1) {
    return(0)
  }
  
  r_norm <- r_km / RMW_km
  
  if (r_norm < 1.0) {
    V_gradient <- sqrt((1/r_norm)^B * exp(1 - (1/r_norm)^B))
    V_site_kt <- Vmax_kt * V_gradient
  } else {
    V_gradient <- sqrt((RMW_km/r_km)^B * exp(1 - (RMW_km/r_km)^B))
    V_site_kt <- Vmax_kt * V_gradient
  }
  
  # ---- Step 4: Calibration adjustment using R34 ----
  if (is.finite(R34_km) && R34_km > 0 && r_km > RMW_km * 1.5) {
    r_norm_34 <- R34_km / RMW_km
    V_at_R34_model <- Vmax_kt * sqrt((RMW_km/R34_km)^B * exp(1 - (RMW_km/R34_km)^B))
    
    if (V_at_R34_model > 25 && V_at_R34_model < 50) {  
      correction <- 34 / V_at_R34_model
      V_site_kt <- V_site_kt * correction
    }
  }
  
  # ---- Step 5: Outer decay cutoff ----
  R_outer <- if (is.finite(R34_km) && R34_km > 0) 1.8 * R34_km else 300
  
  if (r_km > R_outer) {
    decay_factor <- exp(-2 * (r_km - R_outer) / R_outer)
    V_site_kt <- V_site_kt * decay_factor
  }
  
  # Final bounds
  V_site_kt <- pmax(0, pmin(Vmax_kt, V_site_kt))
  
  return(V_site_kt)
}

################################################################################
# IMPROVED WIND FIELD FUNCTIONS - ADD AFTER EXISTING HELPERS
################################################################################


add_forward_motion_asymmetry <- function(
    V_site_base_kt,      # Wind from Holland model (symmetric)
    storm_speed_kt,      # Forward speed of storm (kt)
    r_km,                # Distance from center (km)
    bearing_to_target,   # Bearing from storm center to target (degrees, 0-360)
    storm_heading,       # Direction storm is moving (degrees, 0-360)
    RMW_km = 40         # Radius of maximum winds (km)
) {
  
  # Input validation
  if (!is.finite(V_site_base_kt) || !is.finite(storm_speed_kt)) {
    return(V_site_base_kt)
  }
  if (!is.finite(bearing_to_target) || !is.finite(storm_heading)) {
    return(V_site_base_kt)
  }
  if (storm_speed_kt < 0.1) {
    # Stationary storm - no asymmetry
    return(V_site_base_kt)
  }
  
  # ---- Step 1: Calculate relative angle ----
  # θ = angle between storm motion vector and vector to target
  # θ = 0° → target directly ahead of storm (maximum enhancement)
  # θ = 90° → target to the right (moderate enhancement, NH)
  # θ = 180° → target directly behind storm (maximum reduction)
  
  angle_diff <- bearing_to_target - storm_heading
  # Normalize to -180 to +180
  angle_diff <- ((angle_diff + 180) %% 360) - 180
  # Convert to radians
  theta_rad <- angle_diff * pi / 180
  
  
  # ---- Step 2: Calculate coupling coefficient K(r) ----
  # K varies with distance from storm center:
  # - Maximum near RMW (~0.6)
  # - Decreases linearly to ~0.4 at R34
  # - Further decreases beyond R34
  
  r_norm <- r_km / RMW_km
  
  if (r_norm < 1.0) {
    # Inside RMW: K increases toward RMW
    K <- 0.4 + 0.25 * r_norm  # K: 0.4 at center → 0.65 at RMW
  } else if (r_norm <= 3.0) {
    # RMW to ~3*RMW: K decreases gradually
    K <- 0.65 - 0.15 * (r_norm - 1.0) / 2.0  # K: 0.65 at RMW → 0.50 at 3*RMW
  } else {
    # Far field: K continues to decrease
    K <- 0.50 * exp(-0.2 * (r_norm - 3.0))  # Exponential decay beyond 3*RMW
  }
  
  K <- pmax(0.1, pmin(0.7, K))  # Bounds: 0.1 - 0.7
  
  
  # ---- Step 3: Calculate motion-induced wind component ----
  # This is the vector component of storm motion that adds/subtracts to rotation
  V_motion_kt <- K * storm_speed_kt * cos(theta_rad)
  
  
  # ---- Step 4: Add to symmetric wind field ----
  V_site_asymmetric_kt <- V_site_base_kt + V_motion_kt
  
  
  # ---- Step 5: Physical bounds ----
  # Wind cannot be negative (stalled region behind very fast storms is rare)
  # Wind cannot exceed reasonable maximum (Vmax + full translation speed)
  V_site_asymmetric_kt <- pmax(0, V_site_asymmetric_kt)
  
  return(V_site_asymmetric_kt)
}


# Compute storm heading (needed for asymmetry)
compute_storm_heading <- function(df) {
  df <- df %>%
    arrange(SID, ISO_TIME) %>%
    group_by(SID) %>%
    mutate(
      LAT_next = lead(LAT),
      LON_next = lead(LON),
      
      storm_heading = geosphere::bearing(
        cbind(LON, LAT), 
        cbind(LON_next, LAT_next)
      ),
      
      storm_heading = (storm_heading + 360) %% 360,
      
      storm_heading = ifelse(is.na(storm_heading), 
                             lag(storm_heading), 
                             storm_heading)
    ) %>%
    ungroup() %>%
    select(-LAT_next, -LON_next)
  
  return(df)
}


################################################################################
# FORWARD MOTION ASYMMETRY CORRECTION
################################################################################

add_forward_motion_asymmetry <- function(
    V_site_base_kt,      # Wind from Holland model (symmetric)
    storm_speed_kt,      # Forward speed of storm (kt)
    r_km,                # Distance from center (km)
    bearing_to_target,   # Bearing from storm center to target (degrees, 0-360)
    storm_heading,       # Direction storm is moving (degrees, 0-360)
    RMW_km = 40         # Radius of maximum winds (km)
) {
  
  # Input validation
  if (!is.finite(V_site_base_kt) || !is.finite(storm_speed_kt)) {
    return(V_site_base_kt)
  }
  if (!is.finite(bearing_to_target) || !is.finite(storm_heading)) {
    return(V_site_base_kt)
  }
  if (storm_speed_kt < 0.1) {
    # Stationary storm - no asymmetry
    return(V_site_base_kt)
  }
  
  # ---- Step 1: Calculate relative angle ----
  # θ = angle between storm motion vector and vector to target
  # θ = 0° → target directly ahead of storm (maximum enhancement)
  # θ = 90° → target to the right (moderate enhancement, NH)
  # θ = 180° → target directly behind storm (maximum reduction)
  
  angle_diff <- bearing_to_target - storm_heading
  # Normalize to -180 to +180
  angle_diff <- ((angle_diff + 180) %% 360) - 180
  # Convert to radians
  theta_rad <- angle_diff * pi / 180
  
  
  # ---- Step 2: Calculate coupling coefficient K(r) ----
  # K varies with distance from storm center:
  # - Maximum near RMW (~0.6)
  # - Decreases linearly to ~0.4 at R34
  # - Further decreases beyond R34
  
  r_norm <- r_km / RMW_km
  
  if (r_norm < 1.0) {
    # Inside RMW: K increases toward RMW
    K <- 0.4 + 0.25 * r_norm  # K: 0.4 at center → 0.65 at RMW
  } else if (r_norm <= 3.0) {
    # RMW to ~3*RMW: K decreases gradually
    K <- 0.65 - 0.15 * (r_norm - 1.0) / 2.0  # K: 0.65 at RMW → 0.50 at 3*RMW
  } else {
    # Far field: K continues to decrease
    K <- 0.50 * exp(-0.2 * (r_norm - 3.0))  # Exponential decay beyond 3*RMW
  }
  
  K <- pmax(0.1, pmin(0.7, K))  # Bounds: 0.1 - 0.7
  
  
  # ---- Step 3: Calculate motion-induced wind component ----
  # This is the vector component of storm motion that adds/subtracts to rotation
  V_motion_kt <- K * storm_speed_kt * cos(theta_rad)
  
  
  # ---- Step 4: Add to symmetric wind field ----
  V_site_asymmetric_kt <- V_site_base_kt + V_motion_kt
  
  
  # ---- Step 5: Physical bounds ----
  # Wind cannot be negative (stalled region behind very fast storms is rare)
  # Wind cannot exceed reasonable maximum (Vmax + full translation speed)
  V_site_asymmetric_kt <- pmax(0, V_site_asymmetric_kt)
  
  return(V_site_asymmetric_kt)
}


################################################################################
# COMPUTE STORM HEADING (bearing between consecutive track points)
################################################################################

compute_storm_heading <- function(df) {
  # df must have: SID, ISO_TIME, LAT, LON
  # Returns: storm_heading (degrees) for each point
  
  df <- df %>%
    arrange(SID, ISO_TIME) %>%
    group_by(SID) %>%
    mutate(
      # Lead to next point
      LAT_next = lead(LAT),
      LON_next = lead(LON),
      
      # Bearing from current point to next point
      storm_heading = geosphere::bearing(
        cbind(LON, LAT), 
        cbind(LON_next, LAT_next)
      ),
      
      # Normalize to 0-360
      storm_heading = (storm_heading + 360) %% 360,
      
      # For last point in track, use same heading as previous point
      storm_heading = ifelse(is.na(storm_heading), 
                             lag(storm_heading), 
                             storm_heading)
    ) %>%
    ungroup() %>%
    select(-LAT_next, -LON_next)
  
  return(df)
}


################################################################################
# INTEGRATED WORKFLOW: Holland + Asymmetry
################################################################################

compute_site_winds_full <- function(df, target_lat, target_lon) {
  # df should have: SID, ISO_TIME, LAT, LON, USA_WIND, USA_R34_*, USA_R50_*, USA_R64_*
  
  # Step 1: Geometric calculations
  df <- df %>%
    mutate(
      dist_km = if_else(is.na(dist_km), dist_to_target(LAT, LON, target_lat, target_lon), dist_km),
      bearing_to_target = calculate_bearing(LAT, LON, target_lat, target_lon),
      quadrant = get_quadrant(bearing_to_target),
      
      # Directional radii
      R34_nm = get_directional_radius(quadrant, USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
      R50_nm = get_directional_radius(quadrant, USA_R50_NE, USA_R50_SE, USA_R50_SW, USA_R50_NW),
      R64_nm = get_directional_radius(quadrant, USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW),
      
      R34_km = R34_nm * 1.852,
      R50_km = R50_nm * 1.852,
      R64_km = R64_nm * 1.852,
      
      Vmax_kt = USA_WIND
    )
  
  # Step 2: Compute storm headings
  df <- compute_storm_heading(df)
  
  # Step 3: Estimate RMW (if not in IBTrACS)
  df <- df %>%
    mutate(
      RMW_km = case_when(
        is.finite(R64_km) & R64_km > 0 ~ 0.35 * R64_km,
        is.finite(R50_km) & R50_km > 0 ~ 0.40 * R50_km,
        is.finite(R34_km) & R34_km > 0 ~ 0.50 * R34_km,
        TRUE ~ pmax(10, pmin(150, (35 - 0.11*Vmax_kt + 0.00043*Vmax_kt^2) * 1.852))
      )
    )
  
  # Step 4: Holland model (symmetric wind field)
  df$V_site_symmetric_kt <- mapply(
    estimate_site_wind_holland,
    Vmax_kt = df$Vmax_kt,
    r_km = df$dist_km,
    R34_km = df$R34_km,
    R50_km = df$R50_km,
    R64_km = df$R64_km,
    RMW_km = df$RMW_km,
    MoreArgs = list(Pn = 1013)
  )
  
  # Step 5: Add forward motion asymmetry
  df$V_site_kt <- mapply(
    add_forward_motion_asymmetry,
    V_site_base_kt = df$V_site_symmetric_kt,
    storm_speed_kt = df$STORM_SPEED,
    r_km = df$dist_km,
    bearing_to_target = df$bearing_to_target,
    storm_heading = df$storm_heading,
    RMW_km = df$RMW_km
  )
  
  return(df)
}



################################################################################
# 3. EVENT AGGREGATION (one row per storm per island)
################################################################################

make_storm_events <- function(df, thr_ts = 34, thr_50 = 50, thr_hur = 64) {
  
  df <- df %>% distinct(SID, ISO_TIME, .keep_all = TRUE)
  
  df %>%
    group_by(SID, SEASON) %>%
    summarise(
      year = as.integer(first(SEASON)),
      name = dplyr::first(na.omit(NAME)),
      min_dist_km = min_na(dist_km),
      
      V_site_max_kt = max_na(V_site_kt),
      V_site_p95_kt = if (all(is.na(V_site_kt))) NA_real_ else
        as.numeric(stats::quantile(V_site_kt, probs = 0.95, na.rm = TRUE, names = FALSE)),
      
      hours_ge34 = sum(V_site_kt >= thr_ts, na.rm = TRUE) * 6,
      hours_ge50 = sum(V_site_kt >= thr_50, na.rm = TRUE) * 6,
      hours_ge64 = sum(V_site_kt >= thr_hur, na.rm = TRUE) * 6,
      
      exposed_to_34kt = any(V_site_kt >= thr_ts, na.rm = TRUE),
      exposed_to_50kt = any(V_site_kt >= thr_50, na.rm = TRUE),
      exposed_to_64kt = any(V_site_kt >= thr_hur, na.rm = TRUE),
      
      any_Vsite_known = any(!is.na(V_site_kt)),
      frac_Vsite_unknown = mean(is.na(V_site_kt)),
      
      storm_speed_mean = mean_na(STORM_SPEED),
      max_center_wind  = max_na(USA_WIND),
      
      .groups = "drop"
    ) %>%
    filter(!is.na(year))
}

classify_severity <- function(V_site_max_kt, thr_ts = 34, thr_hur = 64) {
  dplyr::case_when(
    is.na(V_site_max_kt) ~ "unknown",
    V_site_max_kt < thr_ts ~ "none",
    V_site_max_kt < thr_hur ~ "TS",
    TRUE ~ "HUR"
  )
}

################################################################################
# 4. DATA INGEST + CLEANING
################################################################################

read_ibtracs_clean <- function(path, cfg) {
  readr::read_csv(path, show_col_types = FALSE) %>%
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
      
      across(starts_with("USA_R34"), ~to_radii_nm(., cap_nm = cfg$cap_r34_nm)),
      across(starts_with("USA_R50"), ~to_radii_nm(., cap_nm = cfg$cap_r50_nm)),
      across(starts_with("USA_R64"), ~to_radii_nm(., cap_nm = cfg$cap_r64_nm))
    ) %>%
    filter(!is.na(LAT), !is.na(LON))
}

################################################################################
# 5. PER-LOCATION PROCESSING (track points -> V_site_kt)
################################################################################

compute_location_trackpoints <- function(ib_sub, loc, cfg) {
  
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
      R64_km = R64_nm * 1.852
    ) %>%
    { 
      mm <- enforce_monotone_radii(.$R34_km, .$R50_km, .$R64_km)
      mutate(., R34_km = mm$R34_km, R50_km = mm$R50_km, R64_km = mm$R64_km)
    } %>%
    mutate(Vmax_kt = USA_WIND)
  
  # Compute storm headings if asymmetry is enabled
  if (cfg$use_asymmetry) {
    dat_loc <- compute_storm_heading(dat_loc)
  }
  
  # Estimate RMW for Holland model and asymmetry
  if (cfg$use_holland || cfg$use_asymmetry) {
    dat_loc <- dat_loc %>%
      mutate(
        RMW_km = case_when(
          is.finite(R64_km) & R64_km > 0 ~ 0.35 * R64_km,
          is.finite(R50_km) & R50_km > 0 ~ 0.40 * R50_km,
          is.finite(R34_km) & R34_km > 0 ~ 0.50 * R34_km,
          TRUE ~ pmax(10, pmin(150, (35 - 0.11*Vmax_kt + 0.00043*Vmax_kt^2) * 1.852))
        )
      )
  }
  
  # ---- WIND FIELD CALCULATION ----
  if (cfg$use_holland) {
    # Step 1: Holland symmetric wind field
    dat_loc$V_site_symmetric_kt <- mapply(
      estimate_site_wind_holland,
      Vmax_kt = dat_loc$Vmax_kt,
      r_km = dat_loc$dist_km,
      R34_km = dat_loc$R34_km,
      R50_km = dat_loc$R50_km,
      R64_km = dat_loc$R64_km,
      RMW_km = dat_loc$RMW_km,
      MoreArgs = list(Pn = 1013)
    )
    
    # Step 2: Add asymmetry if enabled
    if (cfg$use_asymmetry) {
      dat_loc$V_site_kt <- mapply(
        add_forward_motion_asymmetry,
        V_site_base_kt = dat_loc$V_site_symmetric_kt,
        storm_speed_kt = dat_loc$STORM_SPEED,
        r_km = dat_loc$dist_km,
        bearing_to_target = dat_loc$bearing_to_target,
        storm_heading = dat_loc$storm_heading,
        RMW_km = dat_loc$RMW_km
      )
    } else {
      dat_loc$V_site_kt <- dat_loc$V_site_symmetric_kt
    }
    
  } else {
    # Original piecewise linear model
    dat_loc$V_site_kt <- mapply(
      estimate_site_wind,
      Vmax = dat_loc$Vmax_kt,
      r_km = dat_loc$dist_km,
      R34_km = dat_loc$R34_km,
      R50_km = dat_loc$R50_km,
      R64_km = dat_loc$R64_km,
      MoreArgs = list(r0_mult = cfg$r0_mult)
    )
  }
  
  dat_loc %>% filter(dist_km <= cfg$gate_km)
}

################################################################################
# 6. RATES + TWO-LEVEL SIMULATION
################################################################################

compute_annual_counts <- function(haz_loc, severities = c("TS", "HUR")) {
  haz_loc %>%
    filter(severity %in% severities) %>%
    distinct(year, severity, SID) %>%
    count(year, severity, name = "n_events") %>%
    tidyr::complete(
      year = full_seq(range(year), 1),
      severity = severities,
      fill = list(n_events = 0)
    )
}

compute_lambda_table <- function(haz_loc_counts) {
  haz_loc_counts %>%
    group_by(severity) %>%
    summarise(
      lambda = mean(n_events),
      n_years = n(),
      p_at_least_one = 1 - exp(-lambda),
      p_zero = exp(-lambda),
      .groups = "drop"
    )
}

estimate_k_hat <- function(haz_loc_counts) {
  annual_total <- haz_loc_counts %>%
    group_by(year) %>%
    summarise(N = sum(n_events), .groups = "drop")
  
  mu <- mean(annual_total$N)
  va <- stats::var(annual_total$N)
  
  k_hat <- if (is.finite(va) && va > mu && mu > 0) mu^2 / (va - mu) else 1e6
  list(k_hat = k_hat, annual_total = annual_total, mu = mu, var = va)
}

simulate_twolevel_counts <- function(lambda_table, k_hat, n_years_sim = 1000) {
  lam_TS  <- lambda_table$lambda[lambda_table$severity == "TS"]
  lam_HUR <- lambda_table$lambda[lambda_table$severity == "HUR"]
  
  A <- rgamma(n_years_sim, shape = k_hat, rate = k_hat)
  
  tibble::tibble(
    sim_year = 1:n_years_sim,
    A = A,
    n_TS = rpois(n_years_sim, lam_TS * A),
    n_HUR = rpois(n_years_sim, lam_HUR * A)
  )
}

################################################################################
# 7. DIAGNOSTICS
################################################################################

diagnostics_trackpoints <- function(df) {
  list(
    wind_cap = df %>%
      summarise(
        max_ratio = max(V_site_kt / Vmax_kt, na.rm = TRUE),
        max_diff  = max(V_site_kt - Vmax_kt, na.rm = TRUE),
        .groups = "drop"
      ),
    radii_order = df %>%
      filter(is.finite(R34_km), is.finite(R50_km), is.finite(R64_km)) %>%
      summarise(
        frac_bad_order = mean(!(R64_km <= R50_km & R50_km <= R34_km)),
        n_checked = n(),
        .groups = "drop"
      ),
    exceedance = df %>%
      summarise(
        frac_ge34 = mean(V_site_kt >= 34, na.rm = TRUE),
        frac_ge64 = mean(V_site_kt >= 64, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

diagnostics_events <- function(haz_loc) {
  list(
    one_row_per_event = haz_loc %>%
      count(year, SID) %>%
      summarise(max_rows_per_sid = max(n), .groups = "drop"),
    storms_per_year = haz_loc %>%
      summarise(
        years = n_distinct(year),
        n_TS = sum(severity == "TS"),
        n_HUR = sum(severity == "HUR"),
        .groups = "drop"
      )
  )
}

################################################################################
# 8. RUN PIPELINE
################################################################################

# 8.1 Read & clean
ib_sub <- read_ibtracs_clean(cfg$ibtracs_path, cfg)

# 8.2 Compute per-location track points
results <- setNames(vector("list", nrow(targets)), targets$name)

for (i in seq_len(nrow(targets))) {
  loc <- targets[i, ]
  
  # Compute distance ONCE for gating
  dat_loc <- ib_sub %>%
    mutate(dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon)) %>%
    filter(dist_km <= cfg$gate_km)
  
  # Then compute the rest (bearing, quadrant, radii, wind model, etc.)
  dat_loc <- compute_site_winds_full(dat_loc, loc$lat, loc$lon)
  
  results[[loc$name]] <- dat_loc
}

# 8.3 Aggregate to storm events + classify severity
hazard_events <- lapply(results, function(df) {
  ev <- make_storm_events(df, thr_ts = cfg$thr_ts, thr_50 = cfg$thr_50, thr_hur = cfg$thr_hur)
  ev %>% mutate(severity = classify_severity(V_site_max_kt, thr_ts = cfg$thr_ts, thr_hur = cfg$thr_hur))
})

# 8.4 Select location
loc_name <- "Saba"
haz_loc <- hazard_events[[loc_name]] %>%
  filter(year >= cfg$min_year) %>%
  filter(any_Vsite_known) %>%
  filter(severity != "unknown")

# 8.5 Annual counts + rate model
haz_loc_counts <- compute_annual_counts(haz_loc, severities = c("TS", "HUR"))
lambda_table <- compute_lambda_table(haz_loc_counts)

k_info <- estimate_k_hat(haz_loc_counts)
k_hat <- k_info$k_hat

# 8.6 Simulation
sim_twolevel <- simulate_twolevel_counts(lambda_table, k_hat, n_years_sim = 1000)

p_any_event <- sim_twolevel %>%
  summarise(p_at_least_one = mean((n_TS + n_HUR) >= 1), .groups = "drop")

corr_TS_HUR <- cor(sim_twolevel$n_TS, sim_twolevel$n_HUR)

################################################################################
# 9. OUTPUTS (minimal prints; expand later to writing files)
################################################################################

print(lambda_table)
print(p_any_event)
print(tibble::tibble(k_hat = k_hat, annual_mean = k_info$mu, annual_var = k_info$var, corr_TS_HUR = corr_TS_HUR))

# Diagnostics
diag_tp <- diagnostics_trackpoints(results[[loc_name]])
diag_ev <- diagnostics_events(haz_loc)

print(diag_tp$wind_cap)
print(diag_tp$radii_order)
print(diag_tp$exceedance)
print(diag_ev$one_row_per_event)
print(diag_ev$storms_per_year)
