
# R/functions_hazard_enhanced.R
# Tropical Cyclone Hazard Model functions (IBTrACS + directional radii + Holland + asymmetry)
# ENHANCEMENTS:
#   - Refined Saffir-Simpson severity classification (5 classes)
#   - Integrated Wind Exposure Index (cumulative damage potential)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(geosphere)
  library(readr)
})

# =============================================================================
# 1) Low-level helpers
# =============================================================================

to_num_quiet <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
  x <- stringr::str_replace_all(
    x,
    stringr::regex("\\b(degrees?_north|degrees?_south|degrees?_east|degrees?_west)\\b", ignore_case = TRUE),
    ""
  )
  x <- stringr::str_replace_all(x, stringr::regex("\\b(deg|degrees?)\\b", ignore_case = TRUE), "")
  x <- stringr::str_trim(x)
  
  out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
  out[out %in% c(-999, -99, 999)] <- NA_real_
  out
}

max_na  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else max(x) }
min_na  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else min(x) }
mean_na <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else mean(x) }
sum_na  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) 0 else sum(x) }

to_radii_nm <- function(x, cap_nm) {
  x <- stringr::str_trim(as.character(x))
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
  R64_km <- if_else(is.finite(R64_km) & is.finite(R50_km), pmin(R64_km, R50_km), R64_km)
  R50_km <- if_else(is.finite(R50_km) & is.finite(R34_km), pmin(R50_km, R34_km), R50_km)
  list(R34_km = R34_km, R50_km = R50_km, R64_km = R64_km)
}

# =============================================================================
# 2) Wind field models
# =============================================================================

estimate_site_wind_piecewise <- function(Vmax, r_km, R34_km, R50_km, R64_km, r0_mult = 1.5) {
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
  if (!is.finite(Vmax_kt) || Vmax_kt <= 0) return(NA_real_)
  if (!is.finite(r_km) || r_km < 0) return(NA_real_)
  
  # RMW estimate
  if (!is.finite(RMW_km)) {
    if (is.finite(R64_km) && R64_km > 0) {
      RMW_km <- 0.35 * R64_km
    } else if (is.finite(R50_km) && R50_km > 0) {
      RMW_km <- 0.40 * R50_km
    } else if (is.finite(R34_km) && R34_km > 0) {
      RMW_km <- 0.50 * R34_km
    } else {
      RMW_nm <- 35 - 0.11 * Vmax_kt + 0.00043 * Vmax_kt^2
      RMW_km <- pmax(10, pmin(150, RMW_nm * 1.852))
    }
  }
  RMW_km <- pmax(5, pmin(200, RMW_km))
  
  # Holland B parameter
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
  
  if (r_km < 0.1) return(0)
  
  r_norm <- r_km / RMW_km
  if (r_norm < 1.0) {
    V_gradient <- sqrt((1 / r_norm)^B * exp(1 - (1 / r_norm)^B))
  } else {
    V_gradient <- sqrt((RMW_km / r_km)^B * exp(1 - (RMW_km / r_km)^B))
  }
  V_site_kt <- Vmax_kt * V_gradient
  
  # Calibration using R34 (outer field)
  if (is.finite(R34_km) && R34_km > 0 && r_km > RMW_km * 1.5) {
    V_at_R34_model <- Vmax_kt * sqrt((RMW_km / R34_km)^B * exp(1 - (RMW_km / R34_km)^B))
    if (V_at_R34_model > 25 && V_at_R34_model < 50) {
      V_site_kt <- V_site_kt * (34 / V_at_R34_model)
    }
  }
  
  # Outer cutoff
  R_outer <- if (is.finite(R34_km) && R34_km > 0) 1.8 * R34_km else 300
  if (r_km > R_outer) {
    V_site_kt <- V_site_kt * exp(-2 * (r_km - R_outer) / R_outer)
  }
  
  pmax(0, pmin(Vmax_kt, V_site_kt))
}

compute_storm_heading <- function(df) {
  df %>%
    arrange(SID, ISO_TIME) %>%
    group_by(SID) %>%
    mutate(
      LAT_next = lead(LAT),
      LON_next = lead(LON),
      storm_heading = geosphere::bearing(cbind(LON, LAT), cbind(LON_next, LAT_next)),
      storm_heading = (storm_heading + 360) %% 360,
      storm_heading = if_else(is.na(storm_heading), lag(storm_heading), storm_heading)
    ) %>%
    ungroup() %>%
    select(-LAT_next, -LON_next)
}

add_forward_motion_asymmetry <- function(
    V_site_base_kt,
    storm_speed_kt,
    r_km,
    bearing_to_target,
    storm_heading,
    RMW_km = 40
) {
  if (!is.finite(V_site_base_kt) || !is.finite(storm_speed_kt)) return(V_site_base_kt)
  if (!is.finite(bearing_to_target) || !is.finite(storm_heading)) return(V_site_base_kt)
  if (storm_speed_kt < 0.1) return(V_site_base_kt)
  
  angle_diff <- bearing_to_target - storm_heading
  angle_diff <- ((angle_diff + 180) %% 360) - 180
  theta_rad <- angle_diff * pi / 180
  
  r_norm <- r_km / RMW_km
  if (!is.finite(r_norm) || r_norm <= 0) return(V_site_base_kt)
  
  if (r_norm < 1.0) {
    K <- 0.4 + 0.25 * r_norm
  } else if (r_norm <= 3.0) {
    K <- 0.65 - 0.15 * (r_norm - 1.0) / 2.0
  } else {
    K <- 0.50 * exp(-0.2 * (r_norm - 3.0))
  }
  K <- pmax(0.1, pmin(0.7, K))
  
  V_motion_kt <- K * storm_speed_kt * cos(theta_rad)
  
  V_site_kt <- V_site_base_kt + V_motion_kt
  pmax(0, V_site_kt)
}

compute_site_winds_full <- function(df, t_lat, t_lon) {
  df <- df %>%
    mutate(
      bearing_to_target = calculate_bearing(LAT, LON, t_lat, t_lon),
      quadrant = get_quadrant(bearing_to_target),
      
      R34_nm = get_directional_radius(quadrant, USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
      R50_nm = get_directional_radius(quadrant, USA_R50_NE, USA_R50_SE, USA_R50_SW, USA_R50_NW),
      R64_nm = get_directional_radius(quadrant, USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW),
      
      R34_km = R34_nm * 1.852,
      R50_km = R50_nm * 1.852,
      R64_km = R64_nm * 1.852
    )
  
  radii_fixed <- enforce_monotone_radii(df$R34_km, df$R50_km, df$R64_km)
  df$R34_km <- radii_fixed$R34_km
  df$R50_km <- radii_fixed$R50_km
  df$R64_km <- radii_fixed$R64_km
  
  df <- compute_storm_heading(df)
  
  df <- df %>%
    mutate(
      Vmax_kt = USA_WIND,
      V_site_kt = mapply(
        estimate_site_wind_piecewise,
        Vmax = Vmax_kt,
        r_km = dist_km,
        R34_km = R34_km,
        R50_km = R50_km,
        R64_km = R64_km,
        MoreArgs = list(r0_mult = 1.5)
      )
    )
  
  df
}

# =============================================================================
# 3) IBTrACS reader
# =============================================================================

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

# =============================================================================
# 4) ENHANCED: Refined severity classification (Saffir-Simpson scale)
# =============================================================================

classify_severity_saffir_simpson <- function(V_site_max_kt) {
  # Refined 5-class Saffir-Simpson hurricane categorization
  # Based on actual wind exposure at site (not center intensity)
  #
  # Thresholds:
  #   none:  < 34 kt (below tropical storm strength)
  #   TS:    34-63 kt (tropical storm)
  #   CAT1:  64-82 kt (Category 1 hurricane)
  #   CAT2:  83-95 kt (Category 2 hurricane)  
  #   CAT3:  96-112 kt (Category 3 hurricane - major)
  #   CAT4:  113-136 kt (Category 4 hurricane - major)
  #   CAT5:  137+ kt (Category 5 hurricane - major)
  
  dplyr::case_when(
    is.na(V_site_max_kt)      ~ "unknown",
    V_site_max_kt < 34        ~ "none",
    V_site_max_kt < 64        ~ "TS",
    V_site_max_kt < 83        ~ "CAT1",
    V_site_max_kt < 96        ~ "CAT2",
    V_site_max_kt < 113       ~ "CAT3",
    V_site_max_kt < 137       ~ "CAT4",
    TRUE                      ~ "CAT5"
  )
}

classify_severity_simple <- function(V_site_max_kt) {
  # Backward-compatible 3-class system
  dplyr::case_when(
    is.na(V_site_max_kt) ~ "unknown",
    V_site_max_kt < 34   ~ "none",
    V_site_max_kt < 64   ~ "TS",
    TRUE                 ~ "HUR64plus"
  )
}

classify_severity_major <- function(V_site_max_kt) {
  # 4-class system: none, TS, HUR (Cat 1-2), MAJOR (Cat 3+)
  dplyr::case_when(
    is.na(V_site_max_kt)      ~ "unknown",
    V_site_max_kt < 34        ~ "none",
    V_site_max_kt < 64        ~ "TS",
    V_site_max_kt < 96        ~ "HUR",
    TRUE                      ~ "MAJOR"
  )
}

# =============================================================================
# 5) ENHANCED: Integrated Wind Exposure Index
# =============================================================================

compute_wind_exposure_index <- function(V_site_kt, dt_hours = 6) {
  # Integrated Wind Exposure Index (IWE)
  # 
  # Accumulates wind exposure over time using a damage-potential formulation.
  # Based on Emanuel's Power Dissipation Index (PDI) concept but adapted for
  # site-level exposure assessment.
  #
  # IWE = sum over time of: V^3 * dt  (where V >= 34 kt)
  #
  # Units: kt³·hours (can normalize by dividing by 1000 for readability)
  #
  # Physical interpretation:
  #   - V³ scaling reflects nonlinear damage relationship with wind speed
  #   - Only winds >= 34 kt contribute (destructive wind threshold)
  #   - Integration over exposure duration captures cumulative impact
  #
  # Returns: scalar IWE value (kt³·hours / 1000 for scaling)
  
  V <- V_site_kt[is.finite(V_site_kt)]
  V <- V[V >= 34]  # Only destructive winds contribute
  
  if (length(V) == 0) return(0)
  
  IWE_raw <- sum(V^3) * dt_hours
  IWE_raw / 1000  # Scale for readability
}

compute_destructive_wind_hours <- function(V_site_kt, threshold_kt = 34) {
  # Simple exposure duration at or above destructive wind threshold
  # More interpretable metric for communication
  sum(V_site_kt >= threshold_kt, na.rm = TRUE) * 6  # 6-hour timesteps
}

compute_peak_gust_estimate <- function(V_site_max_kt, gust_factor = 1.3) {
  # Estimate 3-sec peak gust from sustained wind
  # Typical gust factors: 1.2-1.4 depending on terrain/exposure
  # Default 1.3 is standard for open water/coastal exposure
  V_site_max_kt * gust_factor
}

# =============================================================================
# 6) Event aggregation with enhanced metrics
# =============================================================================

make_storm_events_enhanced <- function(df, thresholds, use_saffir_simpson = TRUE) {
  # Enhanced event aggregation with:
  #   - Refined Saffir-Simpson severity classification
  #   - Integrated Wind Exposure Index
  #   - Peak gust estimates
  #   - Multiple exposure duration metrics
  
  thr_ts  <- thresholds$thr_ts
  thr_50  <- thresholds$thr_50
  thr_hur <- thresholds$thr_hur
  
  # optional island-specific disruption thresholds
  thr_port  <- thresholds$thr_port
  thr_infra <- thresholds$thr_infra
  
  df <- df %>% distinct(SID, ISO_TIME, .keep_all = TRUE)
  
  out <- df %>%
    group_by(SID, SEASON) %>%
    summarise(
      year = as.integer(first(SEASON)),
      name = dplyr::first(na.omit(NAME)),
      min_dist_km = min_na(dist_km),
      
      # Peak wind metrics
      V_site_max_kt = max_na(V_site_kt),
      V_site_p95_kt = if (all(is.na(V_site_kt))) NA_real_ else
        as.numeric(stats::quantile(V_site_kt, probs = 0.95, na.rm = TRUE, names = FALSE)),
      V_site_mean_kt = mean_na(V_site_kt),
      
      # Peak gust estimate
      V_gust_max_kt = compute_peak_gust_estimate(max_na(V_site_kt), gust_factor = 1.3),
      
      # ENHANCED: Integrated Wind Exposure Index
      wind_exposure_index = compute_wind_exposure_index(V_site_kt, dt_hours = 6),
      
      # Exposure duration at key thresholds
      hours_ge34 = sum(V_site_kt >= thr_ts, na.rm = TRUE) * 6,
      hours_ge50 = sum(V_site_kt >= thr_50, na.rm = TRUE) * 6,
      hours_ge64 = sum(V_site_kt >= thr_hur, na.rm = TRUE) * 6,
      hours_ge96 = sum(V_site_kt >= 96, na.rm = TRUE) * 6,  # Major hurricane threshold
      
      # Threshold exceedance flags
      exposed_to_34kt = any(V_site_kt >= thr_ts, na.rm = TRUE),
      exposed_to_50kt = any(V_site_kt >= thr_50, na.rm = TRUE),
      exposed_to_64kt = any(V_site_kt >= thr_hur, na.rm = TRUE),
      exposed_to_96kt = any(V_site_kt >= 96, na.rm = TRUE),
      
      # QA metrics
      any_Vsite_known = any(!is.na(V_site_kt)),
      frac_Vsite_unknown = mean(is.na(V_site_kt)),
      n_observations = n(),
      
      # Storm metadata
      storm_speed_mean = mean_na(STORM_SPEED),
      max_center_wind  = max_na(USA_WIND),
      
      .groups = "drop"
    ) %>%
    filter(!is.na(year))
  
  # Add severity classification (user-selectable scheme)
  if (use_saffir_simpson) {
    out <- out %>%
      mutate(
        severity = classify_severity_saffir_simpson(V_site_max_kt),
        severity_simple = classify_severity_simple(V_site_max_kt),
        severity_major = classify_severity_major(V_site_max_kt)
      )
  } else {
    out <- out %>%
      mutate(
        severity = classify_severity_simple(V_site_max_kt),
        severity_saffir = classify_severity_saffir_simpson(V_site_max_kt),
        severity_major = classify_severity_major(V_site_max_kt)
      )
  }
  
  # Add optional disruption channels (port/infra), if provided
  if (is.finite(thr_port)) {
    out <- out %>%
      mutate(
        thr_port_kt = thr_port,
        port_disrupt = V_site_max_kt >= thr_port
      )
  } else {
    out <- out %>% mutate(thr_port_kt = NA_real_, port_disrupt = NA)
  }
  
  if (is.finite(thr_infra)) {
    out <- out %>%
      mutate(
        thr_infra_kt = thr_infra,
        infra_disrupt = V_site_max_kt >= thr_infra
      )
  } else {
    out <- out %>% mutate(thr_infra_kt = NA_real_, infra_disrupt = NA)
  }
  
  out
}

# =============================================================================
# 7) Rate model helpers (flexible severity classes)
# =============================================================================

compute_annual_counts <- function(events, severity_classes = NULL) {
  # Flexible annual counts for any severity classification scheme
  # If severity_classes is NULL, uses all non-"unknown"/"none" classes found
  
  if (is.null(severity_classes)) {
    severity_classes <- setdiff(unique(events$severity), c("unknown", "none"))
  }
  
  events %>%
    filter(severity %in% severity_classes) %>%
    distinct(year, severity, SID) %>%
    count(year, severity, name = "n_events") %>%
    tidyr::complete(
      year = full_seq(range(year), 1),
      severity = severity_classes,
      fill = list(n_events = 0)
    )
}

compute_lambda_table <- function(annual_counts) {
  annual_counts %>%
    group_by(severity) %>%
    summarise(
      lambda = mean(n_events),
      n_years = n(),
      p_at_least_one = 1 - exp(-lambda),
      p_zero = exp(-lambda),
      .groups = "drop"
    )
}

estimate_k_hat <- function(annual_counts) {
  annual_total <- annual_counts %>%
    group_by(year) %>%
    summarise(N = sum(n_events), .groups = "drop")
  
  mu <- mean(annual_total$N)
  va <- stats::var(annual_total$N)
  
  k_hat <- if (is.finite(va) && va > mu && mu > 0) mu^2 / (va - mu) else 1e6
  list(k_hat = k_hat, annual_total = annual_total, mu = mu, var = va)
}

simulate_multilevel_counts <- function(lambda_table, k_hat, n_years_sim = 1000) {
  # Generalized simulation for multiple severity classes
  # Uses shared annual activity factor across all classes
  
  A <- rgamma(n_years_sim, shape = k_hat, rate = k_hat)
  
  sim <- tibble(sim_year = 1:n_years_sim, A = A)
  
  for (sev in lambda_table$severity) {
    lam <- lambda_table$lambda[lambda_table$severity == sev]
    col_name <- paste0("n_", sev)
    sim[[col_name]] <- rpois(n_years_sim, lam * A)
  }
  
  sim
}

# =============================================================================
# 8) Orchestrator: run hazard model for all islands + return tidy outputs
# =============================================================================

run_hazard_model_enhanced <- function(
    cfg,
    targets,
    per_target_cfg = list(),
    severity_scheme = "saffir_simpson",  # "saffir_simpson", "simple", or "major"
    severity_classes = NULL
) {
  # ENHANCED orchestrator with:
  #   - Selectable severity classification schemes
  #   - Integrated wind exposure metrics
  #   - Flexible rate modeling for any severity classes
  
  # Determine which severity classes to model
  if (is.null(severity_classes)) {
    severity_classes <- switch(
      severity_scheme,
      "saffir_simpson" = c("TS", "CAT1", "CAT2", "CAT3", "CAT4", "CAT5"),
      "simple" = c("TS", "HUR64plus"),
      "major" = c("TS", "HUR", "MAJOR"),
      c("TS", "CAT1", "CAT2", "CAT3", "CAT4", "CAT5")  # default
    )
  }
  
  use_saffir_simpson <- (severity_scheme == "saffir_simpson")
  
  # read once
  ib_sub <- read_ibtracs_clean(cfg$ibtracs_path, cfg)
  
  # storage
  trackpoints_list <- setNames(vector("list", nrow(targets)), targets$name)
  events_list      <- setNames(vector("list", nrow(targets)), targets$name)
  annual_counts_list <- setNames(vector("list", nrow(targets)), targets$name)
  lambda_list        <- setNames(vector("list", nrow(targets)), targets$name)
  sim_list           <- setNames(vector("list", nrow(targets)), targets$name)
  kinfo_list         <- setNames(vector("list", nrow(targets)), targets$name)
  
  for (i in seq_len(nrow(targets))) {
    loc <- targets[i, ]
    island <- as.character(loc$name)
    
    # merge global thresholds with per-island overrides (if any)
    thr <- list(
      thr_ts  = cfg$thr_ts,
      thr_50  = cfg$thr_50,
      thr_hur = cfg$thr_hur,
      thr_port  = NA_real_,
      thr_infra = NA_real_
    )
    if (!is.null(per_target_cfg[[island]])) {
      thr[names(per_target_cfg[[island]])] <- per_target_cfg[[island]]
    }
    
    # gate: compute distance once
    dat_loc <- ib_sub %>%
      mutate(dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon)) %>%
      filter(dist_km <= cfg$gate_km)
    
    # compute winds + radii etc
    dat_loc <- compute_site_winds_full(dat_loc, loc$lat, loc$lon)
    
    trackpoints_list[[island]] <- dat_loc
    
    # ENHANCED: events with wind exposure index
    ev <- make_storm_events_enhanced(
      dat_loc,
      thresholds = thr,
      use_saffir_simpson = use_saffir_simpson
    ) %>%
      mutate(island = island) %>%
      relocate(island, .before = SID)
    
    # restrict years and unknowns
    ev <- ev %>%
      filter(year >= cfg$min_year) %>%
      filter(any_Vsite_known) %>%
      filter(severity != "unknown")
    
    events_list[[island]] <- ev
    
    # annual counts + rate model + multilevel sim (per island)
    ac <- compute_annual_counts(ev, severity_classes = severity_classes)
    lt <- compute_lambda_table(ac)
    kinfo <- estimate_k_hat(ac)
    sim <- simulate_multilevel_counts(lt, kinfo$k_hat, n_years_sim = cfg$n_years_sim)
    
    annual_counts_list[[island]] <- ac %>% mutate(island = island) %>% relocate(island, .before = year)
    lambda_list[[island]]        <- lt %>% mutate(island = island) %>% relocate(island, .before = severity)
    sim_list[[island]]           <- sim %>% mutate(island = island) %>% relocate(island, .before = sim_year)
    kinfo_list[[island]]         <- tibble(
      island = island,
      k_hat = kinfo$k_hat,
      annual_mean = kinfo$mu,
      annual_var = kinfo$var
    )
  }
  
  # tidy outputs across islands
  events_all <- dplyr::bind_rows(events_list)
  annual_counts_all <- dplyr::bind_rows(annual_counts_list)
  lambda_all <- dplyr::bind_rows(lambda_list)
  sim_all <- dplyr::bind_rows(sim_list)
  kinfo_all <- dplyr::bind_rows(kinfo_list)
  
  list(
    trackpoints = trackpoints_list,
    events_by_island = events_list,
    events_all = events_all,
    annual_counts_all = annual_counts_all,
    lambda_all = lambda_all,
    sim_all = sim_all,
    kinfo_all = kinfo_all,
    severity_scheme = severity_scheme,
    severity_classes = severity_classes
  )
}