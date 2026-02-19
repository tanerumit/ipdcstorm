# R/functions_hazard.R
# Tropical Cyclone Hazard Model functions (IBTrACS + directional radii + Holland + asymmetry)

#### TODO 

# Check units of STORM_SPEED is actually m/s or km/h, your asymmetry term will be wrong by a large factor.

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(geosphere)
library(readr)

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
    RMW_km,              # <-- REQUIRED
    Pn = 1013,
    Pc = NA
) {
  if (!is.finite(Vmax_kt) || Vmax_kt <= 0) return(NA_real_)
  if (!is.finite(r_km) || r_km < 0) return(NA_real_)
  if (!is.finite(RMW_km) || RMW_km <= 0) return(NA_real_)   # <-- enforce
  
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
  
  # center behavior: do NOT return 0
  if (r_km < 0.1) return(pmax(0, Vmax_kt))
  
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
      storm_heading = geosphere::bearing(
        cbind(LON, LAT),
        cbind(LON_next, LAT_next)
      ),
      storm_heading = (storm_heading + 360) %% 360,
      storm_heading = if_else(
        is.finite(LAT_next) & is.finite(LON_next),
        geosphere::bearing(cbind(LON, LAT), cbind(LON_next, LAT_next)),
        NA_real_
      )
    ) %>%
    # fill missing headings within each SID (handles last point, missing coords, etc.)
    tidyr::fill(storm_heading, .direction = "downup") %>%
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
  pmax(0, V_site_base_kt + V_motion_kt)
}

# =============================================================================
# 3) Core per-location computation (trackpoints -> V_site_kt)
# =============================================================================

compute_site_winds_full <- function(df, target_lat, target_lon) {
  
  df <- df %>%
    mutate(
      dist_km = if_else(is.na(dist_km), dist_to_target(LAT, LON, target_lat, target_lon), dist_km),
      bearing_to_target = calculate_bearing(LAT, LON, target_lat, target_lon),
      quadrant = get_quadrant(bearing_to_target),
      
      R34_nm = get_directional_radius(quadrant, USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
      R50_nm = get_directional_radius(quadrant, USA_R50_NE, USA_R50_SE, USA_R50_SW, USA_R50_NW),
      R64_nm = get_directional_radius(quadrant, USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW),
      
      R34_km = R34_nm * 1.852,
      R50_km = R50_nm * 1.852,
      R64_km = R64_nm * 1.852,
      
      Vmax_kt = USA_WIND
    )
  
  # enforce monotone radii
  mm <- enforce_monotone_radii(df$R34_km, df$R50_km, df$R64_km)
  df <- df %>% mutate(R34_km = mm$R34_km, R50_km = mm$R50_km, R64_km = mm$R64_km)
  
  # storm heading + RMW
  df <- compute_storm_heading(df)
  
  df <- df %>%
    mutate(
      RMW_km = case_when(
        is.finite(R64_km) & R64_km > 0 ~ 0.35 * R64_km,
        is.finite(R50_km) & R50_km > 0 ~ 0.40 * R50_km,
        is.finite(R34_km) & R34_km > 0 ~ 0.50 * R34_km,
        TRUE ~ pmax(10, pmin(150, (35 - 0.11 * Vmax_kt + 0.00043 * Vmax_kt^2) * 1.852))
      )
    )
  
  df$V_site_symmetric_kt <- mapply(
    estimate_site_wind_holland,
    Vmax_kt = df$Vmax_kt,
    r_km    = df$dist_km,
    R34_km  = df$R34_km,
    R50_km  = df$R50_km,
    R64_km  = df$R64_km,
    RMW_km  = df$RMW_km,
    MoreArgs = list(Pn = 1013)
  )
  
  df$V_site_kt <- mapply(
    add_forward_motion_asymmetry,
    V_site_base_kt = df$V_site_symmetric_kt,
    storm_speed_kt = df$STORM_SPEED,
    r_km = df$dist_km,
    bearing_to_target = df$bearing_to_target,
    storm_heading = df$storm_heading,
    RMW_km = df$RMW_km
  )
  
  df
}


#' Read and Clean IBTrACS CSV (North Atlantic / USA fields + pressure)
#'
#' @description
#' Reads an IBTrACS "ALL" or basin CSV export and produces a cleaned tibble with
#' consistent numeric types and derived pressure fields that are coherent with
#' the `USA_*` wind/radii fields (Pc from USA_PRES with WMO fallback; optional
#' structure fields USA_POCI/USA_ROCI/USA_RMW).
#'
#' @param ibtracs_csv Character scalar. Path to IBTrACS CSV file.
#' @param basin Character vector or NULL. If provided, keep only BASIN in this set
#'   (e.g., "NA" for North Atlantic). Use NULL to keep all basins.
#' @param season Integer vector or NULL. If provided, keep only SEASON in this set.
#' @param keep_all Logical. If TRUE, keep all original IBTrACS columns in output
#'   and append cleaned/derived columns. If FALSE, return a compact set used by
#'   hazard workflows.
#' @param verbose Logical. If TRUE, prints basic read/filter messages.
#'
#' @return A tibble with cleaned core fields and added:
#'   - `pres_hpa` (Pc; USA_PRES then WMO_PRES)
#'   - `pres_source` ("USA", "WMO", or NA)
#'   - `usa_pres_hpa`, `wmo_pres_hpa` (numeric)
#'   - `poci_hpa`, `roci_km`, `rmw_km` (numeric; from USA_* structure fields)
#'
#' @export
read_ibtracs_clean <- function(ibtracs_csv,
                               basin = "NA",
                               season = NULL,
                               keep_all = FALSE,
                               verbose = TRUE) {
  
  if (!is.character(ibtracs_csv) || length(ibtracs_csv) != 1L || !nzchar(ibtracs_csv)) {
    stop("`ibtracs_csv` must be a non-empty character scalar path.")
  }
  if (!file.exists(ibtracs_csv)) {
    stop("File not found: ", ibtracs_csv)
  }
  
  # Minimal dependency footprint: readr + dplyr + lubridate are typical in this workflow
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package `readr` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package `lubridate` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required()")

  to_num <- function(x) {
    suppressWarnings(as.numeric(x))
  }

  # Read
  if (isTRUE(verbose)) message("[IBTrACS] Reading CSV: ", ibtracs_csv)

  df <- readr::read_csv(
    file = ibtracs_csv,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 200000
  )

  # Basic column presence checks (fail early, clearly)
  required_cols <- c("SID", "SEASON", "BASIN", "ISO_TIME",
                     "USA_LAT", "USA_LON", "USA_WIND",
                     "USA_R34_NE", "USA_R34_SE", "USA_R34_SW", "USA_R34_NW")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("IBTrACS file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filters
  if (!is.null(basin)) {
    df <- dplyr::filter(df, .data$BASIN %in% basin)
  }
  if (!is.null(season)) {
    df <- dplyr::filter(df, .data$SEASON %in% season)
  }

  # Clean core fields + derive pressure/structure fields
  df2 <- df |>
    dplyr::mutate(
      # Time
      iso_time = lubridate::ymd_hms(.data$ISO_TIME, tz = "UTC", quiet = TRUE),

      # Core position/intensity (USA track)
      lat = to_num(.data$USA_LAT),
      lon = to_num(.data$USA_LON),
      wind_kt = to_num(.data$USA_WIND),

      # Radii (nm in IBTrACS; keep as-is, or convert later explicitly)
      r34_ne_nm = to_num(.data$USA_R34_NE),
      r34_se_nm = to_num(.data$USA_R34_SE),
      r34_sw_nm = to_num(.data$USA_R34_SW),
      r34_nw_nm = to_num(.data$USA_R34_NW),

      r50_ne_nm = dplyr::if_else("USA_R50_NE" %in% names(df), to_num(.data$USA_R50_NE), NA_real_),
      r50_se_nm = dplyr::if_else("USA_R50_SE" %in% names(df), to_num(.data$USA_R50_SE), NA_real_),
      r50_sw_nm = dplyr::if_else("USA_R50_SW" %in% names(df), to_num(.data$USA_R50_SW), NA_real_),
      r50_nw_nm = dplyr::if_else("USA_R50_NW" %in% names(df), to_num(.data$USA_R50_NW), NA_real_),

      r64_ne_nm = dplyr::if_else("USA_R64_NE" %in% names(df), to_num(.data$USA_R64_NE), NA_real_),
      r64_se_nm = dplyr::if_else("USA_R64_SE" %in% names(df), to_num(.data$USA_R64_SE), NA_real_),
      r64_sw_nm = dplyr::if_else("USA_R64_SW" %in% names(df), to_num(.data$USA_R64_SW), NA_real_),
      r64_nw_nm = dplyr::if_else("USA_R64_NW" %in% names(df), to_num(.data$USA_R64_NW), NA_real_),

      # ---- Pressure fields (hPa) ----
      usa_pres_hpa = dplyr::if_else("USA_PRES" %in% names(df), to_num(.data$USA_PRES), NA_real_),
      wmo_pres_hpa = dplyr::if_else("WMO_PRES" %in% names(df), to_num(.data$WMO_PRES), NA_real_),

      pres_hpa = dplyr::case_when(
        is.finite(.data$usa_pres_hpa) ~ .data$usa_pres_hpa,
        is.finite(.data$wmo_pres_hpa) ~ .data$wmo_pres_hpa,
        TRUE ~ NA_real_
      ),
      pres_source = dplyr::case_when(
        is.finite(.data$usa_pres_hpa) ~ "USA",
        is.finite(.data$wmo_pres_hpa) ~ "WMO",
        TRUE ~ NA_character_
      ),

      # ---- Structure fields (optional, USA only) ----
      poci_hpa = dplyr::if_else("USA_POCI" %in% names(df), to_num(.data$USA_POCI), NA_real_),
      roci_km  = dplyr::if_else("USA_ROCI" %in% names(df), to_num(.data$USA_ROCI), NA_real_),
      rmw_km   = dplyr::if_else("USA_RMW"  %in% names(df), to_num(.data$USA_RMW),  NA_real_),

      # Useful passthroughs
      storm_status = dplyr::if_else("USA_STATUS" %in% names(df), as.character(.data$USA_STATUS), NA_character_),
      storm_name   = dplyr::if_else("NAME" %in% names(df), as.character(.data$NAME), NA_character_)
    ) |>
    dplyr::arrange(.data$SID, .data$iso_time)

  # Compact vs keep-all output
  if (isTRUE(keep_all)) {
    out <- df2
  } else {
    out <- df2 |>
      dplyr::select(
        .data$SID, .data$SEASON, .data$BASIN,
        iso_time,
        storm_name, storm_status,
        lat, lon, wind_kt,
        pres_hpa, pres_source, usa_pres_hpa, wmo_pres_hpa,
        poci_hpa, roci_km, rmw_km,
        r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm,
        r50_ne_nm, r50_se_nm, r50_sw_nm, r50_nw_nm,
        r64_ne_nm, r64_se_nm, r64_sw_nm, r64_nw_nm
      )
  }

  # Light QA message
  if (isTRUE(verbose)) {
    n <- nrow(out)
    n_pres <- sum(is.finite(out$pres_hpa))
    message("[IBTrACS] Rows: ", n, " | Pc available (pres_hpa): ", n_pres)
  }

  tibble::as_tibble(out)
}






#################################################################################









read_ibtracs_clean_old <- function(path, cfg) {
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
# 4) Event aggregation + severity + extra threshold channels
# =============================================================================

classify_severity <- function(V_site_max_kt, thr_ts = 34, thr_hur = 64) {
  dplyr::case_when(
    is.na(V_site_max_kt) ~ "unknown",
    V_site_max_kt < thr_ts ~ "none",
    V_site_max_kt < thr_hur ~ "TS",
    TRUE ~ "HUR64plus"
  )
}

make_storm_events <- function(df, thresholds) {
  # thresholds: list with at least thr_ts, thr_50, thr_hur
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
  
  # add optional disruption channels (port/infra), if provided
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
  
  out %>%
    mutate(severity = classify_severity(V_site_max_kt, thr_ts = thr_ts, thr_hur = thr_hur))
}

# =============================================================================
# 5) Rate model helpers
# =============================================================================

compute_annual_counts <- function(events, severities = c("TS", "HUR64plus")) {
  events %>%
    filter(severity %in% severities) %>%
    distinct(year, severity, SID) %>%
    count(year, severity, name = "n_events") %>%
    tidyr::complete(
      year = full_seq(range(year), 1),
      severity = severities,
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

simulate_twolevel_counts <- function(lambda_table, k_hat, n_years_sim = 1000) {
  lam_TS  <- lambda_table$lambda[lambda_table$severity == "TS"]
  lam_HUR <- lambda_table$lambda[lambda_table$severity == "HUR64plus"]
  
  A <- rgamma(n_years_sim, shape = k_hat, rate = k_hat)
  
  tibble(
    sim_year = 1:n_years_sim,
    A = A,
    n_TS = rpois(n_years_sim, lam_TS * A),
    n_HUR64plus = rpois(n_years_sim, lam_HUR * A)
  )
}

# =============================================================================
# 6) Orchestrator: run hazard model for all islands + return tidy outputs
# =============================================================================

run_hazard_model <- function(cfg, targets, per_target_cfg = list(), severities = c("TS", "HUR64plus")) {
  
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
    
    # events
    ev <- make_storm_events(dat_loc, thresholds = thr) %>%
      mutate(island = island) %>%
      relocate(island, .before = SID)
    
    # restrict years and unknowns
    ev <- ev %>%
      filter(year >= cfg$min_year) %>%
      filter(any_Vsite_known) %>%
      filter(severity != "unknown")
    
    events_list[[island]] <- ev
    
    # annual counts + rate model + 2-level sim (per island)
    ac <- compute_annual_counts(ev, severities = severities)
    lt <- compute_lambda_table(ac)
    kinfo <- estimate_k_hat(ac)
    sim <- simulate_twolevel_counts(lt, kinfo$k_hat, n_years_sim = cfg$n_years_sim)
    
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
    kinfo_all = kinfo_all
  )
}
