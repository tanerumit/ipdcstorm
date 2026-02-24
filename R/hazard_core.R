# =============================================================================
# Script overview: core hazard computations
# - dist_to_target(): great-circle distance from track points to a fixed site.
# - calculate_bearing(): bearing from track points to the site.
# - .get_quadrant(): map bearing to NE/SE/SW/NW quadrants.
# - .get_directional_radius(): select quadrant-specific wind radii.
# - .enforce_monotone_radii(): ensure R64 <= R50 <= R34.
# - .estimate_site_wind_holland(): Holland-type radial wind profile.
# - compute_storm_heading(): storm motion heading per track point.
# - .add_forward_motion_asymmetry(): apply forward-motion wind asymmetry.
# - compute_site_winds_full(): end-to-end site wind estimation for track points.
# - classify_severity(): storm class from peak site wind.
# - make_storm_events(): aggregate track points into storm-level events.
# - compute_annual_counts(): annual event counts by storm class.
# - compute_lambda_table(): Poisson rate summaries by storm class.
# - estimate_k_hat(): overdispersion estimate for annual totals.
# =============================================================================

# =============================================================================
# 1) Geometry helpers (distance, bearing, quadrants, directional radii)
# =============================================================================

#' Distance from storm center to target (km), robust
#'
#' @description
#' Computes Haversine distance (km) from each storm track point to a fixed
#' target coordinate. Missing coordinates produce NA, and empty inputs return
#' numeric(0).
#'
#' @param lat,lon Numeric vectors of storm lat/lon in degrees.
#' @param t_lat,t_lon Numeric scalars; target lat/lon in degrees.
#'
#' @return Numeric vector of distances in km (same length as lat/lon).
#'
#' @importFrom geosphere distHaversine
#' @export
dist_to_target <- function(lat, lon, t_lat, t_lon) {
  lat <- as.numeric(lat)
  lon <- as.numeric(lon)

  n <- length(lat)
  if (n == 0L) return(numeric(0))

  ok <- is.finite(lat) & is.finite(lon)
  out <- rep(NA_real_, n)

  if (any(ok)) {
    p1 <- cbind(lon[ok], lat[ok])
    p2 <- cbind(rep(t_lon, sum(ok)), rep(t_lat, sum(ok)))
    out[ok] <- geosphere::distHaversine(p1, p2) / 1000
  }

  out
}

#' Bearing from storm center to target (degrees), robust
#'
#' @description
#' Computes great-circle initial bearing (degrees) from storm points to a fixed
#' target coordinate. Missing coordinates produce NA, and empty inputs return
#' numeric(0).
#'
#' @param lat,lon Numeric vectors of storm lat/lon in degrees.
#' @param t_lat,t_lon Numeric scalars; target lat/lon in degrees.
#'
#' @return Numeric vector of bearings in degrees (same length as lat/lon).
#'
#' @importFrom geosphere bearing
#' @export

calculate_bearing <- function(lat, lon, t_lat, t_lon) {
  lat <- as.numeric(lat)
  lon <- as.numeric(lon)

  n <- length(lat)
  if (n == 0L) return(numeric(0))

  ok <- is.finite(lat) & is.finite(lon)
  out <- rep(NA_real_, n)

  if (any(ok)) {
    p1 <- cbind(lon[ok], lat[ok])
    p2 <- cbind(rep(t_lon, sum(ok)), rep(t_lat, sum(ok)))
    out[ok] <- geosphere::bearing(p1, p2)
  }

  out
}

#' Convert bearing to meteorological quadrant
#'
#' @param bearing Numeric vector of bearings in degrees.
#' @return Character vector in \code{c("NE", "SE", "SW", "NW")} (or \code{NA}).
#' @keywords internal
.get_quadrant <- function(bearing) {
  b <- (bearing + 360) %% 360
  dplyr::case_when(
    b >= 0   & b < 90  ~ "NE",
    b >= 90  & b < 180 ~ "SE",
    b >= 180 & b < 270 ~ "SW",
    b >= 270 & b < 360 ~ "NW",
    TRUE ~ NA_character_
  )
}

#' Select directional wind radius based on quadrant
#'
#' @param quadrant Character vector: NE/SE/SW/NW.
#' @param r_ne,r_se,r_sw,r_nw Numeric vectors/scalars of radii for each quadrant (nm or km).
#'
#' @return Numeric vector of selected radii.
#' @keywords internal

.get_directional_radius <- function(quadrant, r_ne, r_se, r_sw, r_nw) {
  r_ne <- suppressWarnings(as.numeric(r_ne))
  r_se <- suppressWarnings(as.numeric(r_se))
  r_sw <- suppressWarnings(as.numeric(r_sw))
  r_nw <- suppressWarnings(as.numeric(r_nw))

  dplyr::case_when(
    quadrant == "NE" ~ r_ne,
    quadrant == "SE" ~ r_se,
    quadrant == "SW" ~ r_sw,
    quadrant == "NW" ~ r_nw,
    TRUE ~ NA_real_
  )
}

#' Enforce monotonicity of wind radii (R64 <= R50 <= R34)
#'
#' @param R34_km,R50_km,R64_km Numeric vectors of radii in km.
#' @return List with corrected R34_km, R50_km, R64_km.
#' @keywords internal
.enforce_monotone_radii <- function(R34_km, R50_km, R64_km) {
  R64_km <- dplyr::if_else(is.finite(R64_km) & is.finite(R50_km), pmin(R64_km, R50_km), R64_km)
  R50_km <- dplyr::if_else(is.finite(R50_km) & is.finite(R34_km), pmin(R50_km, R34_km), R50_km)
  list(R34_km = R34_km, R50_km = R50_km, R64_km = R64_km)
}

# =============================================================================
# 2) Wind field models (symmetric + forward-motion asymmetry)
# =============================================================================

# =============================================================================
# 2a) Climatological R34 infill
# =============================================================================

#' Estimate climatological R34 (km) from maximum wind using Knaff et al. (2015)
#'
#' @description
#' Provides an empirical estimate of the 34-kt wind radius for storms without
#' observed radii data (primarily pre-2004). Based on the Knaff, Sampson & Chirokova
#' (2015) statistical relationships for the North Atlantic.
#'
#' The relationship captures the well-known expansion of the wind field with
#' intensity up to ~80 kt, then contraction for the most intense compact storms.
#'
#' @param Vmax_kt Numeric vector; maximum sustained wind (kt).
#' @param lat Numeric vector; storm latitude (degrees N). Used for size correction.
#'
#' @return Numeric vector of estimated R34 in km. NA where Vmax < 34 kt.
#' @export
estimate_R34_climo <- function(Vmax_kt, lat = 18) {
  # Knaff et al. (2015) approximate fit for Atlantic basin:
  #   R34_nm ÃƒÂ¢Ã¢â‚¬Â°Ã‹â€  a + bÃƒâ€šÃ‚Â·(Vmax-34) + cÃƒâ€šÃ‚Â·(Vmax-34)^2 + dÃƒâ€šÃ‚Â·|lat-25|
  # Coefficients estimated from Fig 4/Table 3 of Knaff et al. 2015
  a <- 68    # base R34 at TS threshold (nm) ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Atlantic mean
  b <- 1.6   # linear expansion with intensity
  c <- -0.012 # quadratic contraction at high intensity (compact storms)
  d <- -0.8  # latitude correction (smaller at low latitudes)

  dV <- pmax(0, Vmax_kt - 34)
  R34_nm <- a + b * dV + c * dV^2 + d * abs(lat - 25)

  # Floor and cap

  R34_nm <- pmax(30, pmin(300, R34_nm))

  # Below TS threshold: no 34-kt radius
  R34_nm[!is.finite(Vmax_kt) | Vmax_kt < 34] <- NA_real_

  R34_nm * 1.852  # convert to km
}



# =============================================================================
# 2b) Knaff & Zehr RMW estimation (with latitude)
# =============================================================================

#' Estimate radius of maximum wind using Knaff & Zehr (2007) climatology
#'
#' @description
#' Latitude-dependent RMW estimation based on Knaff & Zehr (2007), which captures
#' the key relationship: intense low-latitude storms (like Irma at 17Ãƒâ€šÃ‚Â°N) have
#' smaller RMW than the same intensity at higher latitudes.
#'
#' @param Vmax_kt Numeric vector; maximum sustained wind (kt).
#' @param lat Numeric vector; storm latitude (degrees N).
#'
#' @return Numeric vector of estimated RMW in km.
#' @export
estimate_RMW_knaff <- function(Vmax_kt, lat = 18) {
  # Knaff & Zehr (2007) Eq. 1 (simplified, Atlantic):
  #   RMW_nm = 66.785 - 0.09102Ãƒâ€šÃ‚Â·Vmax + 1.0619Ãƒâ€šÃ‚Â·(lat - 25)
  # Valid for Vmax in [30, 185] kt

  Vmax_clamped <- pmax(30, pmin(185, Vmax_kt))
  RMW_nm <- 66.785 - 0.09102 * Vmax_clamped + 1.0619 * (lat - 25)

  # Floor/cap in nm: [5, 100]
  RMW_nm <- pmax(5, pmin(100, RMW_nm))

  RMW_nm[!is.finite(Vmax_kt)] <- NA_real_

  RMW_nm * 1.852  # convert to km
}







#' Estimate site wind using a Holland-type radial wind profile (patched)
#'
#' @description
#' Computes a gradient-wind-like radial profile. This version fixes:
#' - Holland B: uses Vickery & Wadhera (2008)-style parameterization that
#'   increases B for compact intense storms (was inverted above 120 kt).
#' - R34 calibration: widens the acceptable V_at_R34 window from \[25,50\] to
#'   \[5,60\] so weak tropical storms get calibrated.
#' - Adds climatological R34 infill when observed R34 is missing.
#'
#' @param Vmax_kt Numeric; storm maximum wind (kt).
#' @param r_km Numeric; distance from storm center to site (km).
#' @param R34_km Numeric; 34-kt wind radius (km). If NA and Vmax >= 34 kt,
#'   a climatological estimate is used.
#' @param R50_km,R64_km Numeric; 50/64-kt wind radii (km), optional.
#' @param RMW_km Numeric; radius of maximum wind (km). Required.
#' @param Pn Numeric; ambient pressure (hPa), default 1013.
#' @param Pc Numeric; central pressure (hPa), optional.
#' @param lat Numeric; storm latitude (degrees N), used for climatological infill.
#'
#' @return Numeric scalar; estimated sustained wind at site (kt).
#' @keywords internal
.estimate_site_wind_holland <- function(
    Vmax_kt,
    r_km,
    R34_km,
    R50_km = NA,
    R64_km = NA,
    RMW_km,
    Pn = 1013,
    Pc = NA,
    lat = 18
) {
  n <- length(Vmax_kt)
  # recycle scalars safely (base R recycling rules assumed by caller)
  out <- rep(NA_real_, n)

  ok <- is.finite(Vmax_kt) & Vmax_kt > 0 &
    is.finite(r_km) & r_km >= 0 &
    is.finite(RMW_km) & RMW_km > 0

  if (!any(ok)) return(out)

  Vmax_kt0 <- Vmax_kt
  r_km0 <- r_km
  RMW_km0 <- pmax(5, pmin(200, RMW_km))

  # --- FIX 1: Holland B parameterization ---
  B <- rep(NA_real_, n)

  okP <- ok & is.finite(Pc) & is.finite(Pn) & (Pc < Pn)
  if (any(okP)) {
    deltaP <- Pn[okP] - Pc[okP]
    RMW_nm <- RMW_km0[okP] / 1.852
    Bp <- 1.881 - 0.00557 * RMW_nm - 0.01097 * lat[okP] + 0.0016 * deltaP
    B[okP] <- pmax(1.0, pmin(2.5, Bp))
  }

  okH <- ok & !okP
  if (any(okH)) {
    V <- Vmax_kt0[okH]
    Bh <- rep(NA_real_, length(V))
    i1 <- V >= 100
    i2 <- (V >= 64) & !i1
    i3 <- !i1 & !i2

    if (any(i1)) Bh[i1] <- 1.4 + (V[i1] - 100) * 0.008
    if (any(i2)) Bh[i2] <- 1.3 + (V[i2] - 64)  * 0.003
    if (any(i3)) Bh[i3] <- 1.1 + (V[i3] - 20)  * 0.005

    B[okH] <- pmax(1.0, pmin(2.5, Bh))
  }

  # --- FIX 3: Climatological R34 infill ---
  R34_eff <- R34_km
  R34_is_climo <- rep(FALSE, n)

  ##################################
  lat0 <- lat
  lat0[!is.finite(lat0)] <- 18

  need_R34 <- ok & (!is.finite(R34_eff) | R34_eff <= 0) & (Vmax_kt0 >= 34)
  if (any(need_R34)) {
    R34_eff[need_R34] <- estimate_R34_climo(Vmax_kt0[need_R34], lat = lat0[need_R34])
    R34_is_climo[need_R34] <- TRUE
  }

  # center handling (keep same semantics)
  at_center <- ok & (r_km0 < 0.1)
  out[at_center] <- pmax(0, Vmax_kt0[at_center])

  # Remaining points
  use <- ok & !at_center & is.finite(B)
  if (!any(use)) return(out)

  r_norm <- r_km0[use] / RMW_km0[use]

  # Holland gradient wind (vector)
  inv_r <- 1 / r_norm
  term1 <- (inv_r) ^ B[use]
  term2 <- (RMW_km0[use] / r_km0[use]) ^ B[use]

  V_gradient <- rep(NA_real_, length(r_norm))
  inside <- r_norm < 1.0
  if (any(inside)) V_gradient[inside] <- sqrt(term1[inside] * exp(1 - term1[inside]))
  if (any(!inside)) V_gradient[!inside] <- sqrt(term2[!inside] * exp(1 - term2[!inside]))

  V_site_kt <- Vmax_kt0[use] * V_gradient

  # --- Gradient-to-surface wind correction (Powell et al. 2003, Kepert 2001) ---
  # The Holland profile follows gradient wind decay, which is slower than surface
  # wind decay. The surface/gradient ratio decreases with radius due to boundary
  # layer friction. Since we anchor at Vmax (already a surface wind), we apply
  # a relative correction normalized to 1.0 at the RMW.
  # alpha=0.20, beta=0.4 gives ~10% at 2Ãƒâ€”RMW, ~14% at 3Ãƒâ€”RMW, ~17% at 4Ãƒâ€”RMW,
  # consistent with Powell et al. (2003) surface/gradient ratios of 0.75-0.85.
  # Previous values (alpha=0.12, beta=0.5) were too gentle, contributing ~3-4 kt
  # systematic overprediction across all islands (see bias decomposition diagnostics).
  srf_alpha <- 0.15
  srf_beta  <- 0.5
  srf <- 1 - srf_alpha * (1 - exp(-srf_beta * pmax(0, r_norm - 1)))
  V_site_kt <- V_site_kt * srf

  # --- Vmax-dependent profile steepening ---
  # The Holland profile systematically overpredicts wind at 1-4x RMW for intense
  # storms. Validation diagnostic: Irma at Saba (r~1.5x RMW) gives 113 kt vs 80 kt
  # observed — a 41% overprediction — even after SRF. The Holland gradient wind
  # decay is too slow because the real boundary layer profile is much steeper
  # than the gradient wind (Kepert & Wang 2001, Powell et al. 2003).
  #
  # Design: V *= 1 - gamma * f_int(Vmax) * f_dist(r/RMW)
  #   f_int: 0 at 40 kt (no effect on weak TS), ramps to 1.0 at Cat 5
  #   f_dist: 0 at RMW (preserves direct-hit winds), ramps fast, caps at 2.0
  #   gamma: 0.35 — gives ~45% reduction for Cat5 at 1.4x+ RMW
  #
  # Worked example — Irma (155 kt) at r/RMW=1.2 (just outside eyewall):
  #   f_int = 1.0, f_dist = min(2.0, (1.2-1)*5) = 1.0
  #   factor = max(0.55, 1 - 0.35*1.0*1.0) = 0.65 → 35% reduction
  #   Holland gives ~150 kt → 150*0.65 = 98 kt (addresses near-eye overprediction)
  #
  # Effect on TS (45 kt) at r/RMW=1.5:
  #   f_int = (45-40)/115 = 0.04, f_dist = 2.0
  #   factor = max(0.55, 1 - 0.35*0.04*2.0) = 0.97 → 3% reduction (preserves rates)
  steep_gamma <- 0.35
  f_int  <- pmax(0, (Vmax_kt0[use] - 40) / 115)    # 0 at 40kt, 0.35 at Cat1, 1.0 at Cat5
  f_dist <- pmin(2.0, pmax(0, (r_norm - 1) * 5))    # 0 at RMW, 1.0 at 1.2x, 2.0 at 1.4x+
  steep_factor <- pmax(0.55, 1 - steep_gamma * f_int * f_dist)
  V_site_kt <- V_site_kt * steep_factor

  # --- FIX 2: R34 calibration (patched for overprediction) ---
  # The R34 calibration adjusts the Holland profile to match 34 kt at the observed
  # R34 distance. However, the calibration factor can be very large (>2x) for
  # compact storms where Holland decays too fast. A linear taper over-applies this
  # correction at intermediate distances (1.2-3Ã— RMW), inflating site winds by
  # 15-40% in the range that matters most for near-miss events.
  #
  # Fixes applied:
  #   (a) Quadratic taper: cal effect grows as (distance/R34)^2 instead of linearly
  #   (b) Hard cap on cal_factor at 1.4 (max 40% inflation)
  #   (c) Intensity-dependent damping: weaker calibration for strong hurricanes
  #       where the Holland inner-core structure is already well-constrained
  can_cal <- is.finite(R34_eff[use]) & (R34_eff[use] > 0) & (r_km0[use] > RMW_km0[use] * 1.2)
  if (any(can_cal)) {
    R34u <- R34_eff[use][can_cal]
    Bu <- B[use][can_cal]
    Vmaxu <- Vmax_kt0[use][can_cal]
    RMWu <- RMW_km0[use][can_cal]

    V_at_R34_model <- Vmaxu * sqrt((RMWu / R34u) ^ Bu * exp(1 - (RMWu / R34u) ^ Bu))
    # Apply same SRF at R34 distance (consistent with profile correction above)
    r_norm_R34 <- R34u / RMWu
    srf_R34 <- 1 - srf_alpha * (1 - exp(-srf_beta * pmax(0, r_norm_R34 - 1)))
    V_at_R34_model <- V_at_R34_model * srf_R34
    good <- is.finite(V_at_R34_model) & (V_at_R34_model > 5) & (V_at_R34_model < 60)

    if (any(good)) {
      cal_factor <- rep(1, length(V_at_R34_model))
      cal_factor[good] <- 34 / V_at_R34_model[good]

      # (b) Cap calibration factor to prevent extreme inflation
      cal_factor <- pmin(cal_factor, 1.4)

      # (a) Quadratic taper: effect concentrates near R34, minimal at intermediate r
      r_site <- r_km0[use][can_cal]
      R34u_safe <- pmax(R34u, RMWu * 1.5)
      taper_linear <- pmin(1.0, pmax(0.0, (r_site - RMWu * 1.2) / (R34u_safe - RMWu * 1.2)))
      taper <- taper_linear^2  # quadratic: 0.55 linear â†’ 0.30 quadratic
      cal_factor[good] <- 1 + taper[good] * (cal_factor[good] - 1)

      # (c) Intensity-dependent damping: reduce calibration for strong hurricanes
      # For Vmax > 96 kt (Cat 3+), the inner core dominates and R34 calibration
      # should not inflate intermediate-distance winds
      intensity_damp <- pmin(1.0, pmax(0.3, 1.0 - (Vmaxu - 64) / 120))
      cal_factor[good] <- 1 + intensity_damp[good] * (cal_factor[good] - 1)

      # Blend if climatological R34: reduce calibration weight since
      # climo R34 tends to overestimate (represents mean, not storm-specific
      # structure). Reduced from 0.7 to 0.5 based on validation showing
      # systematic overprediction for pre-2004 storms lacking radii data.
      is_climo_u <- R34_is_climo[use][can_cal]
      blend <- is_climo_u & good
      if (any(blend)) {
        cal_factor[blend] <- 1 + 0.5 * (cal_factor[blend] - 1)
      }

      idx <- which(can_cal)
      V_site_kt[idx] <- V_site_kt[idx] * cal_factor
    }
  }

  # Outer cutoff â€” tightened to reduce wind exposure footprint
  # Previous: 1.8Ã— observed, 2.4Ã— climo; typical R34~200km â†’ 360-480 km exposure
  # New: 1.5Ã— observed, 1.8Ã— climo; â†’ 300-360 km, closer to NHC wind field extent
  R_outer <- rep(300, length(r_norm))
  has_R34 <- is.finite(R34_eff[use]) & (R34_eff[use] > 0)
  R_outer[has_R34] <- 1.8 * R34_eff[use][has_R34]

  # climo extension (stable indexing)
  if (any(has_R34)) {
    idx_has <- which(has_R34)
    is_cl  <- R34_is_climo[use][has_R34]
    if (any(is_cl)) {
      R_outer[idx_has[is_cl]] <- 2.4 * R34_eff[use][has_R34][is_cl]
    }
  }

  beyond <- r_km0[use] > R_outer
  if (any(beyond)) {
    # Exponential decay beyond R_outer
    V_site_kt[beyond] <- V_site_kt[beyond] * exp(-2 * (r_km0[use][beyond] - R_outer[beyond]) / R_outer[beyond])
  }

  out[use] <- pmax(0, pmin(Vmax_kt0[use], V_site_kt))
  out
}


#' Compute storm heading (track motion bearing) per track point
#'
#' @description
#' Adds `heading_deg` in [0, 360) where 0=N, 90=E, 180=S, 270=W.
#' Heading is computed from point i to i+1 (forward bearing). The last point
#' uses the previous segment if available.
#'
#' @param df Data frame with columns: SID, iso_time (POSIXct), lat, lon.
#' @return Same data frame with `heading_deg` (numeric).
#' @export
compute_storm_heading <- function(df) {
  stopifnot(all(c("SID", "iso_time", "lat", "lon") %in% names(df)))

  deg2rad <- function(x) x * pi / 180
  rad2deg <- function(x) x * 180 / pi

  bearing_gc <- function(lat1, lon1, lat2, lon2) {
    phi1 <- deg2rad(lat1); phi2 <- deg2rad(lat2)
    lam1 <- deg2rad(lon1); lam2 <- deg2rad(lon2)
    dlam <- lam2 - lam1
    dlam <- (dlam + pi) %% (2 * pi) - pi

    y <- sin(dlam) * cos(phi2)
    x <- cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlam)

    brng <- atan2(y, x)
    (rad2deg(brng) + 360) %% 360
  }

  df |>
    dplyr::arrange(.data$SID, .data$iso_time) |>
    dplyr::group_by(.data$SID) |>
    dplyr::mutate(
      lat_next = dplyr::lead(.data$lat),
      lon_next = dplyr::lead(.data$lon),
      lat_prev = dplyr::lag(.data$lat),
      lon_prev = dplyr::lag(.data$lon),

      heading_fwd = dplyr::if_else(
        is.finite(.data$lat) & is.finite(.data$lon) &
          is.finite(.data$lat_next) & is.finite(.data$lon_next),
        bearing_gc(.data$lat, .data$lon, .data$lat_next, .data$lon_next),
        NA_real_
      ),

      heading_bwd = dplyr::if_else(
        is.finite(.data$lat_prev) & is.finite(.data$lon_prev) &
          is.finite(.data$lat) & is.finite(.data$lon),
        bearing_gc(.data$lat_prev, .data$lon_prev, .data$lat, .data$lon),
        NA_real_
      ),

      heading_deg = dplyr::case_when(
        is.finite(.data$heading_fwd) ~ .data$heading_fwd,
        is.finite(.data$heading_bwd) ~ .data$heading_bwd,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::any_of(c("lat_next","lon_next","lat_prev","lon_prev","heading_fwd","heading_bwd")))
}

#' Add forward-motion asymmetry to a symmetric site-wind estimate
#'
#' @description
#' Applies a direction-dependent additive term proportional to storm translation
#' speed and cos(angle) between storm heading and bearing-to-target.
#'
#' @param V_site_base_kt Numeric; symmetric site wind (kt).
#' @param storm_speed_kt Numeric; storm translation speed (kt).
#' @param r_km Numeric; distance from storm center to site (km).
#' @param bearing_to_target Numeric; bearing from storm center to site (deg).
#' @param storm_heading Numeric; storm motion heading (deg).
#' @param RMW_km Numeric; radius of maximum wind (km), used to taper asymmetry with radius.
#'
#' @return Numeric scalar; adjusted site wind (kt).
#'
#' @keywords internal
.add_forward_motion_asymmetry <- function(
    V_site_base_kt,
    storm_speed_kt,
    r_km,
    bearing_to_target,
    storm_heading,
    RMW_km = 40
) {
  n <- length(V_site_base_kt)
  out <- V_site_base_kt

  ok <- is.finite(V_site_base_kt) & is.finite(storm_speed_kt) &
    is.finite(bearing_to_target) & is.finite(storm_heading) &
    is.finite(r_km) & is.finite(RMW_km) & (RMW_km > 0)

  if (!any(ok)) return(out)

  slow <- ok & (storm_speed_kt < 0.1)
  if (any(slow)) out[slow] <- V_site_base_kt[slow]

  use <- ok & !slow
  if (!any(use)) return(out)

  angle_diff <- bearing_to_target[use] - storm_heading[use]
  angle_diff <- ((angle_diff + 180) %% 360) - 180
  theta_rad <- angle_diff * pi / 180

  r_norm <- r_km[use] / RMW_km[use]
  bad_r <- !is.finite(r_norm) | (r_norm <= 0)
  if (any(bad_r)) {
    # unchanged semantics: return base wind for bad r_norm
    idx_use <- which(use)
    out[idx_use[bad_r]] <- V_site_base_kt[idx_use[bad_r]]
  }

  good <- is.finite(r_norm) & (r_norm > 0)
  if (!any(good)) return(out)

  rn <- r_norm[good]
  K <- rep(NA_real_, length(rn))

  i1 <- rn < 1.0
  i2 <- (rn >= 1.0) & (rn <= 3.0)
  i3 <- rn > 3.0

  # K profile: peaks ~0.40 at RMW (Lin & Chavas 2012 surface-level estimates)
  # Previous values (0.50 peak) overestimated due to not accounting for the
  # gradient-to-surface correction that the asymmetry component should also
  # undergo. Reduced from 0.65Ã¢â€ â€™0.50Ã¢â€ â€™0.40 based on validation diagnostics
  # showing persistent positive intensity bias in right-front quadrant.
  if (any(i1)) K[i1] <- 0.25 + 0.15 * rn[i1]           # 0.25 to 0.40 at RMW
  if (any(i2)) K[i2] <- 0.40 - 0.08 * (rn[i2] - 1.0) / 2.0  # 0.40 to 0.32 at 3Ãƒâ€”RMW
  if (any(i3)) K[i3] <- 0.32 * exp(-0.2 * (rn[i3] - 3.0))   # decays from 0.32

  K <- pmax(0.05, pmin(0.45, K))

  idx_use <- which(use)
  idx_good <- idx_use[good]

  V_motion_kt <- K * storm_speed_kt[idx_good] * cos(theta_rad[good])
  out[idx_good] <- pmax(0, V_site_base_kt[idx_good] + V_motion_kt)
  out
}



# =============================================================================
# 3) Core per-location wind computation (trackpoints -> V_site_kt)
# =============================================================================

# =============================================================================
# Patched compute_site_winds_full (uses v2 Holland + improved RMW)
# =============================================================================

#' Compute site-level winds for track points (patched version)
#'
#' @description
#' Drop-in replacement for \code{compute_site_winds_full()} that uses the patched
#' Holland profile, improved RMW estimation, and climatological R34 infill.
#'
#' @param df Track-point data frame from \code{read_ibtracs_clean()}.
#' @param target_lat,target_lon Numeric scalars for the site location (degrees).
#'
#' @return The input data frame with added site-wind columns.
#' @export
compute_site_winds_full <- function(df, target_lat, target_lon) {

  df <- df |>
    dplyr::mutate(
      dist_km = dplyr::if_else(
        is.na(.data$dist_km),
        dist_to_target(.data$lat, .data$lon, target_lat, target_lon),
        .data$dist_km
      ),
      bearing_to_target = calculate_bearing(.data$lat, .data$lon, target_lat, target_lon),
      quadrant = .get_quadrant(.data$bearing_to_target),

      R34_nm = .get_directional_radius(.data$quadrant, .data$r34_ne_nm, .data$r34_se_nm, .data$r34_sw_nm, .data$r34_nw_nm),
      R50_nm = .get_directional_radius(.data$quadrant, .data$r50_ne_nm, .data$r50_se_nm, .data$r50_sw_nm, .data$r50_nw_nm),
      R64_nm = .get_directional_radius(.data$quadrant, .data$r64_ne_nm, .data$r64_se_nm, .data$r64_sw_nm, .data$r64_nw_nm),

      R34_km = .data$R34_nm * 1.852,
      R50_km = .data$R50_nm * 1.852,
      R64_km = .data$R64_nm * 1.852,

      Vmax_kt = .data$wind_kt
    )


  mm <- .enforce_monotone_radii(df$R34_km, df$R50_km, df$R64_km)
  df <- df |>
    dplyr::mutate(R34_km = mm$R34_km, R50_km = mm$R50_km, R64_km = mm$R64_km)

  df <- compute_storm_heading(df)

  if (!("storm_speed_kt" %in% names(df))) df$storm_speed_kt <- NA_real_

# --- Knaff & Zehr RMW with latitude dependence ---
  df <- df |>
    dplyr::mutate(
      R34_missing = !is.finite(.data$R34_km) | .data$R34_km <= 0,
      RMW_km = dplyr::case_when(
        # Prefer observed radii-derived RMW
        is.finite(.data$R64_km) & .data$R64_km > 0 ~ 0.35 * .data$R64_km,
        is.finite(.data$R50_km) & .data$R50_km > 0 ~ 0.40 * .data$R50_km,
        is.finite(.data$R34_km) & .data$R34_km > 0 & !.data$R34_missing ~ 0.50 * .data$R34_km,
        # Fallback: Knaff & Zehr (2007) with latitude
        TRUE ~ estimate_RMW_knaff(.data$Vmax_kt, .data$lat)
      )
    )

  stopifnot(nrow(df) == length(df$storm_speed_kt), nrow(df) == length(df$heading_deg))

  # --- Use patched Holland profile (vectorized; removes mapply bottleneck) ---
  df |>
    dplyr::mutate(
      V_site_symmetric_kt = .estimate_site_wind_holland(
        Vmax_kt = .data$Vmax_kt,
        r_km    = .data$dist_km,
        R34_km  = .data$R34_km,
        R50_km  = .data$R50_km,
        R64_km  = .data$R64_km,
        RMW_km  = .data$RMW_km,
        lat     = .data$lat,
        Pn      = 1013
      ),
      V_site_kt = .add_forward_motion_asymmetry(
        V_site_base_kt    = .data$V_site_symmetric_kt,
        storm_speed_kt    = .data$storm_speed_kt,
        r_km              = .data$dist_km,
        bearing_to_target = .data$bearing_to_target,
        storm_heading     = .data$heading_deg,
        RMW_km            = .data$RMW_km
      )
    )
}

# =============================================================================
# 4) Event classification and summarization
# =============================================================================

#' Classify storm class from peak site wind
#'
#' @param V_site_max_kt Numeric vector of peak site winds (kt).
#' @param ts_threshold_kt Threshold (kt) for Tropical Storm.
#' @param hurricane_threshold_kt Threshold (kt) for Hurricane.
#'
#' @return Character vector with values \code{c("TD", "TS", "HUR", "unknown")}.
#' @export
classify_severity <- function(V_site_max_kt,
                              ts_threshold_kt = 34,
                              hurricane_threshold_kt = 64) {
  v <- V_site_max_kt
  out <- rep("unknown", length(v))
  ok <- is.finite(v)
  if (!any(ok)) return(out)
  out[ok & v < ts_threshold_kt] <- "TD"
  out[ok & v >= ts_threshold_kt & v < hurricane_threshold_kt] <- "TS"
  out[ok & v >= hurricane_threshold_kt] <- "HUR"
  out
}



#' Aggregate track points into storm-level events
#'
#' @param track_df Track-point tibble/data.frame with at least SID, iso_time.
#'
#' @return Tibble with one row per storm and key event attributes.
#' @export
make_storm_events <- function(track_df) {
  if (!requireNamespace('lubridate', quietly = TRUE)) stop('Package `lubridate` is required.')

  df <- track_df

  if (!("V_site_kt" %in% names(df))) df$V_site_kt <- NA_real_
  if (!("wind_kt"   %in% names(df))) df$wind_kt   <- NA_real_
  if (!("pres_hpa"  %in% names(df))) df$pres_hpa  <- NA_real_
  if (!("poci_hpa"  %in% names(df))) df$poci_hpa  <- NA_real_
  if (!("rmw_km"    %in% names(df))) df$rmw_km    <- NA_real_

  df <- df |> dplyr::filter(!is.na(.data$iso_time))

  out <- df |>
    dplyr::mutate(
      dP_hpa = dplyr::if_else(
        is.finite(.data$poci_hpa) & is.finite(.data$pres_hpa),
        .data$poci_hpa - .data$pres_hpa,
        NA_real_
      )
    ) |>
    dplyr::group_by(.data$SID) |>
    dplyr::summarise(
      start_time = suppressWarnings(min(.data$iso_time, na.rm = TRUE)),
      end_time   = suppressWarnings(max(.data$iso_time, na.rm = TRUE)),
      n_points   = dplyr::n(),

      peak_wind_kt = suppressWarnings(max(.data$V_site_kt, na.rm = TRUE)),

      storm_intensity_kt = suppressWarnings(max(.data$wind_kt, na.rm = TRUE)),

      min_pressure_hpa  = suppressWarnings(min(.data$pres_hpa, na.rm = TRUE)),
      pressure_deficit_hpa  = suppressWarnings(max(.data$dP_hpa, na.rm = TRUE)),

      rmw_min_km  = suppressWarnings(min(.data$rmw_km, na.rm = TRUE)),
      rmw_mean_km = suppressWarnings(mean(.data$rmw_km, na.rm = TRUE)),

      .groups = "drop"
    ) |>
    dplyr::rename(storm_id = "SID") |>
    dplyr::mutate(
      peak_wind_kt = dplyr::if_else(is.finite(.data$peak_wind_kt), .data$peak_wind_kt, NA_real_),
      storm_intensity_kt   = dplyr::if_else(is.finite(.data$storm_intensity_kt),   .data$storm_intensity_kt,   NA_real_),
      min_pressure_hpa    = dplyr::if_else(is.finite(.data$min_pressure_hpa),    .data$min_pressure_hpa,    NA_real_),
      pressure_deficit_hpa    = dplyr::if_else(is.finite(.data$pressure_deficit_hpa),    .data$pressure_deficit_hpa,    NA_real_),
      rmw_min_km    = dplyr::if_else(is.finite(.data$rmw_min_km),    .data$rmw_min_km,    NA_real_),
      rmw_mean_km   = dplyr::if_else(is.finite(.data$rmw_mean_km),   .data$rmw_mean_km,   NA_real_),
      year = lubridate::year(.data$start_time)
    ) |>
    tibble::as_tibble()

  out
}

# =============================================================================
# 5) Rate model helpers (per-severity annual counts)
# =============================================================================

#' Compute annual counts of unique storm events by storm class
#'
#' @param events Tibble with at least columns: year, storm_class, storm_id.
#' @param severities Character vector of classes to include.
#'
#' @return Tibble with columns year, storm_class, n_events, completed to include
#'   all years and classes with zeros.
#'
#' @export
compute_annual_counts <- function(events, severities = c("TS", "HUR")) {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` is required.")

  events |>
    dplyr::filter(.data$storm_class %in% severities) |>
    dplyr::distinct(.data$year, .data$storm_class, .data$storm_id) |>
    dplyr::count(.data$year, .data$storm_class, name = "n_events") |>
    tidyr::complete(
      year = tidyr::full_seq(range(.data$year, na.rm = TRUE), 1),
      storm_class = severities,
      fill = list(n_events = 0)
    )
}

#' Compute Poisson rate table by storm class from annual counts
#'
#' @param annual_counts Output from \code{compute_annual_counts()}.
#' @return Tibble with lambda, n_years, prob_annual, prob_none by storm_class.
#' @export
compute_lambda_table <- function(annual_counts) {

  annual_counts |>
    dplyr::group_by(.data$storm_class) |>
    dplyr::summarise(
      lambda = mean(.data$n_events),
      n_years = dplyr::n(),
      prob_annual = 1 - exp(-.data$lambda),
      prob_none = exp(-.data$lambda),
      .groups = "drop"
    )
}

#' Derive annual counts from model output events
#'
#' @description
#' Builds annual storm counts from \code{out$events}, zero-filled for all
#' year Ã— storm_class combinations by location.
#'
#' @param out List returned by \code{run_hazard_model()}.
#'
#' @return Tibble with columns \code{location}, \code{year}, \code{storm_class},
#'   \code{n_events}.
#' @export
get_annual_counts <- function(out) {
  if (is.null(out$events)) stop("out$events is required.", call. = FALSE)
  events <- dplyr::as_tibble(out$events)
  if (nrow(events) == 0L) {
    return(tibble::tibble(
      location = character(0),
      year = integer(0),
      storm_class = character(0),
      n_events = integer(0)
    ))
  }
  classes <- sort(unique(stats::na.omit(events$storm_class)))
  events |>
    dplyr::filter(is.finite(.data$year), !is.na(.data$location), !is.na(.data$storm_class)) |>
    dplyr::group_by(.data$location) |>
    dplyr::group_modify(~ compute_annual_counts(.x, severities = classes)) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$location, .data$year, .data$storm_class)
}

#' Estimate overdispersion (k-hat) for a Poisson-Gamma annual activity factor
#'
#' @description
#' Uses total annual counts N = sum_severity n_events and the NegBin identity
#' Var(N) = mu + mu^2/k to estimate k. If Var <= mu, returns a large k (ÃƒÂ¢Ã¢â‚¬Â°Ã‹â€ Poisson).
#'
#' @param annual_counts Output from \code{compute_annual_counts()}.
#' @return List with k_hat, annual_total (year,N), mu, var.
#' @export
estimate_k_hat <- function(annual_counts) {

  annual_total <- annual_counts |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(N = sum(.data$n_events), .groups = "drop")

  mu <- mean(annual_total$N)
  va <- stats::var(annual_total$N)

  k_hat <- if (is.finite(va) && va > mu && mu > 0) mu^2 / (va - mu) else 1e6
  list(k_hat = k_hat, annual_total = annual_total, mu = mu, var = va)
}

