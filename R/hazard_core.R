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

#' @title Compute great-circle distance from storm center to target.
#' @description Computes Haversine distance in kilometres from each storm track point to a fixed target coordinate.
#' @param lat Numeric vector of storm latitudes in decimal degrees.
#' @param lon Numeric vector of storm longitudes in decimal degrees.
#' @param t_lat Numeric scalar target latitude in decimal degrees.
#' @param t_lon Numeric scalar target longitude in decimal degrees.
#' @return Numeric vector of distances in kilometres with the same length as `lat` and `lon`.
#' @examples
#' dist_to_target(lat = c(18, 18.1), lon = c(-63, -63.1), t_lat = 18.05, t_lon = -63.05)
#' @seealso \code{\link{calculate_bearing}}
#' @importFrom geosphere distHaversine
#' @keywords internal
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

#' @title Compute bearing from storm center to target.
#' @description Computes great-circle initial bearing in degrees from each storm track point to a fixed target coordinate.
#' @param lat Numeric vector of storm latitudes in decimal degrees.
#' @param lon Numeric vector of storm longitudes in decimal degrees.
#' @param t_lat Numeric scalar target latitude in decimal degrees.
#' @param t_lon Numeric scalar target longitude in decimal degrees.
#' @return Numeric vector of bearings in degrees with the same length as `lat` and `lon`.
#' @examples
#' calculate_bearing(lat = c(18, 18.1), lon = c(-63, -63.1), t_lat = 18.05, t_lon = -63.05)
#' @seealso \code{\link{dist_to_target}}
#' @importFrom geosphere bearing
#' @keywords internal
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

#' @title Convert bearing to quadrant.
#' @description Maps a bearing in degrees to the corresponding meteorological quadrant.
#' @param bearing Numeric vector of bearings in degrees.
#' @return Character vector containing `"NE"`, `"SE"`, `"SW"`, `"NW"`, or `NA`.
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

#' @title Select directional wind radius.
#' @description Returns the quadrant-specific wind radius matching the bearing-derived quadrant.
#' @param quadrant Character vector of quadrant labels such as `"NE"`, `"SE"`, `"SW"`, or `"NW"`.
#' @param r_ne Numeric vector or scalar radius for the northeast quadrant.
#' @param r_se Numeric vector or scalar radius for the southeast quadrant.
#' @param r_sw Numeric vector or scalar radius for the southwest quadrant.
#' @param r_nw Numeric vector or scalar radius for the northwest quadrant.
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

#' @title Enforce monotone wind radii.
#' @description Adjusts directional radii so that `R64_km <= R50_km <= R34_km` wherever paired values are finite.
#' @param R34_km Numeric vector of 34-kt radii in kilometres.
#' @param R50_km Numeric vector of 50-kt radii in kilometres.
#' @param R64_km Numeric vector of 64-kt radii in kilometres.
#' @return List with elements `R34_km`, `R50_km`, and `R64_km` after monotonicity enforcement.
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

#' @title Estimate climatological R34.
#' @description Provides an empirical estimate of the 34-kt wind radius for storms without observed radii data using a North Atlantic climatological relationship.
#' @param Vmax_kt Numeric vector of maximum sustained wind speeds in knots.
#' @param lat Numeric vector of storm latitudes in decimal degrees north.
#' @return Numeric vector of estimated 34-kt wind radii in kilometres, with `NA` where `Vmax_kt < 34` or inputs are non-finite.
#' @examples
#' estimate_R34_climo(Vmax_kt = c(30, 50, 100), lat = c(18, 18, 20))
#' @seealso \code{\link{estimate_RMW_knaff}}, \code{\link{compute_site_winds_full}}
#' @note This is used as a fallback when directional 34-kt radii are unavailable.
#' @keywords internal
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

#' @title Return Holland outer-cutoff multipliers.
#' @description Returns the deterministic multipliers used to convert `R34_km` into the outer cutoff radius for the Holland profile.
#' @return Named list with `observed` and `climo` multipliers.
#' @keywords internal
.holland_outer_cutoff_multipliers <- function() {
  list(
    observed = 1.5,
    climo = 1.5
  )
}

#' @title Resolve Holland outer cutoff radius.
#' @description Computes the outer cutoff radius in kilometres used to taper winds beyond the resolved 34-kt radius.
#' @param R34_km Numeric vector of 34-kt wind radii in kilometres.
#' @param R34_is_fallback Logical vector indicating whether each `R34_km` value comes from climatological fallback rather than observations.
#' @return Numeric vector of outer cutoff radii in kilometres.
#' @keywords internal
.resolve_holland_outer_cutoff_km <- function(R34_km, R34_is_fallback = FALSE) {
  mult_cfg <- .holland_outer_cutoff_multipliers()
  n <- length(R34_km)
  R_outer_km <- rep(300, n)
  R34_is_fallback <- rep_len(R34_is_fallback, n)
  has_R34 <- is.finite(R34_km) & (R34_km > 0)

  if (any(has_R34)) {
    mult <- rep(mult_cfg$observed, sum(has_R34))
    mult[R34_is_fallback[has_R34]] <- mult_cfg$climo
    R_outer_km[has_R34] <- mult * R34_km[has_R34]
  }

  R_outer_km[!is.finite(R_outer_km)] <- 300
  R_outer_km <- pmax(0, R_outer_km)

  R_outer_km
}



# =============================================================================
# 2b) Knaff & Zehr RMW estimation (with latitude)
# =============================================================================

#' @title Validate observed radius of maximum wind.
#' @description Flags observed radius-of-maximum-wind values that fall inside the accepted physical range.
#' @param rmw_km Numeric vector of observed radius of maximum wind values in kilometres.
#' @return Logical vector indicating whether each value is finite and within the accepted bounds.
#' @keywords internal
.is_valid_observed_rmw_km <- function(rmw_km) {
  is.finite(rmw_km) & rmw_km > 5 & rmw_km < 150
}

#' @title Cap inferred radius of maximum wind.
#' @description Applies intensity-dependent lower and upper bounds to inferred radius-of-maximum-wind values.
#' @param rmw_km Numeric vector of inferred radius of maximum wind values in kilometres.
#' @param Vmax_kt Numeric vector of maximum sustained wind speeds in knots.
#' @return Numeric vector of bounded radius of maximum wind values in kilometres.
#' @keywords internal
.cap_inferred_rmw_km <- function(rmw_km, Vmax_kt) {
  out <- rmw_km
  ok <- is.finite(out)
  if (!any(ok)) return(out)

  Vmax_clamped <- Vmax_kt
  Vmax_clamped[!is.finite(Vmax_clamped)] <- 34
  Vmax_clamped <- pmax(34, pmin(185, Vmax_clamped))

  # Empirical upper guard from the packaged North Atlantic IBTrACS sample:
  # the 95th-percentile observed RMW decreases with intensity, so non-observed
  # estimates are capped on a simple monotone curve instead of a flat ceiling.
  max_km <- pmax(50, 140 - 0.75 * (Vmax_clamped - 34))
  out[ok] <- pmax(8, pmin(max_km[ok], out[ok]))
  out
}

#' @title Estimate radius of maximum wind from mean wind radii.
#' @description Infers the radius of maximum wind from storm-wide mean 64-, 50-, or 34-kt radii using fixed regression coefficients.
#' @param R64_mean_km Numeric vector of mean 64-kt radii in kilometres.
#' @param R50_mean_km Numeric vector of mean 50-kt radii in kilometres.
#' @param R34_mean_km Numeric vector of mean 34-kt radii in kilometres.
#' @return Numeric vector of inferred radius of maximum wind values in kilometres.
#' @keywords internal
.estimate_rmw_from_mean_radii <- function(R64_mean_km, R50_mean_km, R34_mean_km) {
  n <- max(length(R64_mean_km), length(R50_mean_km), length(R34_mean_km))
  out <- rep(NA_real_, n)

  has_r64 <- is.finite(R64_mean_km) & R64_mean_km >= 10 & R64_mean_km <= 305.58
  has_r50 <- is.finite(R50_mean_km) & R50_mean_km >= 10 & R50_mean_km <= 481.52
  has_r34 <- is.finite(R34_mean_km) & R34_mean_km >= 10 & R34_mean_km <= 888.96

  # Fixed no-intercept coefficients fit once from the packaged North Atlantic
  # IBTrACS sample (rows with valid USA_RMW, using storm-wide mean quadrant radii):
  #   R64 tier  = 0.6517550 * mean(R64)
  #   R50 tier  = 0.6676334 * mean(R50) when R64 is unavailable
  #   R34 tier  = 0.4106665 * mean(R34) when R64/R50 are unavailable
  out[has_r64] <- 0.6517550 * R64_mean_km[has_r64]

  use_r50 <- !has_r64 & has_r50
  out[use_r50] <- 0.6676334 * R50_mean_km[use_r50]

  use_r34 <- !has_r64 & !has_r50 & has_r34
  out[use_r34] <- 0.4106665 * R34_mean_km[use_r34]

  out
}

#' @title Resolve track-point radius of maximum wind.
#' @description Resolves the radius of maximum wind at each track point by prioritising valid observations, then radii-based estimates, then climatological estimates.
#' @param rmw_obs_km Numeric vector of observed radius of maximum wind values in kilometres.
#' @param R64_mean_km Numeric vector of mean 64-kt radii in kilometres.
#' @param R50_mean_km Numeric vector of mean 50-kt radii in kilometres.
#' @param R34_mean_km Numeric vector of mean 34-kt radii in kilometres.
#' @param Vmax_kt Numeric vector of maximum sustained wind speeds in knots.
#' @param lat Numeric vector of storm latitudes in decimal degrees north.
#' @return Numeric vector of resolved radius of maximum wind values in kilometres.
#' @keywords internal
.resolve_trackpoint_rmw_km <- function(rmw_obs_km,
                                       R64_mean_km,
                                       R50_mean_km,
                                       R34_mean_km,
                                       Vmax_kt,
                                       lat) {
  obs_ok <- .is_valid_observed_rmw_km(rmw_obs_km)

  rmw_radii_km <- .estimate_rmw_from_mean_radii(R64_mean_km, R50_mean_km, R34_mean_km)
  rmw_radii_km <- .cap_inferred_rmw_km(rmw_radii_km, Vmax_kt)

  rmw_knaff_km <- estimate_RMW_knaff(Vmax_kt, lat)

  out <- rmw_knaff_km
  use_radii <- is.finite(rmw_radii_km)
  out[use_radii] <- rmw_radii_km[use_radii]
  out[obs_ok] <- rmw_obs_km[obs_ok]
  out
}

#' @title Estimate climatological radius of maximum wind.
#' @description Estimates radius of maximum wind from storm intensity and latitude using a Knaff and Zehr climatological relationship.
#' @param Vmax_kt Numeric vector of maximum sustained wind speeds in knots.
#' @param lat Numeric vector of storm latitudes in decimal degrees north.
#' @return Numeric vector of estimated radius of maximum wind values in kilometres.
#' @examples
#' estimate_RMW_knaff(Vmax_kt = c(40, 80, 120), lat = c(15, 18, 22))
#' @seealso \code{\link{estimate_R34_climo}}, \code{\link{compute_site_winds_full}}
#' @keywords internal
#' @export
estimate_RMW_knaff <- function(Vmax_kt, lat = 18) {
  # Knaff & Zehr (2007) Eq. 1 (simplified, Atlantic):
  #   RMW_nm = 66.785 - 0.09102Ãƒâ€šÃ‚Â·Vmax + 1.0619Ãƒâ€šÃ‚Â·(lat - 25)
  # Valid for Vmax in [30, 185] kt

  Vmax_clamped <- pmax(30, pmin(185, Vmax_kt))
  RMW_nm <- 66.785 - 0.09102 * Vmax_clamped + 1.0619 * (lat - 25)

  RMW_km <- RMW_nm * 1.852
  RMW_km[!is.finite(Vmax_kt)] <- NA_real_
  .cap_inferred_rmw_km(RMW_km, Vmax_kt)
}



#' @title Estimate site wind with Holland profile.
#' @description Computes site wind speed from storm intensity, storm size, and site distance using a Holland-type radial wind profile with internal fallback logic for missing radii.
#' @param Vmax_kt Numeric vector of maximum sustained wind speeds in knots.
#' @param r_km Numeric vector of distances from storm centre to site in kilometres.
#' @param R34_km Numeric vector of 34-kt radii in kilometres.
#' @param R50_km Numeric vector of 50-kt radii in kilometres.
#' @param R64_km Numeric vector of 64-kt radii in kilometres.
#' @param RMW_km Numeric vector of radius of maximum wind values in kilometres.
#' @param Pn Numeric vector of ambient pressures in hPa.
#' @param Pc Numeric vector of central pressures in hPa.
#' @param lat Numeric vector of storm latitudes in decimal degrees north.
#' @return Numeric vector of estimated sustained site winds in knots.
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
  disable_r34_calibration <- isTRUE(getOption("ipdcstorm.disable_r34_calibration", FALSE))

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

  rmw_over_r34_cap <- 4.0
  clamp_rmw <- ok & R34_is_climo & is.finite(R34_eff) & (R34_eff > 0)
  if (any(clamp_rmw)) {
    RMW_km0[clamp_rmw] <- pmax(5, pmin(RMW_km0[clamp_rmw], R34_eff[clamp_rmw] / rmw_over_r34_cap))
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
  can_cal <- !disable_r34_calibration &
    is.finite(R34_eff[use]) &
    (R34_eff[use] > 0) &
    (r_km0[use] > RMW_km0[use] * 1.2)
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
      # structure). Reduced from 0.5 to 0.3 for tighter tail geometry.
      is_climo_u <- R34_is_climo[use][can_cal]
      blend <- is_climo_u & good
      if (any(blend)) {
        cal_factor[blend] <- 1 + 0.3 * (cal_factor[blend] - 1)
      }

      idx <- which(can_cal)
      V_site_kt[idx] <- V_site_kt[idx] * cal_factor
    }
  }

  # Outer cutoff is deterministic and pathway-specific.
  # Both observed and fallback/climatological pathways currently use 1.5x.
  R_outer <- .resolve_holland_outer_cutoff_km(
    R34_km = R34_eff[use],
    R34_is_fallback = R34_is_climo[use]
  )

  beyond <- r_km0[use] > R_outer
  if (any(beyond)) {
    # Exponential decay beyond R_outer
    V_site_kt[beyond] <- V_site_kt[beyond] * exp(-2 * (r_km0[use][beyond] - R_outer[beyond]) / R_outer[beyond])
  }

  out[use] <- pmax(0, pmin(Vmax_kt0[use], V_site_kt))
  out
}


#' @title Compute storm heading by track point.
#' @description Adds a `heading_deg` column based on great-circle motion between successive track points within each storm.
#' @param df Data frame with columns `SID`, `iso_time`, `lat`, and `lon`.
#' @return Data frame with the original columns plus numeric `heading_deg` in degrees clockwise from north.
#' @examples
#' df <- data.frame(
#'   SID = c("AL012000", "AL012000"),
#'   iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
#'   lat = c(18, 18.2),
#'   lon = c(-63, -63.2)
#' )
#' compute_storm_heading(df)
#' @seealso \code{\link{compute_site_winds_full}}
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

#' @title Add forward-motion asymmetry.
#' @description Applies a translation-speed asymmetry adjustment to a symmetric site wind estimate using the angle between storm heading and site bearing.
#' @param V_site_base_kt Numeric vector of symmetric site wind estimates in knots.
#' @param storm_speed_kt Numeric vector of storm translation speeds in knots.
#' @param r_km Numeric vector of distances from storm centre to site in kilometres.
#' @param bearing_to_target Numeric vector of bearings from storm centre to site in degrees.
#' @param storm_heading Numeric vector of storm motion headings in degrees.
#' @param RMW_km Numeric vector of radius of maximum wind values in kilometres.
#' @return Numeric vector of asymmetry-adjusted site winds in knots.
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

#' @title Compute site winds for storm track points.
#' @description Computes distance, bearing, storm-motion diagnostics, resolved wind radii, and site wind estimates for each track point relative to a fixed target location.
#' @param df Track-point data frame containing the fields required for wind-field reconstruction.
#' @param target_lat Numeric scalar target latitude in decimal degrees.
#' @param target_lon Numeric scalar target longitude in decimal degrees.
#' @return Data frame with the original track-point fields plus derived wind-field and site-wind columns.
#' @examples
#' df <- data.frame(
#'   SID = c("AL012000", "AL012000"),
#'   iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
#'   lat = c(18, 18.2),
#'   lon = c(-63, -63.2),
#'   dist_km = c(NA_real_, NA_real_),
#'   wind_kt = c(60, 65),
#'   rmw_km = c(25, 25),
#'   r34_ne_nm = c(60, 60), r34_se_nm = c(55, 55), r34_sw_nm = c(50, 50), r34_nw_nm = c(55, 55),
#'   r50_ne_nm = c(30, 30), r50_se_nm = c(25, 25), r50_sw_nm = c(20, 20), r50_nw_nm = c(25, 25),
#'   r64_ne_nm = c(15, 15), r64_se_nm = c(10, 10), r64_sw_nm = c(10, 10), r64_nw_nm = c(10, 10),
#'   storm_speed_kt = c(12, 12)
#' )
#' compute_site_winds_full(df, target_lat = 18.05, target_lon = -63.05)
#' @seealso \code{\link{compute_storm_heading}}, \code{\link{estimate_R34_climo}}, \code{\link{estimate_RMW_knaff}}
#' @note Input data must already contain the expected IBTrACS-derived fields used in the internal wind solver.
#' @keywords internal
#' @export
compute_site_winds_full <- function(df, target_lat, target_lon) {
  mean_radius_nm <- function(r_ne, r_se, r_sw, r_nw) {
    m <- cbind(r_ne, r_se, r_sw, r_nw)
    m[!is.finite(m)] <- NA_real_
    out <- rowMeans(m, na.rm = TRUE)
    out[rowSums(is.finite(m)) == 0] <- NA_real_
    out
  }

  if (!("rmw_km" %in% names(df))) df$rmw_km <- NA_real_

  df <- df |>
    dplyr::mutate(
      dist_km = dplyr::if_else(
        is.na(.data$dist_km),
        dist_to_target(.data$lat, .data$lon, target_lat, target_lon),
        .data$dist_km
      ),
      bearing_to_target = calculate_bearing(.data$lat, .data$lon, target_lat, target_lon),
      quadrant = .get_quadrant(.data$bearing_to_target),

      # R34 radii quality gate:
      # If fewer than 3 quadrants are present, directional radii are often unreliable.
      # Fall back to mean radius when >=2 quadrants exist; otherwise leave missing
      # and allow climo infill downstream.
      nq34 = rowSums(is.finite(cbind(.data$r34_ne_nm, .data$r34_se_nm, .data$r34_sw_nm, .data$r34_nw_nm))),
      R34_nm_dir = .get_directional_radius(.data$quadrant, .data$r34_ne_nm, .data$r34_se_nm, .data$r34_sw_nm, .data$r34_nw_nm),
      R34_nm_mean = mean_radius_nm(.data$r34_ne_nm, .data$r34_se_nm, .data$r34_sw_nm, .data$r34_nw_nm),
      R34_nm = dplyr::case_when(
        .data$nq34 >= 3 ~ .data$R34_nm_dir,
        .data$nq34 >= 2 ~ .data$R34_nm_mean,
        TRUE ~ NA_real_
      ),
      R50_nm = .get_directional_radius(.data$quadrant, .data$r50_ne_nm, .data$r50_se_nm, .data$r50_sw_nm, .data$r50_nw_nm),
      R64_nm = .get_directional_radius(.data$quadrant, .data$r64_ne_nm, .data$r64_se_nm, .data$r64_sw_nm, .data$r64_nw_nm),

      R34_km = .data$R34_nm * 1.852,
      R50_km = .data$R50_nm * 1.852,
      R64_km = .data$R64_nm * 1.852,
      R34_mean_km = .data$R34_nm_mean * 1.852,
      R50_mean_km = mean_radius_nm(.data$r50_ne_nm, .data$r50_se_nm, .data$r50_sw_nm, .data$r50_nw_nm) * 1.852,
      R64_mean_km = mean_radius_nm(.data$r64_ne_nm, .data$r64_se_nm, .data$r64_sw_nm, .data$r64_nw_nm) * 1.852,

      Vmax_kt = .data$wind_kt
    )


  mm <- .enforce_monotone_radii(df$R34_km, df$R50_km, df$R64_km)
  df <- df |>
    dplyr::mutate(R34_km = mm$R34_km, R50_km = mm$R50_km, R64_km = mm$R64_km)

  df <- compute_storm_heading(df)

  if (!("storm_speed_kt" %in% names(df))) df$storm_speed_kt <- NA_real_

  # RMW precedence is deterministic:
  # 1) observed USA_RMW parsed to `rmw_km` when valid
  # 2) calibrated mapping from storm-wide mean wind radii
  # 3) guarded Knaff fallback
  df <- df |>
    dplyr::mutate(
      R34_missing = !is.finite(.data$R34_km) | .data$R34_km <= 0,
      RMW_km = .resolve_trackpoint_rmw_km(
        rmw_obs_km = .data$rmw_km,
        R64_mean_km = .data$R64_mean_km,
        R50_mean_km = .data$R50_mean_km,
        R34_mean_km = .data$R34_mean_km,
        Vmax_kt = .data$Vmax_kt,
        lat = .data$lat
      )
    )

  # ---------------------------------------------------------------------------
  # Diagnostics/traceability + robustness:
  # Compute the effective R34 and the effective RMW actually used in the wind
  # solver, so event-max diagnostics can see whether the missing-radii clamp is
  # activated for tail-driving points.
  # This also ensures the clamp affects the wind solve even if internal solver
  # state is not reflected back into the returned tibble.
  # ---------------------------------------------------------------------------
  rmw_over_r34_cap <- 4.0
  df <- df |>
    dplyr::mutate(
      lat0 = dplyr::if_else(is.finite(.data$lat), .data$lat, 18),
      R34_is_climo = .data$R34_missing & is.finite(.data$Vmax_kt) & (.data$Vmax_kt >= 34),
      R34_eff_km = dplyr::if_else(.data$R34_is_climo, estimate_R34_climo(.data$Vmax_kt, lat = .data$lat0), .data$R34_km),
      RMW_used_km = pmax(5, pmin(200, .data$RMW_km)),
      RMW_used_km = dplyr::if_else(
        .data$R34_is_climo & is.finite(.data$R34_eff_km) & (.data$R34_eff_km > 0),
        pmax(5, pmin(.data$RMW_used_km, .data$R34_eff_km / rmw_over_r34_cap)),
        .data$RMW_used_km
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
        RMW_km  = .data$RMW_used_km,
        lat     = .data$lat,
        Pn      = 1013
      ),
      V_site_kt = .add_forward_motion_asymmetry(
        V_site_base_kt    = .data$V_site_symmetric_kt,
        storm_speed_kt    = .data$storm_speed_kt,
        r_km              = .data$dist_km,
        bearing_to_target = .data$bearing_to_target,
        storm_heading     = .data$heading_deg,
        RMW_km            = .data$RMW_used_km
      ),
      # Physical guard: final site wind should not exceed storm Vmax.
      # Symmetric Holland profile already caps at Vmax; apply same guard
      # after adding forward-motion asymmetry.
      V_site_kt = dplyr::if_else(
        is.finite(.data$V_site_kt) & is.finite(.data$Vmax_kt),
        pmin(.data$V_site_kt, .data$Vmax_kt),
        .data$V_site_kt
      )
    )
}

# =============================================================================
# 4) Event classification and summarization
# =============================================================================

#' @title Classify storm class from peak site wind.
#' @description Classifies peak site wind into tropical depression, tropical storm, hurricane, or unknown.
#' @param V_site_max_kt Numeric vector of peak site wind speeds in knots.
#' @param ts_threshold_kt Numeric scalar tropical-storm threshold in knots.
#' @param hurricane_threshold_kt Numeric scalar hurricane threshold in knots.
#' @return Character vector containing `"TD"`, `"TS"`, `"HUR"`, or `"unknown"`.
#' @examples
#' classify_severity(c(20, 40, 80, NA_real_))
#' @seealso \code{\link{make_storm_events}}, \code{\link{compute_annual_counts}}
#' @keywords internal
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



#' @title Aggregate track points into storm events.
#' @description Summarises track-point data to one row per storm with peak site wind, storm intensity, pressure, timing, and radius-of-maximum-wind diagnostics.
#' @param track_df Track-point tibble or data frame with at least `SID` and `iso_time`, plus optional wind and pressure fields.
#' @return Tibble with one row per storm and event-level summary attributes.
#' @examples
#' track_df <- data.frame(
#'   SID = c("AL012000", "AL012000"),
#'   iso_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-08-01 06:00:00"), tz = "UTC"),
#'   V_site_kt = c(40, 55),
#'   wind_kt = c(60, 65)
#' )
#' make_storm_events(track_df)
#' @seealso \code{\link{classify_severity}}, \code{\link{compute_annual_counts}}
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

#' @title Normalise storm-class selection.
#' @description Resolves a user-supplied storm-class vector to a unique character vector with defaults when missing.
#' @param storm_classes Character vector of storm classes to retain.
#' @return Character vector of unique storm classes.
#' @keywords internal
.normalize_storm_classes <- function(storm_classes = NULL) {
  if (is.null(storm_classes)) {
    storm_classes <- c("TS", "HUR")
  }
  unique(as.character(storm_classes))
}

#' @title Compute annual storm-event counts.
#' @description Counts unique storm events per year and storm class, then completes the series so missing year-class combinations are filled with zeros.
#' @param events Tibble with at least `year`, `storm_class`, and `storm_id` columns.
#' @param storm_classes Character vector of storm classes to include.
#' @return Tibble with columns `year`, `storm_class`, and `n_events`, completed over all represented years and classes.
#' @examples
#' events <- data.frame(
#'   year = c(2000, 2000, 2001),
#'   storm_class = c("TS", "TS", "HUR"),
#'   storm_id = c("A", "B", "C")
#' )
#' compute_annual_counts(events)
#' @seealso \code{\link{compute_lambda_table}}, \code{\link{get_annual_counts}}
#' @keywords internal
#' @export
compute_annual_counts <- function(events,
                                  storm_classes = c("TS", "HUR")) {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` is required.")
  storm_classes <- .normalize_storm_classes(storm_classes = storm_classes)

  events |>
    dplyr::filter(.data$storm_class %in% storm_classes) |>
    dplyr::distinct(.data$year, .data$storm_class, .data$storm_id) |>
    dplyr::count(.data$year, .data$storm_class, name = "n_events") |>
    tidyr::complete(
      year = tidyr::full_seq(range(.data$year, na.rm = TRUE), 1),
      storm_class = storm_classes,
      fill = list(n_events = 0)
    )
}

#' @title Compute Poisson rate table by storm class.
#' @description Summarises annual event counts into storm-class-specific Poisson rates and annual exceedance probabilities.
#' @param annual_counts Tibble returned by \code{\link{compute_annual_counts}}.
#' @return Tibble with `storm_class`, `lambda`, `n_years`, `prob_annual`, and `prob_none` columns.
#' @examples
#' annual_counts <- data.frame(
#'   year = c(2000, 2001, 2000, 2001),
#'   storm_class = c("TS", "TS", "HUR", "HUR"),
#'   n_events = c(1, 0, 0, 1)
#' )
#' compute_lambda_table(annual_counts)
#' @seealso \code{\link{compute_annual_counts}}, \code{\link{estimate_k_hat}}
#' @keywords internal
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

#' @title Derive annual counts from model output.
#' @description Builds location-specific annual storm counts from `out$events` and zero-fills missing year-class combinations within each location.
#' @param out List returned by `run_hazard_model()` containing an `events` component.
#' @return Tibble with columns `location`, `year`, `storm_class`, and `n_events`.
#' @examples
#' out <- list(events = data.frame(
#'   location = c("Saba", "Saba"),
#'   year = c(2000, 2001),
#'   storm_class = c("TS", "HUR"),
#'   storm_id = c("A", "B")
#' ))
#' get_annual_counts(out)
#' @seealso \code{\link{compute_annual_counts}}, \code{\link{compute_lambda_table}}
#' @keywords internal
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
    dplyr::group_modify(~ compute_annual_counts(.x, storm_classes = classes)) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$location, .data$year, .data$storm_class)
}

#' @title Estimate annual-count overdispersion.
#' @description Estimates the negative-binomial overdispersion parameter from total annual storm counts using the identity `Var(N) = mu + mu^2 / k`.
#' @param annual_counts Tibble returned by \code{\link{compute_annual_counts}}.
#' @return List with elements `k_hat`, `annual_total`, `mu`, and `var`.
#' @examples
#' annual_counts <- data.frame(
#'   year = c(2000, 2000, 2001, 2001),
#'   storm_class = c("TS", "HUR", "TS", "HUR"),
#'   n_events = c(1, 0, 2, 1)
#' )
#' estimate_k_hat(annual_counts)
#' @seealso \code{\link{compute_annual_counts}}, \code{\link{compute_lambda_table}}
#' @keywords internal
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
