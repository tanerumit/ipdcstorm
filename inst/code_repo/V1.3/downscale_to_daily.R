



#' Build an empirical event library for daily hazard simulation
#'
#' Constructs an event library for a single location from historical
#' storm-event outputs and trackpoint-level site winds. The library
#' contains (i) storm-level intensity and duration summaries by severity
#' class and (ii) an empirical seasonal distribution of storm timing
#' (day-of-year) at the site.
#'
#' This library is used to stochastically generate synthetic events
#' when downscaling annual hazard frequencies to daily time series.
#'
#' @param out List returned by \code{run_hazard_model()}, containing at least
#'   \code{events_all} and \code{trackpoints}.
#' @param island Character string giving the target island/location name.
#' @param severities Character vector of severity classes to include
#'   (default: \code{c("TS", "HUR")}).
#' @param thr_ts Numeric wind threshold (kt) used to identify the first
#'   site-impact timing for seasonality estimation (default: 34 kt).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{events}{Data frame of historical storm events at the site,
#'     including empirical durations (in days) and peak site wind speeds.}
#'   \item{doy}{Integer vector of day-of-year values representing the
#'     empirical seasonal distribution of storm impacts at the site.}
#' }
#'
#' @details
#' Event durations are derived from exposure hours (\code{hours_ge34},
#' \code{hours_ge50}, \code{hours_ge64}) reported by the hazard model and
#' converted to whole days. Seasonality is estimated from the first time
#' each storm exceeds \code{thr_ts} at the site, falling back to the first
#' available trackpoint if no exceedance occurs.
#'
#' @export
build_event_library <- function(out, island,
                                severities = c("TS", "HUR"),
                                thr_ts = 34) {
  stopifnot(is.list(out))
  
  ev <- out$events_all |>
    dplyr::filter(island == !!island, severity %in% severities) |>
    dplyr::mutate(
      dur_days_34 = pmax(1L, as.integer(ceiling(hours_ge34 / 24))),
      dur_days_50 = pmax(1L, as.integer(ceiling(hours_ge50 / 24))),
      dur_days_64 = pmax(1L, as.integer(ceiling(hours_ge64 / 24)))
    )
  
  tp <- out$trackpoints[[island]]
  if (is.null(tp)) stop("No trackpoints found for island: ", island)
  
  # For seasonality: day-of-year of first time V_site_kt >= thr_ts for each storm
  doy_tbl <- tp |>
    dplyr::filter(is.finite(V_site_kt), is.finite(ISO_TIME)) |>
    dplyr::group_by(SID) |>
    dplyr::summarise(
      first_time = {
        idx <- which(V_site_kt >= thr_ts)
        if (length(idx) > 0) min(ISO_TIME[idx]) else min(ISO_TIME)
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(doy = lubridate::yday(first_time)) |>
    dplyr::filter(is.finite(doy), doy >= 1, doy <= 366)
  
  list(
    events = ev,
    doy = doy_tbl$doy
  )
}

#' Sample synthetic storm events for a single year
#'
#' Randomly samples tropical storm and hurricane events for a given year
#' using an empirical event library. Each sampled event is assigned a
#' start date (based on observed seasonality), duration, and peak
#' site wind speed.
#'
#' @param lib Event library produced by \code{build_event_library()}.
#' @param year Integer calendar year for which events are generated.
#' @param n_ts Integer number of tropical storm (TS) events to sample.
#' @param n_hur Integer number of hurricane (64+ kt) events to sample.
#' @param use_duration Character string indicating which duration field
#'   to use (currently reserved for future extensions).
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A tibble with one row per sampled event and columns:
#' \describe{
#'   \item{severity}{Event severity class.}
#'   \item{start_date}{Calendar date of event onset.}
#'   \item{dur_days}{Event duration in whole days.}
#'   \item{V_peak}{Peak sustained wind speed at the site (kt).}
#' }
#'
#' @details
#' Event timing is sampled from the empirical day-of-year distribution
#' in the library. Peak wind speeds and durations are resampled from
#' observed site-level storm events of the same severity class.
#'
#' @export
sample_events_for_year <- function(lib, year,
                                   n_ts, n_hur,
                                   use_duration = c("dur_days_34", "dur_days_64"),
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  use_duration <- match.arg(use_duration)
  
  ev <- lib$events
  doy <- lib$doy
  if (length(doy) == 0) stop("No DOY seasonality sample available.")
  
  sample_one <- function(sev) {
    pool <- ev |> dplyr::filter(severity == sev)
    if (nrow(pool) == 0) stop("No events in library for severity: ", sev)
    
    row <- pool[sample.int(nrow(pool), 1), , drop = FALSE]
    
    doy0 <- sample(doy, 1)                           # seasonal timing
    start_date <- as.Date(sprintf("%d-01-01", year)) + (doy0 - 1)
    
    # duration choice
    dur_days <- if (sev == "HUR") row$dur_days_64 else row$dur_days_34
    dur_days <- max(1L, as.integer(dur_days))
    
    # peak choice
    V_peak <- as.numeric(row$V_site_max_kt)
    if (!is.finite(V_peak) || V_peak <= 0) V_peak <- if (sev == "HUR") 80 else 40
    
    list(
      severity = sev,
      start_date = start_date,
      dur_days = dur_days,
      V_peak = V_peak
    )
  }
  
  events <- c(
    replicate(n_ts,  sample_one("TS"),        simplify = FALSE),
    replicate(n_hur, sample_one("HUR"), simplify = FALSE)
  )
  
  tibble::tibble(
    severity   = vapply(events, `[[`, character(1), "severity"),
    start_date = as.Date(vapply(events, function(x) as.character(x$start_date), character(1))),
    dur_days   = as.integer(vapply(events, `[[`, numeric(1), "dur_days")),
    V_peak     = as.numeric(vapply(events, `[[`, numeric(1), "V_peak"))
  )
}

#' Generate a parametric wind pulse for a storm event
#'
#' Creates a smooth daily wind-speed profile for a storm event,
#' scaled to a specified peak wind speed and duration.
#'
#' @param dur_days Integer event duration in days.
#' @param V_peak Numeric peak sustained wind speed (kt).
#' @param shape Character string specifying the pulse shape:
#'   \code{"cosine"} (default) or \code{"triangle"}.
#'
#' @return Numeric vector of length \code{dur_days} giving daily
#'   sustained wind speeds (kt).
#'
#' @details
#' The cosine pulse produces a smooth rise and decay in intensity,
#' while the triangular pulse produces a linear ramp-up and ramp-down.
#' Both are deterministic given duration and peak wind.
#'
#' @export
event_pulse <- function(dur_days, V_peak, shape = c("cosine", "triangle")) {
  shape <- match.arg(shape)
  d <- as.integer(dur_days)
  if (d <= 0) return(numeric(0))
  
  t <- seq_len(d)
  if (shape == "triangle") {
    mid <- (d + 1) / 2
    w <- 1 - abs(t - mid) / mid
  } else {
    # cosine hump from 0..pi
    w <- sin(pi * (t - 0.5) / d)
  }
  pmax(0, V_peak * w)
}




#' Generate a daily wind time series for a single calendar year
#'
#' Converts a set of sampled storm events into a daily time series of
#' maximum sustained wind speeds at a site. Overlapping events are
#' combined conservatively by taking the daily maximum wind.
#'
#' @param year Integer calendar year.
#' @param sampled_events Tibble produced by \code{sample_events_for_year()}.
#' @param thr_port Optional numeric wind threshold (kt) for port disruption.
#' @param thr_infra Optional numeric wind threshold (kt) for infrastructure disruption.
#' @param pulse_shape Character string specifying the within-event wind
#'   profile shape (passed to \code{event_pulse()}).
#'
#' @return A tibble with one row per day and columns:
#' \describe{
#'   \item{date}{Calendar date.}
#'   \item{wind_kt}{Daily maximum sustained wind speed (kt).}
#'   \item{port_disrupt}{Logical indicator of port disruption (if threshold provided).}
#'   \item{infra_disrupt}{Logical indicator of infrastructure disruption (if threshold provided).}
#' }
#'
#' @details
#' This function does not introduce new stochasticity; all randomness
#' enters through event sampling. The resulting daily series is suitable
#' for direct coupling to system-dynamics or impact models.
#'
#' @export
generate_daily_year <- function(year, sampled_events,
                                thr_port = NA_real_,
                                thr_infra = NA_real_,
                                pulse_shape = "cosine") {
  start <- as.Date(sprintf("%d-01-01", year))
  end   <- as.Date(sprintf("%d-12-31", year))
  dates <- seq.Date(start, end, by = "day")
  wind  <- rep(0, length(dates))
  
  for (k in seq_len(nrow(sampled_events))) {
    s <- sampled_events$start_date[k]
    d <- sampled_events$dur_days[k]
    V <- sampled_events$V_peak[k]
    
    idx0 <- as.integer(s - start) + 1L
    idx1 <- idx0 + d - 1L
    if (idx1 < 1L || idx0 > length(dates)) next
    
    idx0_clip <- max(1L, idx0)
    idx1_clip <- min(length(dates), idx1)
    
    pulse <- event_pulse(d, V, shape = pulse_shape)
    
    # align pulse to clipped indices
    pulse_start <- idx0_clip - idx0 + 1L
    pulse_end   <- pulse_start + (idx1_clip - idx0_clip)
    wind[idx0_clip:idx1_clip] <- pmax(wind[idx0_clip:idx1_clip], pulse[pulse_start:pulse_end])
  }
  
  out <- tibble::tibble(date = dates, wind_kt = wind)
  
  # Optional disruption channels (binary daily flags)
  out <- out |>
    dplyr::mutate(
      port_disrupt  = if (is.finite(thr_port))  wind_kt >= thr_port  else NA,
      infra_disrupt = if (is.finite(thr_infra)) wind_kt >= thr_infra else NA
    )
  
  out
}

#' Generate daily hazard time series from a stochastic hazard model
#'
#' Produces synthetic daily wind time series for a single location by
#' downscaling annual tropical cyclone event counts into daily sequences.
#' The procedure preserves empirical seasonality, intensity, and duration
#' characteristics inferred from historical site-level storm impacts.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param island Character string giving the target island/location name.
#' @param sim_years Integer vector of synthetic simulation years to generate
#'   (default: \code{1:1000}).
#' @param year0 Integer base calendar year corresponding to \code{sim_year == 1}.
#' @param thr_port Optional numeric wind threshold (kt) for port disruption.
#' @param thr_infra Optional numeric wind threshold (kt) for infrastructure disruption.
#' @param seed Integer random seed for reproducibility.
#'
#' @return A tibble with daily resolution containing:
#' \describe{
#'   \item{island}{Location name.}
#'   \item{sim_year}{Synthetic simulation year index.}
#'   \item{date}{Calendar date.}
#'   \item{wind_kt}{Daily maximum sustained wind speed (kt).}
#'   \item{port_disrupt}{Logical port disruption indicator (if defined).}
#'   \item{infra_disrupt}{Logical infrastructure disruption indicator (if defined).}
#' }
#'
#' @details
#' Annual event frequencies are taken from the hazard simulation
#' (\code{out$sim_all}). Individual events are generated by resampling
#' from historical site-level storm characteristics. The output is
#' consistent with the annual hazard statistics by construction.
#'
#' @export
generate_daily_from_hazard <- function(out, island,
                                       sim_years = 1:1000,
                                       year0 = 2000,
                                       thr_port = NA_real_,
                                       thr_infra = NA_real_,
                                       seed = 1) {
  set.seed(seed)
  
  lib <- build_event_library(out, island = island)
  
  sim <- out$sim_all |>
    dplyr::filter(island == !!island, sim_year %in% sim_years)
  
  if (nrow(sim) == 0) stop("No sim years found for island in out$sim_all.")
  
  daily_list <- vector("list", nrow(sim))
  
  for (i in seq_len(nrow(sim))) {
    yr <- year0 + (sim$sim_year[i] - 1L)
    
    sampled <- sample_events_for_year(
      lib = lib,
      year = yr,
      n_ts = sim$n_TS[i],
      n_hur = sim$n_HUR[i]
    )
    
    daily_list[[i]] <- generate_daily_year(
      year = yr,
      sampled_events = sampled,
      thr_port = thr_port,
      thr_infra = thr_infra
    ) |>
      dplyr::mutate(sim_year = sim$sim_year[i], island = island)
  }
  
  dplyr::bind_rows(daily_list) |>
    dplyr::relocate(island, sim_year, date)
}



#' Map daily wind speed to bounded hazard intensity and damage rate
#'
#' @param daily tibble/data.frame with at least `date` and `wind_kt`.
#' @param V0 Wind threshold (kt) below which damage is zero (default 34).
#' @param V1 Wind (kt) at which intensity saturates to 1 (default 120).
#' @param p  Nonlinearity exponent (default 3).
#' @param dmax Maximum daily damage fraction at intensity 1 (default 0.02).
#' @param thr_port,thr_infra Optional disruption thresholds (kt) to compute service flags.
#'
#' @return Input table with added columns:
#'   - haz_intensity (0..1)
#'   - damage_rate (0..dmax)
#'   - damage_increment (same as damage_rate for daily timestep)
#'   - port_service / infra_service (if thresholds provided)
#'   - recovery_blocked (0/1)
add_damage_forcing <- function(daily,
                               V0 = 34, V1 = 120,
                               p = 3,
                               dmax = 0.02,
                               thr_port = NA_real_,
                               thr_infra = NA_real_) {
  
  stopifnot(is.data.frame(daily))
  if (!("wind_kt" %in% names(daily))) stop("daily must contain `wind_kt`.", call. = FALSE)
  
  clip01 <- function(x) pmax(0, pmin(1, x))
  
  x <- (daily$wind_kt - V0) / (V1 - V0)
  I <- clip01(x)^p
  
  d <- dmax * I
  
  out <- dplyr::as_tibble(daily) |>
    dplyr::mutate(
      haz_intensity   = I,
      damage_rate     = d,
      damage_increment = d,                 # daily timestep
      recovery_blocked = as.integer(wind_kt >= V0)
    )
  
  if (is.finite(thr_port)) {
    out <- out |>
      dplyr::mutate(port_service = 1 - as.integer(wind_kt >= thr_port))
  }
  if (is.finite(thr_infra)) {
    out <- out |>
      dplyr::mutate(infra_service = 1 - as.integer(wind_kt >= thr_infra))
  }
  
  out
}




damage_rate_from_wind <- function(wind_kt,
                                  thr = 34,
                                  V_ref = 80,
                                  d_ref = 0.03,
                                  p = 3,
                                  d_max = 0.10) {
  stopifnot(is.numeric(wind_kt))
  denom <- (V_ref - thr)
  if (!is.finite(denom) || denom <= 0) stop("V_ref must be > thr")
  
  x <- pmax(0, (wind_kt - thr) / denom)
  rate <- d_ref * (x ^ p)
  rate <- pmin(rate, d_max)
  rate[!is.finite(rate)] <- NA_real_
  rate
}


