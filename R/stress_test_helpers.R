# Helpers for stress-test experiment scripts (inst/extcode/05-stress-test-*).
# All helpers are internal (dot-prefix, @noRd) and accessed from external
# scripts via ipdcstorm:::.<name>().

# .weibull_scale_from_mean ----------------------------------------------------

#' Solve Weibull scale from a target mean
#'
#' Returns the Weibull scale parameter that yields the requested mean wind
#' speed at a fixed shape parameter. Used by
#' `.make_stress_test_background_wind_cfg()` to derive per-location Weibull
#' marginals from literature-based mean trade-wind values.
#'
#' @param mean_kt Numeric. Target mean wind speed in knots.
#' @param shape Numeric. Weibull shape parameter.
#' @return Numeric. Weibull scale parameter.
#' @keywords internal
#' @noRd
.weibull_scale_from_mean <- function(mean_kt, shape) {
  mean_kt / gamma(1 + 1 / shape)
}

# .make_stress_test_background_wind_cfg ---------------------------------------

#' Build a background-wind config for the stress-test experiment
#'
#' Constructs a 12-month Weibull marginal per Caribbean target location
#' (St_Martin, Saba, Statia, Puerto_Rico, Miami) using literature-based mean
#' trade-wind values, then wraps the result in a background-wind config via
#' [make_background_wind_cfg()].
#'
#' @return A background-wind configuration list as returned by
#'   [make_background_wind_cfg()].
#' @keywords internal
#' @noRd
.make_stress_test_background_wind_cfg <- function() {
  bg_means <- tibble::tribble(
    ~location,      ~mean_kt, ~shape,
    "St_Martin",      9.5,      2.0,
    "Saba",          10.5,      2.0,
    "Statia",        10.0,      2.0,
    "Puerto_Rico",    8.5,      2.0,
    "Miami",          5.5,      2.2
  )

  weibull_params <- stats::setNames(
    lapply(seq_len(nrow(bg_means)), function(i) {
      tibble::tibble(
        month = 1:12,
        shape = bg_means$shape[i],
        scale = .weibull_scale_from_mean(bg_means$mean_kt[i], bg_means$shape[i])
      )
    }),
    bg_means$location
  )

  make_background_wind_cfg(
    weibull_params = weibull_params,
    ar1 = 0.4
  )
}

# .remap_replay_years ---------------------------------------------------------

#' Remap replay-slot sim_years onto baseline-year IDs
#'
#' Replay runs use sequential sim_year slots `1..length(baseline_ids)`. This
#' helper rewrites the `sim_year` column of each daily tibble to the
#' corresponding baseline-year ID, so scenario output is keyed identically
#' to the baseline for paired comparison.
#'
#' @param daily_list Named list of daily tibbles (one element per location).
#' @param baseline_ids Integer vector of baseline year IDs, aligned by
#'   position with replay slots `seq_along(baseline_ids)`.
#' @return A named list with the same shape as `daily_list`, with `sim_year`
#'   remapped.
#' @keywords internal
#' @noRd
.remap_replay_years <- function(daily_list, baseline_ids) {
  replay_slots <- seq_along(baseline_ids)
  lapply(daily_list, function(df) {
    out <- df
    slot_idx <- match(out$sim_year, replay_slots)
    out$sim_year <- baseline_ids[slot_idx]
    out
  })
}

# .build_pin_map --------------------------------------------------------------

#' Sample one focal SID per sim_year and return a pinned_sids map
#'
#' Returns the named-list format expected by the `pinned_sids` argument of
#' [generate_daily_hazard_impact_spatial()]. Seeding is local to this call,
#' so pool assignment is reproducible regardless of outside RNG state;
#' [run_hazard_model()] reseeds internally before drawing storm counts.
#'
#' @param sim_years Integer vector of sim_year slots to assign.
#' @param pool Character vector of candidate SIDs sampled with replacement.
#' @param seed Integer seed for the pool sampler.
#' @return A named list. Names are `sim_years` as character; values are the
#'   single SID assigned to each sim_year slot.
#' @keywords internal
#' @noRd
.build_pin_map <- function(sim_years, pool, seed) {
  if (length(pool) == 0L) {
    stop("`pool` must contain at least one SID.", call. = FALSE)
  }
  set.seed(seed)
  assigned <- sample(pool, length(sim_years), replace = TRUE)
  stats::setNames(as.list(assigned), as.character(sim_years))
}

# .build_pin_jitter -----------------------------------------------------------

#' Build per-year jitter specs for pinned focal events
#'
#' Returns a named list of jitter specifications, one per sim_year. Each
#' specification is a list with fields:
#'   - `doy_offset` : integer day offset applied to the historical DOY of the
#'     pinned SID (used in Gaussian mode).
#'   - `doy_abs`    : absolute day-of-year that overrides the historical DOY
#'     of the pinned SID (used in uniform-season mode; `NA` otherwise).
#'   - `v_scale`    : peak-wind multiplier.
#'   - `r_scale`    : RMW multiplier.
#'
#' Two timing modes:
#'
#'   - `"gaussian"` (default) - timing jitter is sampled from
#'     \eqn{N(0, \mathrm{doy\_sd})} and applied as an offset around the
#'     historical DOY. Preserves the climatological timing of each SID.
#'
#'   - `"uniform_season"` - timing is resampled uniformly within the Atlantic
#'     hurricane season (defaults: June 1 - November 30). The historical DOY
#'     of the pinned SID is discarded in favour of the drawn date. Each
#'     sim_year is therefore a distinct "what if Irma hit on August 10 /
#'     October 5 / November 15" realisation. Because the jitter map is
#'     reused verbatim by the scenario replays, the same year has the
#'     same focal date under baseline and KNMI scenarios.
#'
#' Scientific rationale:
#'   - `doy_sd` (Gaussian mode) captures within-peak-season timing
#'     variability; 7 days is ~1/3 of the Atlantic peak width.
#'   - `v_scale_sd` represents best-track peak-intensity uncertainty for
#'     Cat-4/5 storms (Torn & Snyder 2012, Landsea & Franklin 2013 report
#'     operational bias and sampling spread of ~5 percent at peak).
#'   - `r_scale_sd` represents uncertainty in estimated RMW (Vickery &
#'     Wadhera 2008 report ~15 percent scatter; 5 percent is a conservative
#'     internal-variability proxy).
#'   - Uniform-season mode assumes that conditional on being a Cat-4/5 event
#'     in the north-eastern Caribbean, timing within the season is roughly
#'     uniform in a stress-test setting. Narrower ranges (e.g., Aug 1 -
#'     Oct 31) can be configured for a peak-season-only stress test.
#'
#' Scales are clamped to `[0.5, Inf)` so the pinned event cannot be scaled
#' below half its catalogued intensity or radius.
#'
#' @param sim_years Integer vector of sim_year slots.
#' @param mode One of `"gaussian"` (default) or `"uniform_season"`.
#' @param doy_sd Numeric SD (days) for Gaussian timing jitter. Default 7.
#' @param doy_min Integer lower bound (day-of-year) for uniform-season
#'   sampling. Default 152 (June 1 in a non-leap year).
#' @param doy_max Integer upper bound (day-of-year) for uniform-season
#'   sampling. Default 334 (November 30).
#' @param v_scale_sd Numeric SD for the peak-wind multiplier. Default 0.05.
#' @param r_scale_sd Numeric SD for the RMW multiplier. Default 0.05.
#' @param seed Integer seed for the jitter sampler.
#' @return A named list. Names are `sim_years` as character; values are
#'   lists with fields `doy_offset`, `doy_abs`, `v_scale`, `r_scale`.
#' @keywords internal
#' @noRd
.build_pin_jitter <- function(sim_years,
                              mode = c("gaussian", "uniform_season"),
                              doy_sd = 7,
                              doy_min = 152L,
                              doy_max = 334L,
                              v_scale_sd = 0.05,
                              r_scale_sd = 0.05,
                              seed) {
  mode <- match.arg(mode)
  if (!is.numeric(doy_sd) || doy_sd < 0 ||
      !is.numeric(v_scale_sd) || v_scale_sd < 0 ||
      !is.numeric(r_scale_sd) || r_scale_sd < 0) {
    stop("jitter SDs must be non-negative numerics", call. = FALSE)
  }
  if (!is.numeric(doy_min) || !is.numeric(doy_max) ||
      doy_min < 1 || doy_max > 366 || doy_min >= doy_max) {
    stop("`doy_min` and `doy_max` must satisfy 1 <= doy_min < doy_max <= 366.",
         call. = FALSE)
  }

  set.seed(seed)
  n <- length(sim_years)

  if (mode == "uniform_season") {
    doy_abs <- as.integer(floor(stats::runif(n, doy_min, doy_max + 1)))
    doy_abs <- pmin(doy_max, doy_abs)
    doy_offsets <- rep(NA_integer_, n)
  } else {
    doy_offsets <- as.integer(round(stats::rnorm(n, 0, doy_sd)))
    doy_abs <- rep(NA_integer_, n)
  }

  v_scales <- pmax(0.5, stats::rnorm(n, 1, v_scale_sd))
  r_scales <- pmax(0.5, stats::rnorm(n, 1, r_scale_sd))

  stats::setNames(
    lapply(seq_len(n), function(i) {
      list(
        doy_offset = doy_offsets[i],
        doy_abs    = doy_abs[i],
        v_scale    = v_scales[i],
        r_scale    = r_scales[i]
      )
    }),
    as.character(sim_years)
  )
}

# .apply_state_dependent_damage ----------------------------------------------

#' Amplify daily damage by cumulative season damage state
#'
#' Post-processes the daily hazard tibble so that each day's damage rate is
#' amplified by the cumulative damage already accrued earlier in the same
#' sim_year at the same location. This captures the HAZUS-HM / progressive-
#' damage phenomenon: already-compromised structures have lower residual
#' capacity, so subsequent wind exposure causes disproportionate new damage.
#'
#' Formulation (simplified from HAZUS-HM sequential-damage module):
#'
#'     V(t) = sum_{s <= t} damage_rate_amplified(s)
#'     damage_rate_amplified(t) = damage_rate_raw(t) * (1 + alpha * min(V(t-1), cap))
#'
#' where `V(t)` is the running damage state at the end of day t. Amplification
#' uses `V(t-1)`, so intra-day feedback is excluded; intra-event feedback
#' across days IS included, which is physical: day-2 wind loads a roof
#' already weakened by day-1 debris and fatigue.
#'
#' References:
#'   - Vickery, P.J., Skerlj, P.F., Lin, J., Twisdale, L.A., Young, M.A.,
#'     Lavelle, F.M. (2006). HAZUS-MH Hurricane Model Methodology. FEMA.
#'   - Pinelli, J.P. et al. (2004). Hurricane damage prediction model for
#'     residential structures. J. Struct. Eng. 130(11): 1685-1691.
#'   - Grossi, P. & Kunreuther, H. (2005). Catastrophe Modeling: A New
#'     Approach to Managing Risk. Springer.
#'
#' @param daily_list Named list of daily tibbles (one element per location).
#'   Each tibble must have columns `sim_year`, `doy`, `damage_rate`,
#'   `cum_damage`.
#' @param alpha Numeric amplification coefficient. Literature range 2-5 for
#'   wind-driven structural damage. Default 3.
#' @param cap Numeric upper bound on the damage state used for amplification,
#'   preventing runaway amplification in saturated years. Default 0.5.
#' @return A named list with the same shape as `daily_list`, with
#'   `damage_rate` and `cum_damage` replaced by amplified values. Adds a
#'   column `damage_rate_raw` preserving the pre-amplification values.
#' @keywords internal
#' @noRd
.apply_state_dependent_damage <- function(daily_list,
                                          alpha = 3,
                                          cap = 0.5) {
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
    stop("`alpha` must be a single non-negative numeric.", call. = FALSE)
  }
  if (!is.numeric(cap) || length(cap) != 1L || cap <= 0 || cap > 1) {
    stop("`cap` must be a single numeric in (0, 1].", call. = FALSE)
  }

  amplify_series <- function(d_raw) {
    n <- length(d_raw)
    if (n == 0L) return(numeric(0))
    d_new <- numeric(n)
    V <- 0
    for (i in seq_len(n)) {
      d_new[i] <- d_raw[i] * (1 + alpha * min(cap, V))
      V <- V + d_new[i]
    }
    d_new
  }

  lapply(daily_list, function(df) {
    if (nrow(df) == 0L) return(df)
    # Sort once, then partition by sim_year so the per-year index lookup is
    # O(N) via split() instead of O(N) per year via which().
    df <- df[order(df$sim_year, df$doy), , drop = FALSE]
    df$damage_rate_raw <- df$damage_rate
    year_idx <- split(seq_len(nrow(df)), df$sim_year)
    for (idx in year_idx) {
      df$damage_rate[idx] <- amplify_series(df$damage_rate_raw[idx])
      df$cum_damage[idx]  <- cumsum(df$damage_rate[idx])
    }
    df
  })
}

# .verify_pins_landed ---------------------------------------------------------

#' Verify that pinned SIDs appear in daily output
#'
#' For each entry in `pin_map`, check that the daily output for that
#' `sim_year` contains at least one row whose `event_id` starts with the
#' pinned SID. Event ids use the format `"{SID}_y{cal_year}_{counter}"`, so
#' a prefix match identifies the pinned event regardless of the calendar-
#' year and counter suffix.
#'
#' @param daily_list Named list of daily tibbles (one element per location).
#' @param pin_map Named list mapping `sim_year` (as character) to the
#'   pinned SID for that year.
#' @return A data frame with columns `sim_year`, `sid`, `landed` for entries
#'   whose pin did not land. Empty if all pins landed.
#' @keywords internal
#' @noRd
.verify_pins_landed <- function(daily_list, pin_map) {
  all_rows <- dplyr::bind_rows(daily_list)
  checks <- lapply(names(pin_map), function(yr_key) {
    yr  <- as.integer(yr_key)
    sid <- pin_map[[yr_key]]
    year_rows <- all_rows[all_rows$sim_year == yr & !is.na(all_rows$event_id), ]
    landed <- any(startsWith(as.character(year_rows$event_id), sid))
    data.frame(sim_year = yr, sid = sid, landed = landed,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, checks)
  res[!res$landed, , drop = FALSE]
}
