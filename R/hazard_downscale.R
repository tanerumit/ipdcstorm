
# =============================================================================
# Script overview: temporal downscaling and impact forcing
# - .assign_severity_simple(): simple TD/TS/HUR class from peak wind.
# - build_event_library(): empirical seasonality + stratified event resampling.
# - build_event_library_from_out(): convenience wrapper using run_hazard_model() output.
# - sample_events_for_year(): sample TS/HUR events with dates, duration, and peak wind.
# - event_pulse(): deterministic within-event daily wind profile.
# - generate_daily_year(): daily wind time series for one year.
# - generate_daily_from_hazard(): daily series from stochastic hazard output.
# - sample_events_for_year_extended(): sample events with additional attributes.
# - generate_daily_year_extended(): daily wind + dominant event attributes.
# - generate_daily_hazard_impact(): daily wind + event + damage forcing series
#     (with Level 3 storm perturbation integration via perturb_event()).
# - add_damage_forcing(): intensity + damage rate from daily wind.
# - damage_rate_from_wind(): bounded power-law damage rate function.
# =============================================================================

# =============================================================================
# 9) Temporal downscaling: event library + sampling + daily series
# =============================================================================

#' Assign simple severity class from peak wind
#'
#' @param wind_max_kt Numeric; peak wind (kt).
#' @return Character scalar in {"TD","TS","HUR"} (or NA).
#' @keywords internal
.assign_severity_simple <- function(wind_max_kt) {
  dplyr::case_when(
    !is.finite(wind_max_kt) ~ NA_character_,
    wind_max_kt >= 64 ~ "HUR",
    wind_max_kt >= 34 ~ "TS",
    TRUE ~ "TD"
  )
}


#' Build an empirical event library for resampling (seasonality + stratified events)
#'
#' @description
#' Constructs (1) empirical day-of-year (DOY) samples by class and
#' (2) stratified bins of historical events for resampling.
#'
#' @param track_df Track-point tibble for a single location, with columns `SID` and `iso_time`.
#' @param event_df Storm-event tibble (one row per SID) with at least `SID`.
#'   Recommended: `peak_wind_kt`, `storm_intensity_kt`, `min_pressure_hpa`,
#'   `rmw_mean_km`, `start_time`.
#' @param sev_levels Character vector of class levels.
#' @param bins List of break vectors for stratification (wind, Pc, RMW).
#' @param seed Optional integer seed.
#' @param resampling_method Character; one of {"stratified","copula_nn"}.
#'   Default "stratified" preserves existing behavior.
#' @param copula_min_n Integer; minimum complete-case sample size per class needed to fit copula.
#' @param copula_k Integer; number of nearest neighbours to sample from after copula proposal (>=1).
#' @param copula_robust_scale Logical; if TRUE use median/MAD scaling for NN distance, else mean/SD.
#'
#' @return A list with DOY samples, stratification bins, the event table, and sampler functions.
#' @export
# =============================================================================
# Drop-in replacement: build_event_library()
# Patches:
# - FIX silent corruption: factor/character severity compare in strat sampler
# - FIX silent corruption: guard match(SID) -> NA row (hard error instead)
# - FIX robustness: drop NA SIDs inside strat bins (prevents NA sampling)
# - FIX consistency: dur_days (copula_nn) computed as inclusive day count (>=1)
# =============================================================================
build_event_library <- function(track_df, event_df,
                                sev_levels = c("TD", "TS", "HUR"),
                                bins = list(
                                  wind = c(0, 34, 64, 83, 96, 113, Inf),
                                  Pc   = c(850, 900, 940, 970, 1000, 1050),
                                  RMW  = c(0, 10, 20, 30, 40, 60, Inf)
                                ),
                                seed = NULL,
                                resampling_method = c("stratified", "copula_nn"),
                                copula_min_n = 30L,
                                copula_k = 1L,
                                copula_robust_scale = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package `lubridate` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  
  resampling_method <- match.arg(resampling_method)
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(all(c("SID", "iso_time") %in% names(track_df)))
  if (!("SID" %in% names(event_df)) && ("storm_id" %in% names(event_df))) {
    event_df <- dplyr::mutate(event_df, SID = .data$storm_id)
  }
  stopifnot("SID" %in% names(event_df))
  
  ev <- tibble::as_tibble(event_df)
  
  if (!("V_site_max_kt" %in% names(ev)) && ("peak_wind_kt" %in% names(ev))) ev$V_site_max_kt <- ev$peak_wind_kt
  if (!("V_site_max_kt" %in% names(ev))) ev$V_site_max_kt <- NA_real_
  if (!("wind_max_kt"   %in% names(ev)) && ("storm_intensity_kt" %in% names(ev))) ev$wind_max_kt <- ev$storm_intensity_kt
  if (!("wind_max_kt"   %in% names(ev))) ev$wind_max_kt   <- NA_real_
  if (!("Pc_min_hPa"    %in% names(ev)) && ("min_pressure_hpa" %in% names(ev))) ev$Pc_min_hPa <- ev$min_pressure_hpa
  if (!("Pc_min_hPa"    %in% names(ev))) ev$Pc_min_hPa    <- NA_real_
  if (!("dP_max_hPa"    %in% names(ev)) && ("pressure_deficit_hpa" %in% names(ev))) ev$dP_max_hPa <- ev$pressure_deficit_hpa
  if (!("dP_max_hPa"    %in% names(ev))) ev$dP_max_hPa    <- NA_real_
  if (!("RMW_mean_km"   %in% names(ev)) && ("rmw_mean_km" %in% names(ev))) ev$RMW_mean_km <- ev$rmw_mean_km
  if (!("RMW_mean_km"   %in% names(ev))) ev$RMW_mean_km   <- NA_real_
  
  if (!("start_time" %in% names(ev))) {
    starts <- track_df |>
      dplyr::filter(!is.na(.data$iso_time)) |>
      dplyr::arrange(.data$SID, .data$iso_time) |>
      dplyr::group_by(.data$SID) |>
      dplyr::summarise(start_time = min(.data$iso_time), .groups = "drop")
    
    ev <- ev |>
      dplyr::left_join(starts, by = "SID")
  }
  
  # Only needed for copula_nn; kept inert otherwise
  if (resampling_method == "copula_nn" && !("dur_days" %in% names(ev))) {
    durs <- track_df |>
      dplyr::filter(!is.na(.data$iso_time)) |>
      dplyr::arrange(.data$SID, .data$iso_time) |>
      dplyr::group_by(.data$SID) |>
      dplyr::summarise(
        # inclusive day-count, min 1 (consistent with downstream duration logic)
        dur_days = max(
          1L,
          as.integer(floor(as.numeric(difftime(max(.data$iso_time), min(.data$iso_time), units = "days"))) + 1L)
        ),
        .groups = "drop"
      )
    ev <- ev |>
      dplyr::left_join(durs, by = "SID")
  }
  
  if (!("doy" %in% names(ev))) {
    ev <- ev |>
      dplyr::mutate(doy = lubridate::yday(.data$start_time))
  }
  
  ev <- ev |>
    dplyr::filter(!is.na(.data$start_time))
  
  if (nrow(ev) == 0) {
    stop("No events left after requiring non-missing start_time. Cannot build library.")
  }
  
  ev <- ev |>
    dplyr::mutate(
      wind_for_sev = dplyr::coalesce(.data$V_site_max_kt, .data$wind_max_kt),
      severity = .assign_severity_simple(.data$wind_for_sev),
      severity = factor(.data$severity, levels = sev_levels),
      doy = as.integer(.data$doy),
      
      wind_bin = cut(.data$wind_for_sev, breaks = bins$wind, include.lowest = TRUE, right = FALSE),
      Pc_bin   = cut(.data$Pc_min_hPa,   breaks = bins$Pc,   include.lowest = TRUE, right = FALSE),
      RMW_bin  = cut(.data$RMW_mean_km,  breaks = bins$RMW,  include.lowest = TRUE, right = FALSE)
    )
  
  doy_by_sev <- ev |>
    dplyr::filter(!is.na(.data$severity), is.finite(.data$doy)) |>
    dplyr::group_by(.data$severity) |>
    dplyr::summarise(doy = list(.data$doy), .groups = "drop")
  
  strat <- ev |>
    dplyr::filter(!is.na(.data$severity)) |>
    dplyr::mutate(
      wind_bin = as.character(.data$wind_bin),
      Pc_bin   = as.character(.data$Pc_bin),
      RMW_bin  = as.character(.data$RMW_bin)
    ) |>
    dplyr::group_by(.data$severity, .data$wind_bin, .data$Pc_bin, .data$RMW_bin) |>
    dplyr::summarise(
      n = dplyr::n(),
      # drop NA SIDs to prevent NA sampling and match() NA corruption
      sid = list(stats::na.omit(.data$SID)),
      .groups = "drop"
    ) |>
    dplyr::filter(lengths(.data$sid) > 0) |>
    dplyr::group_by(.data$severity) |>
    dplyr::mutate(w = .data$n / sum(.data$n)) |>
    dplyr::ungroup()
  
  if (nrow(strat) == 0) {
    stop("No stratification bins created. Check severity levels and available event data.")
  }
  
  # ---------------------------------------------------------------------------
  # copula_nn helpers (local, base R only; inert unless resampling_method="copula_nn")
  # ---------------------------------------------------------------------------
  .rank_to_norm <- function(x) {
    ok <- is.finite(x)
    z <- rep(NA_real_, length(x))
    if (!any(ok)) return(z)
    r <- rank(x[ok], ties.method = "average")
    u <- (r - 0.5) / length(r)
    z[ok] <- qnorm(u)
    z
  }
  
  .norm_to_emp_quantile <- function(u, x_ref) {
    xr <- x_ref[is.finite(x_ref)]
    if (length(xr) == 0) return(NA_real_)
    u <- pmin(1 - 1e-12, pmax(1e-12, u))
    as.numeric(stats::quantile(xr, probs = u, names = FALSE, type = 8))
  }
  
  .robust_scale_vec <- function(x) {
    xok <- x[is.finite(x)]
    if (length(xok) < 2) return(list(center = NA_real_, scale = NA_real_))
    if (isTRUE(copula_robust_scale)) {
      cen <- stats::median(xok)
      sc <- stats::mad(xok, center = cen, constant = 1, na.rm = TRUE)
    } else {
      cen <- mean(xok)
      sc <- stats::sd(xok)
    }
    if (!is.finite(sc) || sc <= 0) sc <- stats::sd(xok)
    if (!is.finite(sc) || sc <= 0) sc <- 1
    list(center = cen, scale = max(sc, 1e-12))
  }
  
  .standardize_matrix <- function(X) {
    p <- ncol(X)
    Z <- X
    cen <- rep(NA_real_, p)
    sc  <- rep(NA_real_, p)
    for (j in seq_len(p)) {
      ss <- .robust_scale_vec(X[, j])
      cen[j] <- ss$center
      sc[j]  <- ss$scale
      Z[, j] <- (X[, j] - cen[j]) / sc[j]
    }
    list(Z = Z, center = cen, scale = sc)
  }
  
  .fit_copula_by_sev <- function(ev_sub) {
    feat <- c("wind_for_sev", "Pc_min_hPa", "RMW_mean_km", "dur_days")
    feat <- feat[feat %in% names(ev_sub)]
    if (length(feat) < 2) return(NULL)
    
    X <- as.data.frame(ev_sub[, feat, drop = FALSE])
    cc <- stats::complete.cases(X)
    if (sum(cc) < as.integer(copula_min_n)) return(NULL)
    
    Xc <- as.matrix(X[cc, , drop = FALSE])
    Zc <- apply(Xc, 2, .rank_to_norm)
    if (any(colSums(is.finite(Zc)) < as.integer(copula_min_n))) return(NULL)
    
    R <- stats::cor(Zc, use = "pairwise.complete.obs")
    
    eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(eig)) || min(eig) <= 1e-10) {
      add <- max(1e-8, 1e-6 - min(eig, na.rm = TRUE))
      R <- R + diag(add, nrow(R))
      d <- sqrt(diag(R))
      R <- R / (d %o% d)
    }
    
    L <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(L)) return(NULL)
    
    std <- .standardize_matrix(as.matrix(ev_sub[, feat, drop = FALSE]))
    
    list(
      feat = feat,
      L = L,
      ev_sub = ev_sub,
      X = as.matrix(ev_sub[, feat, drop = FALSE]),
      Zstd = std$Z
    )
  }
  
  .sample_event_copula_nn <- function(fit) {
    p <- length(fit$feat)
    z <- as.numeric(crossprod(fit$L, rnorm(p)))
    u <- stats::pnorm(z)
    
    x_star <- rep(NA_real_, p)
    for (j in seq_len(p)) {
      x_star[j] <- .norm_to_emp_quantile(u[j], fit$X[, j])
    }
    
    X <- fit$X
    x_std <- rep(NA_real_, p)
    for (j in seq_len(p)) {
      ss <- .robust_scale_vec(X[, j])
      x_std[j] <- (x_star[j] - ss$center) / ss$scale
    }
    
    ok_dim <- is.finite(x_std)
    if (!any(ok_dim)) {
      i <- sample.int(nrow(fit$ev_sub), 1)
      return(fit$ev_sub[i, , drop = FALSE])
    }
    
    Z <- fit$Zstd[, ok_dim, drop = FALSE]
    dx <- sweep(Z, 2, x_std[ok_dim], FUN = "-")
    d2 <- rowSums(dx * dx)
    d2[!is.finite(d2)] <- Inf
    
    k <- max(1L, as.integer(copula_k))
    ord <- order(d2)
    ord <- ord[is.finite(d2[ord])]
    if (length(ord) == 0) {
      i <- sample.int(nrow(fit$ev_sub), 1)
      return(fit$ev_sub[i, , drop = FALSE])
    }
    top <- ord[seq_len(min(k, length(ord)))]
    i <- sample(top, 1)
    fit$ev_sub[i, , drop = FALSE]
  }
  
  copula_fits <- NULL
  if (resampling_method == "copula_nn") {
    copula_fits <- list()
    for (sev in sev_levels) {
      ev_sub <- ev[as.character(ev$severity) == sev, , drop = FALSE]
      if (nrow(ev_sub) > 0) {
        copula_fits[[sev]] <- .fit_copula_by_sev(ev_sub)
      } else {
        copula_fits[[sev]] <- NULL
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Samplers
  # ---------------------------------------------------------------------------
  sample_doy <- function(sev) {
    sev <- as.character(sev)
    idx <- match(sev, as.character(doy_by_sev$severity))
    if (is.na(idx)) stop("No DOY data for severity: ", sev)
    v <- doy_by_sev$doy[[idx]]
    if (length(v) == 0) stop("No DOY data for severity: ", sev)
    sample(v, size = 1)
  }
  
  sample_event <- function(sev) {
    sev <- as.character(sev)
    
    if (resampling_method == "copula_nn") {
      fit <- copula_fits[[sev]]
      if (!is.null(fit)) {
        return(.sample_event_copula_nn(fit))
      }
      # fall back to stratified
    }
    
    # Stratified behavior (patched)
    sub <- strat[as.character(strat$severity) == sev, , drop = FALSE]
    if (nrow(sub) == 0) stop("No events for severity: ", sev)
    
    p <- sub$w
    p[!is.finite(p) | p < 0] <- 0
    if (sum(p) == 0) p <- rep(1 / nrow(sub), nrow(sub)) else p <- p / sum(p)
    
    k <- sample.int(nrow(sub), size = 1, prob = p)
    sid <- sample(sub$sid[[k]], size = 1)
    if (is.na(sid)) stop("Sampled SID is NA (strat bin contains NA SIDs).", call. = FALSE)
    
    idx <- match(sid, ev$SID)
    if (is.na(idx)) stop("Sampled SID not found in event table: ", sid, call. = FALSE)
    ev[idx, , drop = FALSE]
  }
  
  list(
    doy_by_severity = doy_by_sev,
    strat_bins = strat,
    events = ev,
    sample_doy = sample_doy,
    sample_event = sample_event
  )
}

#' Compute disruption flags from daily hazard output
#'
#' @param daily Tibble from \code{generate_daily_hazard_impact()}.
#' @param thr_port Wind threshold (kt) for port disruption.
#' @param thr_infra Wind threshold (kt) for infrastructure disruption.
#' @param thr_surge Surge threshold (m) for surge disruption.
#' @param use_gust Logical; if TRUE, use \code{wind_gust_kt} for port/infra thresholds.
#'
#' @return Input tibble with additional logical columns:
#'   \code{port_disrupt}, \code{infra_disrupt}, \code{surge_disrupt}.
#' @export
disruption_flags <- function(daily,
                             thr_port = NA_real_,
                             thr_infra = NA_real_,
                             thr_surge = NA_real_,
                             use_gust = FALSE) {
  wind_col <- if (use_gust) "wind_gust_kt" else "wind_kt"
  out <- daily
  if (is.finite(thr_port)) out$port_disrupt <- out[[wind_col]] >= thr_port
  if (is.finite(thr_infra)) out$infra_disrupt <- out[[wind_col]] >= thr_infra
  if (is.finite(thr_surge)) out$surge_disrupt <- out$surge_m >= thr_surge
  out
}

#' Flag tropical cyclone or hurricane days
#'
#' @param daily Tibble from \code{generate_daily_hazard_impact()}.
#'
#' @return Logical vector.
#' @export
is_tc_day <- function(daily) {
  daily$event_class %in% c("TC", "HUR")
}

#' @rdname is_tc_day
#' @export
is_hur_day <- function(daily) {
  daily$event_class == "HUR"
}

#' Compute daily exposure hours above a wind threshold
#'
#' @param daily Tibble from \code{generate_daily_hazard_impact()}.
#' @param threshold_kt Wind speed threshold (kt).
#' @param use_gust Logical; if TRUE, use \code{wind_gust_kt}.
#'
#' @return Numeric vector of hours (0 or 24 at daily resolution).
#' @export
exposure_hours <- function(daily, threshold_kt, use_gust = FALSE) {
  wind_col <- if (use_gust) "wind_gust_kt" else "wind_kt"
  24 * (daily[[wind_col]] >= threshold_kt)
}

#' Summarise peak wind per simulation year
#'
#' @param daily Tibble from \code{generate_daily_hazard_impact()}.
#'
#' @return Tibble with \code{location}, \code{sim_year}, \code{scenario}, \code{peak_wind_kt}.
#' @export
peak_wind_by_year <- function(daily) {
  daily |>
    dplyr::group_by(.data$location, .data$sim_year, .data$scenario) |>
    dplyr::summarise(peak_wind_kt = max(.data$wind_kt, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(peak_wind_kt = dplyr::if_else(is.finite(.data$peak_wind_kt), .data$peak_wind_kt, NA_real_))
}









#' Convenience wrapper: build an event library from run_hazard_model() output
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param location Character; target location name (must exist in out lists).
#' @param ... Passed to \code{build_event_library()}.
#' @param seed Optional integer seed.
#'
#' @return Event library as returned by \code{build_event_library()}.
#'
#' @export
build_event_library_from_out <- function(out, location, ..., seed = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  if (is.null(out$trackpoints[[location]])) {
    stop("out$trackpoints has no entry for location='", location, "'.")
  }
  if (is.null(out$events)) {
    stop("out$events is required.", call. = FALSE)
  }
  
  track_df <- out$trackpoints[[location]]
  event_df <- out$events |>
    dplyr::filter(.data$location == location)
  
  build_event_library(
    track_df = track_df,
    event_df = event_df,
    ...,
    seed = seed
  )
}

#' Sample synthetic storm events for a single year
#'
#' @description
#' Randomly samples TS and HUR events for a calendar year using an empirical
#' event library. Each sampled event is assigned a start date (seasonality),
#' a duration, and a peak site wind speed.
#'
#' @param lib Event library produced by \code{build_event_library()}.
#' @param year Integer calendar year.
#' @param n_ts Integer number of tropical storms to sample.
#' @param n_hur Integer number of hurricanes to sample.
#' @param seed Optional integer seed.
#'
#' @return Tibble with one row per event: severity, start_date, dur_days, V_peak.
#' @export
sample_events_for_year <- function(lib, year, n_ts, n_hur, seed = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.list(lib), is.function(lib$sample_doy), is.function(lib$sample_event))
  year <- as.integer(year)
  
  get_dur_days <- function(row) {
    if ("dur_days" %in% names(row) && is.finite(row$dur_days) && row$dur_days > 0) {
      return(as.integer(row$dur_days))
    }
    if (all(c("start_time", "end_time") %in% names(row)) &&
        !is.na(row$start_time) && !is.na(row$end_time)) {
      d <- as.numeric(difftime(row$end_time, row$start_time, units = "days"))
      return(max(1L, as.integer(floor(d) + 1L)))
    }
    if ("n_points" %in% names(row) && is.finite(row$n_points) && row$n_points > 0) {
      return(max(1L, as.integer(ceiling(row$n_points / 4))))
    }
    1L
  }
  
  get_V_peak <- function(row, sev) {
    v <- NA_real_
    if ("peak_wind_kt" %in% names(row)) v <- row$peak_wind_kt
    if (!is.finite(v) && "V_site_max_kt" %in% names(row)) v <- row$V_site_max_kt
    if (!is.finite(v) && "storm_intensity_kt" %in% names(row)) v <- row$storm_intensity_kt
    if (!is.finite(v) && "wind_max_kt" %in% names(row)) v <- row$wind_max_kt
    if (!is.finite(v) || v <= 0) v <- if (sev == "HUR") 80 else if (sev == "TS") 40 else 25
    as.numeric(v)
  }
  
  sample_one <- function(sev) {
    doy0 <- as.integer(lib$sample_doy(sev))
    if (!is.finite(doy0) || doy0 < 1L || doy0 > 366L) stop("Invalid DOY sampled: ", doy0)
    
    row <- dplyr::as_tibble(lib$sample_event(sev))
    
    start_date <- as.Date(sprintf("%d-01-01", year)) + (doy0 - 1L)
    dur_days <- get_dur_days(row)
    V_peak <- get_V_peak(row, sev)
    
    tibble::tibble(
      severity   = sev,
      start_date = start_date,
      dur_days   = as.integer(dur_days),
      V_peak     = as.numeric(V_peak)
    )
  }
  
  out <- dplyr::bind_rows(
    if (n_ts  > 0) dplyr::bind_rows(replicate(n_ts,  sample_one("TS"),  simplify = FALSE)) else NULL,
    if (n_hur > 0) dplyr::bind_rows(replicate(n_hur, sample_one("HUR"), simplify = FALSE)) else NULL
  )
  
  if (nrow(out) == 0) {
    tibble::tibble(severity = character(0), start_date = as.Date(character(0)),
                   dur_days = integer(0), V_peak = numeric(0))
  } else {
    out
  }
}


#' Sample synthetic storm events for a year with extended attributes
#'
#' @description
#' Extended version of \code{sample_events_for_year()} that also carries
#' atmospheric attributes (Pc, dP, RMW) and event metadata needed by
#' \code{generate_daily_year_extended()}.
#'
#' For each sampled event, extracts from the event library row:
#' \itemize{
#'   \item \code{event_id}: unique identifier (SID or generated).
#'   \item \code{event_class}: "TC" or "HUR" (for the daily dominant-event tracker).
#'   \item \code{Pc_min_hPa}: minimum central pressure.
#'   \item \code{dP_max_hPa}: maximum pressure deficit.
#'   \item \code{RMW_mean_km}: mean radius of maximum wind.
#' }
#'
#' @param lib Event library produced by \code{build_event_library()}.
#' @param year Integer calendar year.
#' @param n_ts Integer number of tropical storms to sample.
#' @param n_hur Integer number of hurricanes to sample.
#' @param seed Optional integer seed.
#'
#' @return Tibble with columns: severity, start_date, dur_days, V_peak,
#'   event_id, event_class, Pc_min_hPa, dP_max_hPa, RMW_mean_km.
#' @export
sample_events_for_year_extended <- function(lib, year, n_ts, n_hur, seed = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")

  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.list(lib), is.function(lib$sample_doy), is.function(lib$sample_event))
  year <- as.integer(year)

  # Helper: extract duration in days from a sampled event row
  get_dur_days <- function(row) {
    if ("dur_days" %in% names(row) && is.finite(row$dur_days) && row$dur_days > 0) {
      return(as.integer(row$dur_days))
    }
    if (all(c("start_time", "end_time") %in% names(row)) &&
        !is.na(row$start_time) && !is.na(row$end_time)) {
      d <- as.numeric(difftime(row$end_time, row$start_time, units = "days"))
      return(max(1L, as.integer(floor(d) + 1L)))
    }
    if ("n_points" %in% names(row) && is.finite(row$n_points) && row$n_points > 0) {
      return(max(1L, as.integer(ceiling(row$n_points / 4))))
    }
    1L
  }

  # Helper: extract peak wind, with severity-aware fallbacks
  get_V_peak <- function(row, sev) {
    v <- NA_real_
    if ("V_site_max_kt" %in% names(row)) v <- row$V_site_max_kt
    if (!is.finite(v) && "wind_max_kt" %in% names(row)) v <- row$wind_max_kt
    if (!is.finite(v) || v <= 0) v <- if (sev == "HUR") 80 else if (sev == "TS") 40 else 25
    as.numeric(v)
  }

  # Helper: safely extract a numeric attribute, returning NA if missing
  safe_num <- function(row, col) {
    if (col %in% names(row)) {
      val <- as.numeric(row[[col]])
      if (is.finite(val)) return(val)
    }
    NA_real_
  }

  # Counter for unique event IDs within this year
  event_counter <- 0L

  sample_one <- function(sev) {
    doy0 <- as.integer(lib$sample_doy(sev))
    if (!is.finite(doy0) || doy0 < 1L || doy0 > 366L) stop("Invalid DOY sampled: ", doy0)

    row <- dplyr::as_tibble(lib$sample_event(sev))

    start_date <- as.Date(sprintf("%d-01-01", year)) + (doy0 - 1L)
    dur_days <- get_dur_days(row)
    V_peak <- get_V_peak(row, sev)

    # Event ID: use SID if available, otherwise generate
    event_counter <<- event_counter + 1L
    eid <- if ("SID" %in% names(row) && !is.na(row$SID[1])) {
      paste0(row$SID[1], "_y", year, "_", event_counter)
    } else {
      paste0("evt_", year, "_", event_counter)
    }

    # Event class for dominant-event tracking
    event_class <- if (sev == "HUR") "HUR" else "TC"

    tibble::tibble(
      severity    = sev,
      start_date  = start_date,
      dur_days    = as.integer(dur_days),
      V_peak      = as.numeric(V_peak),
      event_id    = eid,
      event_class = event_class,
      Pc_min_hPa  = dplyr::coalesce(safe_num(row, "min_pressure_hpa"), safe_num(row, "Pc_min_hPa")),
      dP_max_hPa  = dplyr::coalesce(safe_num(row, "pressure_deficit_hpa"), safe_num(row, "dP_max_hPa")),
      RMW_mean_km = dplyr::coalesce(safe_num(row, "rmw_mean_km"), safe_num(row, "RMW_mean_km"))
    )
  }

  out <- dplyr::bind_rows(
    if (n_ts  > 0) dplyr::bind_rows(replicate(n_ts,  sample_one("TS"),  simplify = FALSE)) else NULL,
    if (n_hur > 0) dplyr::bind_rows(replicate(n_hur, sample_one("HUR"), simplify = FALSE)) else NULL
  )

  if (nrow(out) == 0) {
    tibble::tibble(
      severity = character(0), start_date = as.Date(character(0)),
      dur_days = integer(0), V_peak = numeric(0),
      event_id = character(0), event_class = character(0),
      Pc_min_hPa = numeric(0), dP_max_hPa = numeric(0),
      RMW_mean_km = numeric(0)
    )
  } else {
    out
  }
}


#' Generate a parametric wind pulse for a storm event
#'
#' @description
#' Creates a deterministic within-event daily wind profile with specified duration
#' and peak wind. Used to convert discrete sampled events into a daily time series.
#'
#' @param dur_days Integer; event duration (days).
#' @param V_peak Numeric; peak wind (kt).
#' @param shape Pulse shape: "cosine" (smooth) or "triangle" (linear).
#'
#' @return Numeric vector length dur_days of daily winds (kt).
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
    w <- sin(pi * (t - 0.5) / d)
  }
  pmax(0, V_peak * w)
}




#' Generate a daily wind time series for a single calendar year
#'
#' @description
#' Converts sampled storm events into a daily time series of maximum sustained
#' wind at a site. Overlapping events are combined by daily maximum wind.
#'
#' @param year Integer calendar year.
#' @param sampled_events Tibble from \code{sample_events_for_year()}.
#' @param thr_port Optional numeric wind threshold (kt) for port disruption.
#' @param thr_infra Optional numeric wind threshold (kt) for infrastructure disruption.
#' @param pulse_shape Pulse shape passed to \code{event_pulse()}.
#'
#' @return Tibble with date, wind_kt, and optional disruption flags.
#' @export
generate_daily_year <- function(year, sampled_events,
                                thr_port = NA_real_,
                                thr_infra = NA_real_,
                                pulse_shape = "cosine") {
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  start <- as.Date(sprintf("%d-01-01", year))
  end   <- as.Date(sprintf("%d-12-31", year))
  dates <- seq.Date(start, end, by = "day")
  
  n <- length(dates)
  wind <- rep(0, n)
  
  if (nrow(sampled_events) > 0) {
    for (k in seq_len(nrow(sampled_events))) {
      s <- sampled_events$start_date[k]
      d <- sampled_events$dur_days[k]
      V <- sampled_events$V_peak[k]
      
      idx0 <- as.integer(s - start) + 1L
      idx1 <- idx0 + d - 1L
      if (idx1 < 1L || idx0 > n) next
      
      idx0c <- max(1L, idx0)
      idx1c <- min(n, idx1)
      
      pulse <- event_pulse(d, V, shape = pulse_shape)
      ps <- idx0c - idx0 + 1L
      pe <- ps + (idx1c - idx0c)
      
      wind[idx0c:idx1c] <- pmax(wind[idx0c:idx1c], pulse[ps:pe])
    }
  }
  
  tibble::tibble(
    date = dates,
    wind_kt = wind
  ) |>
    dplyr::mutate(
      port_disrupt  = if (is.finite(thr_port))  .data$wind_kt >= thr_port  else NA,
      infra_disrupt = if (is.finite(thr_infra)) .data$wind_kt >= thr_infra else NA
    )
}


#' Generate a daily wind + dominant-event attribute series for a single calendar year
#'
#' @description
#' Extended daily series used for SD coupling: tracks dominant event per day
#' and carries key storm attributes (pressure, pressure deficit, RMW).
#'
#' @param year Integer calendar year.
#' @param sampled_events Tibble from \code{sample_events_for_year_extended()}.
#' @param pulse_shape Pulse shape passed to \code{event_pulse()}.
#'
#' @return Tibble with columns: date, wind_kt, event_id, event_class, pressure_hpa,
#'   pressure_deficit_hpa, rmw_km.
#' @export
generate_daily_year_extended <- function(year, sampled_events,
                                         pulse_shape = "cosine") {
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  start <- as.Date(sprintf("%d-01-01", year))
  end   <- as.Date(sprintf("%d-12-31", year))
  dates <- seq.Date(start, end, by = "day")
  
  n <- length(dates)
  wind <- rep(0, n)
  
  # Dominant event per day
  event_id <- rep(NA_character_, n)
  event_class <- rep(NA_character_, n)
  pressure_hpa  <- rep(NA_real_, n)
  pressure_deficit_hpa  <- rep(NA_real_, n)
  rmw_km <- rep(NA_real_, n)
  
  best_wind_contrib <- rep(-Inf, n)
  
  if (nrow(sampled_events) > 0) {
    for (k in seq_len(nrow(sampled_events))) {
      s  <- sampled_events$start_date[k]
      d  <- sampled_events$dur_days[k]
      V  <- sampled_events$V_peak[k]
      id <- sampled_events$event_id[k]
      
      idx0 <- as.integer(s - start) + 1L
      idx1 <- idx0 + d - 1L
      if (idx1 < 1L || idx0 > n) next
      
      idx0c <- max(1L, idx0)
      idx1c <- min(n, idx1)
      
      pulse <- event_pulse(d, V, shape = pulse_shape)
      
      ps <- idx0c - idx0 + 1L
      pe <- ps + (idx1c - idx0c)
      
      contrib <- pulse[ps:pe]
      
      wind[idx0c:idx1c] <- pmax(wind[idx0c:idx1c], contrib)
      
      take <- contrib > best_wind_contrib[idx0c:idx1c]
      if (any(take)) {
        ii <- (idx0c:idx1c)[take]
        best_wind_contrib[ii] <- contrib[take]
        event_id[ii] <- id
        event_class[ii] <- sampled_events$event_class[k]
        pressure_hpa[ii]  <- sampled_events$Pc_min_hPa[k]
        pressure_deficit_hpa[ii]  <- sampled_events$dP_max_hPa[k]
        rmw_km[ii] <- sampled_events$RMW_mean_km[k]
      }
    }
  }
  
  tibble::tibble(
    date = dates,
    wind_kt = wind,
    event_id = event_id,
    event_class = event_class,
    pressure_hpa = pressure_hpa,
    pressure_deficit_hpa = pressure_deficit_hpa,
    rmw_km = rmw_km
  )
}




#' Generate daily synthetic hazard + impact series from hazard model output
#'
#' @description
#' For a given location and simulated years, creates daily hazard and damage
#' series with a fixed 15-column schema.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param location Character; target location name.
#' @param sim_years Integer vector of simulation-year indices to generate.
#' @param year0 Integer base calendar year corresponding to sim_year == 1.
#' @param gust_factor Numeric scalar to derive gust wind from sustained wind.
#' @param damage_method One of: "intensity" (uses add_damage_forcing) or "powerlaw".
#' @param damage_params List of parameters for the chosen method.
#' @param pulse_shape Passed to \code{event_pulse()}.
#' @param scenario Optional character scenario label carried to output.
#' @param seed Integer seed.
#'
#' @return Tibble with columns:
#'   \code{location}, \code{sim_year}, \code{scenario}, \code{date},
#'   \code{wind_kt}, \code{wind_gust_kt}, \code{surge_m},
#'   \code{event_id}, \code{event_class}, \code{pressure_hpa},
#'   \code{pressure_deficit_hpa}, \code{rmw_km},
#'   \code{damage_intensity}, \code{damage_rate}, \code{cum_damage}.
#' @export
generate_daily_hazard_impact <- function(
    out,
    location,
    sim_years = 1:1000,
    year0 = 2000,
    gust_factor = 1,
    damage_method = c("intensity", "powerlaw"),
    damage_params = list(),
    pulse_shape = "cosine",
    scenario = NA_character_,
    seed = 1) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  damage_method <- match.arg(damage_method)
  set.seed(seed)
  
  method <- if (!is.null(out$cfg) && !is.null(out$cfg$resampling_method)) out$cfg$resampling_method else NULL
  if (is.null(method)) method <- "stratified"
  copula_min_n <- if (!is.null(out$cfg) && !is.null(out$cfg$copula_min_n)) out$cfg$copula_min_n else 30L
  copula_k <- if (!is.null(out$cfg) && !is.null(out$cfg$copula_k)) out$cfg$copula_k else 1L
  copula_robust_scale <- if (!is.null(out$cfg) && !is.null(out$cfg$copula_robust_scale)) out$cfg$copula_robust_scale else TRUE
  
  lib <- build_event_library_from_out(
    out,
    location = location,
    seed = seed,
    resampling_method = method,
    copula_min_n = copula_min_n,
    copula_k = copula_k,
    copula_robust_scale = copula_robust_scale
  )
  
  if (is.null(out$sim)) stop("out$sim is NULL.", call. = FALSE)
  
  sim <- out$sim |>
    dplyr::filter(.data$location == location, .data$sim_year %in% sim_years)
  
  if (nrow(sim) == 0) stop("No sim years found for location in out$sim.", call. = FALSE)
  
  .get_sim_col <- function(sim, candidates) {
    hit <- candidates[candidates %in% names(sim)][1]
    if (is.na(hit)) {
      stop("Missing expected column in out$sim. Tried: ",
           paste(candidates, collapse = ", "), call. = FALSE)
    }
    sim[[hit]]
  }
  
  n_ts_vec  <- .get_sim_col(sim, c("n_ts"))
  n_hur_vec <- .get_sim_col(sim, c("n_hur"))
  
  # --- Level 3: resolve SST anomaly vector for storm perturbation ---
  # sst_scenario can be passed explicitly, or extracted from out$sst_info$scenario
  sst_anom_vec <- if (!is.null(out$cfg) && !is.null(out$cfg$sst_scenario)) out$cfg$sst_scenario$sst_anomaly else NULL
  cc_params <- if (!is.null(out$fit)) attr(out$fit, "cc_params") else NULL

  l3_enabled <- !is.null(cc_params) && !is.null(sst_anom_vec)

  daily_list <- vector("list", nrow(sim))
  
  for (i in seq_len(nrow(sim))) {
    yr <- year0 + (sim$sim_year[i] - 1L)
    
    sampled <- sample_events_for_year_extended(
      lib = lib,
      year = yr,
      n_ts = n_ts_vec[i],
      n_hur = n_hur_vec[i]
    )
    
    # --- Level 3: perturb storm characteristics ---
    if (l3_enabled && nrow(sampled) > 0) {
      sy <- sim$sim_year[i]
      delta_sst_i <- if (sy >= 1L && sy <= length(sst_anom_vec)) {
        sst_anom_vec[sy]
      } else {
        0
      }
      sampled <- perturb_event(sampled, delta_sst = delta_sst_i, cc_params = cc_params)
    }
    
    daily0 <- generate_daily_year_extended(
      year = yr,
      sampled_events = sampled,
      pulse_shape = pulse_shape
    ) |>
      dplyr::mutate(
        sim_year = sim$sim_year[i],
        location = location,
        scenario = scenario
      )
    
    V0 <- if (!is.null(damage_params$V0)) damage_params$V0 else 34
    V1 <- if (!is.null(damage_params$V1)) damage_params$V1 else 120
    p_exp <- if (!is.null(damage_params$p)) damage_params$p else 3
    daily0 <- daily0 |>
      dplyr::mutate(
        damage_intensity = pmax(0, pmin(1, (.data$wind_kt - V0) / (V1 - V0)))^p_exp
      )

    if (damage_method == "intensity") {
      pars <- modifyList(
        list(V0 = 34, V1 = 120, p = 3, dmax = 0.02),
        damage_params
      )
      daily1 <- do.call(add_damage_forcing, c(list(daily = daily0), pars))
      daily1 <- daily1 |>
        dplyr::mutate(cum_damage = cumsum(dplyr::coalesce(.data$damage_rate, 0)))
    } else {
      pars <- modifyList(
        list(thr = 34, V_ref = 80, d_ref = 0.03, p = 3, d_max = 0.10),
        damage_params
      )
      daily1 <- daily0 |>
        dplyr::mutate(
          damage_rate = do.call(damage_rate_from_wind, c(list(wind_kt = .data$wind_kt), pars)),
          cum_damage = cumsum(dplyr::coalesce(.data$damage_rate, 0))
        )
    }
    daily_list[[i]] <- daily1
  }
  result <- dplyr::bind_rows(daily_list) |>
    dplyr::mutate(
      wind_gust_kt = .data$wind_kt * gust_factor,
      surge_m = ifelse(is.finite(.data$pressure_hpa), 0.14 * (1013 - .data$pressure_hpa), NA_real_)
    ) |>
    dplyr::select(
      .data$location, .data$sim_year, .data$scenario, .data$date,
      .data$wind_kt, .data$wind_gust_kt, .data$surge_m,
      .data$event_id, .data$event_class, .data$pressure_hpa,
      .data$pressure_deficit_hpa, .data$rmw_km, .data$damage_intensity,
      .data$damage_rate, .data$cum_damage
    ) |>
    dplyr::relocate(.data$location, .data$sim_year, .data$scenario, .data$date)
  attr(result, "gust_factor") <- gust_factor
  result
}




# =============================================================================
# 10) Damage forcing
# =============================================================================

#' Add hazard intensity and damage forcing from daily wind
#'
#' @description
#' Maps daily wind speed to a bounded hazard intensity (0..1) and a bounded daily
#' damage rate (0..dmax).
#'
#' @param daily Tibble/data.frame with at least `wind_kt`.
#' @param V0 Numeric; threshold wind (kt) below which intensity is zero.
#' @param V1 Numeric; wind (kt) at which intensity saturates at 1.
#' @param p Numeric; nonlinearity exponent.
#' @param dmax Numeric; maximum daily damage fraction.
#'
#' @return Tibble with added columns:
#'   damage_intensity, damage_rate.
#'
#' @export
add_damage_forcing <- function(daily,
                               V0 = 34, V1 = 120,
                               p = 3,
                               dmax = 0.02) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  stopifnot(is.data.frame(daily))
  if (!("wind_kt" %in% names(daily))) stop("daily must contain `wind_kt`.", call. = FALSE)
  
  clip01 <- function(x) pmax(0, pmin(1, x))
  
  x <- (daily$wind_kt - V0) / (V1 - V0)
  I <- clip01(x)^p
  d <- dmax * I
  
  out <- dplyr::as_tibble(daily) |>
    dplyr::mutate(
      damage_intensity = I,
      damage_rate = d
    )
  out
}





#' Bounded power-law damage rate from wind speed
#'
#' @description
#' Converts wind speed to a daily damage fraction using a power law above a
#' threshold, scaled such that damage equals d_ref at V_ref and capped at d_max.
#'
#' @param wind_kt Numeric vector of wind speeds (kt).
#' @param thr Numeric threshold wind (kt) below which damage is zero.
#' @param V_ref Numeric reference wind (kt) where damage equals d_ref.
#' @param d_ref Numeric damage fraction at V_ref.
#' @param p Numeric exponent controlling nonlinearity.
#' @param d_max Numeric cap on daily damage fraction.
#'
#' @return Numeric vector of damage rates (fractions per day).
#' @export
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
