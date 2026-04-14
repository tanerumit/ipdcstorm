# =============================================================================
# Temporal downscaling and impact forcing
# - Event classification helpers.
# - Event-library construction and resampling.
# - Daily event-to-time-series conversion.
# - Daily hazard-impact generation (independent and spatially coherent).
# - Damage forcing helpers.
# =============================================================================

# =============================================================================
# 1) Event classification helpers
# =============================================================================

#' Assign simple severity class from peak wind
#'
#' @description
#' Applies fixed tropical-depression, tropical-storm, and hurricane wind
#' thresholds to a realized peak wind value.
#'
#' @param wind_max_kt Numeric; peak wind (kt).
#' @return Character scalar in \code{c("TD", "TS", "HUR")} (or \code{NA}).
#' @keywords internal
.assign_severity_simple <- function(wind_max_kt) {
  dplyr::case_when(
    !is.finite(wind_max_kt) ~ NA_character_,
    wind_max_kt >= 64 ~ "HUR",
    wind_max_kt >= 34 ~ "TS",
    TRUE ~ "TD"
  )
}

#' Classify realized downscaled event peaks
#'
#' @description
#' Classifies realized downscaled event peaks using the canonical
#' \code{classify_severity()} thresholds when available and otherwise falls back
#' to the local simple classifier.
#'
#' @param peak_kt Numeric vector of realized event peak winds (kt).
#' @return Character vector of storm classes with \code{NA} for unknown peaks.
#' @keywords internal
.classify_downscaled_event_peak <- function(peak_kt) {
  peak_kt <- as.numeric(peak_kt)

  if (any(is.finite(peak_kt) & peak_kt < 0, na.rm = TRUE)) {
    stop("`peak_kt` must not contain negative realized event peaks.", call. = FALSE)
  }

  classify_fun <- get0("classify_severity", mode = "function")
  if (is.function(classify_fun)) {
    out <- classify_fun(
      peak_kt,
      ts_threshold_kt = 34,
      hurricane_threshold_kt = 64
    )
    out[out == "unknown"] <- NA_character_
    return(out)
  }

  .assign_severity_simple(peak_kt)
}


# =============================================================================
# 2) Event library construction
# =============================================================================

#' Build an empirical event library for resampling
#'
#' @description
#' Builds a resampling library from historical storm events for one target
#' location. The library stores empirical day-of-year samples by storm class,
#' stratified historical event bins, and sampler closures used by the daily
#' downscaling workflow.
#'
#' @param track_df Tibble/data.frame of track points for one location; must
#'   contain columns \code{SID} and \code{iso_time}.
#' @param event_df Tibble/data.frame of storm events with one row per storm; must
#'   contain \code{SID} or \code{storm_id}.
#' @param storm_classes Character vector of storm classes retained in the
#'   library.
#' @param bins Named list of numeric break vectors for wind, pressure, and RMW
#'   stratification.
#' @param seed Optional integer seed for deterministic library construction.
#' @param resampling_method Character scalar; one of \code{"stratified"} or
#'   \code{"copula_nn"}.
#' @param copula_min_n Integer scalar; minimum complete-case sample size needed
#'   to fit a class-specific copula sampler.
#' @param copula_k Integer scalar; number of nearest neighbours sampled from a
#'   copula proposal.
#' @param copula_robust_scale Logical scalar; if \code{TRUE}, use median/MAD
#'   scaling in nearest-neighbour distance calculations.
#'
#' @return List with empirical day-of-year samples, stratification bins, the
#'   processed event table, and sampler functions.
#' @examples
#' track_df <- tibble::tibble(
#'   SID = c("A", "A", "B", "B"),
#'   iso_time = as.POSIXct(c(
#'     "2000-08-01 00:00:00", "2000-08-01 06:00:00",
#'     "2000-09-10 00:00:00", "2000-09-10 06:00:00"
#'   ), tz = "UTC")
#' )
#' event_df <- tibble::tibble(
#'   SID = c("A", "B"),
#'   peak_wind_kt = c(45, 80),
#'   storm_intensity_kt = c(45, 80),
#'   min_pressure_hpa = c(995, 970),
#'   rmw_mean_km = c(40, 25),
#'   start_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-09-10 00:00:00"), tz = "UTC")
#' )
#' lib <- build_event_library(track_df, event_df, storm_classes = c("TS", "HUR"))
#' names(lib)
#' @seealso \code{\link{build_event_library_from_out}}
#' @export
build_event_library <- function(track_df, event_df,
                                storm_classes = c("TD", "TS", "HUR"),
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
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package `lubridate` is required.")

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
      severity = factor(.data$severity, levels = storm_classes),
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
      feat   = feat,
      L      = L,
      ev_sub = ev_sub,
      X      = as.matrix(ev_sub[, feat, drop = FALSE]),
      Zstd   = std$Z,
      center = std$center,
      scale  = std$scale
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

    # Use pre-computed center/scale from fit rather than recomputing per draw
    x_std <- (x_star - fit$center) / fit$scale

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
    for (sev in storm_classes) {
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

# =============================================================================
# 3) Daily summary and disruption helpers
# =============================================================================

#' Compute disruption flags from daily hazard output
#'
#' @description
#' Adds logical disruption indicators to a daily hazard table by comparing daily
#' wind or surge values with user-supplied thresholds.
#'
#' @param daily Tibble/data.frame returned by
#'   \code{generate_daily_hazard_impact()}.
#' @param thr_port Numeric scalar wind threshold (kt) for port disruption.
#' @param thr_infra Numeric scalar wind threshold (kt) for infrastructure
#'   disruption.
#' @param thr_surge Numeric scalar surge threshold (m) for surge disruption.
#' @param use_gust Logical scalar; if \code{TRUE}, use \code{wind_gust_kt}
#'   instead of \code{wind_kt}.
#'
#' @return Input tibble with added logical disruption columns.
#' @examples
#' daily <- tibble::tibble(
#'   wind_kt = c(20, 45),
#'   wind_gust_kt = c(24, 54),
#'   surge_m = c(0.1, 0.8)
#' )
#' disruption_flags(daily, thr_port = 34, thr_infra = 50, thr_surge = 0.5)
#' @seealso \code{\link{generate_daily_hazard_impact}}
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
#' @description
#' Flags days in a daily hazard table that are associated with any tropical
#' storm or hurricane event.
#'
#' @param daily Tibble/data.frame returned by
#'   \code{generate_daily_hazard_impact()}.
#'
#' @return Logical vector identifying tropical-storm or hurricane days.
#' @examples
#' daily <- tibble::tibble(event_class = c(NA, "TS", "HUR"))
#' is_tc_day(daily)
#' is_hur_day(daily)
#' @seealso \code{\link{generate_daily_hazard_impact}}
#' @export
is_tc_day <- function(daily) {
  daily$event_class %in% c("TS", "HUR")
}

#' @rdname is_tc_day
#' @export
is_hur_day <- function(daily) {
  daily$event_class == "HUR"
}

#' Compute daily exposure hours above a wind threshold
#'
#' @description
#' Converts a daily wind threshold exceedance test into daily exposure hours at
#' package daily resolution.
#'
#' @param daily Tibble/data.frame returned by
#'   \code{generate_daily_hazard_impact()}.
#' @param threshold_kt Numeric scalar wind threshold (kt).
#' @param use_gust Logical scalar; if \code{TRUE}, use \code{wind_gust_kt}.
#'
#' @return Numeric vector of daily exposure hours at package daily resolution.
#' @examples
#' daily <- tibble::tibble(wind_kt = c(20, 40), wind_gust_kt = c(24, 48))
#' exposure_hours(daily, threshold_kt = 34)
#' @export
exposure_hours <- function(daily, threshold_kt, use_gust = FALSE) {
  wind_col <- if (use_gust) "wind_gust_kt" else "wind_kt"
  24 * (daily[[wind_col]] >= threshold_kt)
}

#' Summarise peak wind per simulation year
#'
#' @description
#' Aggregates a daily hazard table to annual peak sustained wind by location,
#' simulation year, and scenario.
#'
#' @param daily Tibble/data.frame returned by
#'   \code{generate_daily_hazard_impact()}.
#'
#' @return Tibble with one row per location and simulation year and a
#'   \code{peak_wind_kt} column.
#' @examples
#' daily <- tibble::tibble(
#'   location = c("Saba", "Saba"),
#'   sim_year = c(1, 1),
#'   scenario = c("baseline", "baseline"),
#'   wind_kt = c(20, 55)
#' )
#' peak_wind_by_year(daily)
#' @export
peak_wind_by_year <- function(daily) {
  daily |>
    dplyr::group_by(.data$location, .data$sim_year, .data$scenario) |>
    dplyr::summarise(peak_wind_kt = max(.data$wind_kt, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(peak_wind_kt = dplyr::if_else(is.finite(.data$peak_wind_kt), .data$peak_wind_kt, NA_real_))
}


#' Convenience wrapper: build an event library from run_hazard_model() output
#'
#' @description
#' Extracts one location's track-point and event tables from
#' \code{run_hazard_model()} output and forwards them to
#' \code{build_event_library()}.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param location Character scalar target location name.
#' @param ... Additional arguments passed to \code{build_event_library()}.
#' @param seed Optional integer seed for deterministic library construction.
#'
#' @return List in the format returned by \code{build_event_library()}.
#' @examples
#' \dontrun{
#' lib <- build_event_library_from_out(out, location = "Saba")
#' lib$events
#' }
#' @seealso \code{\link{build_event_library}}
#' @export
build_event_library_from_out <- function(out, location, ..., seed = NULL) {

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


# =============================================================================
# 4) Event sampling and daily pulse generation
# =============================================================================

#' Sample synthetic storm events for a year with extended attributes
#'
#' @description
#' Samples synthetic storm events for one calendar year from an event library
#' and carries forward the atmospheric attributes and identifiers needed by
#' \code{generate_daily_year_extended()}.
#'
#' For each sampled event, extracts from the event library row:
#' \itemize{
#'   \item \code{event_id}: unique identifier (SID or generated).
#'   \item \code{event_class}: "TS" or "HUR" (for the daily dominant-event tracker).
#'   \item \code{Pc_min_hPa}: minimum central pressure.
#'   \item \code{dP_max_hPa}: maximum pressure deficit.
#'   \item \code{RMW_mean_km}: mean radius of maximum wind.
#' }
#'
#' @param lib List event library produced by \code{build_event_library()}.
#' @param year Integer scalar calendar year.
#' @param n_ts Integer scalar number of tropical storms to sample.
#' @param n_hur Integer scalar number of hurricanes to sample.
#' @param seed Optional integer seed for deterministic sampling.
#'
#' @return Tibble with one row per sampled event and event-level attributes used
#'   by the daily downscaling workflow.
#' @examples
#' sample_lib <- list(
#'   sample_doy = function(sev) if (sev == "TS") 220L else 250L,
#'   sample_event = function(sev) {
#'     tibble::tibble(
#'       SID = paste0(sev, "_1"),
#'       dur_days = if (sev == "TS") 3L else 4L,
#'       V_site_max_kt = if (sev == "TS") 45 else 80,
#'       Pc_min_hPa = if (sev == "TS") 995 else 970,
#'       dP_max_hPa = if (sev == "TS") 18 else 40,
#'       RMW_mean_km = if (sev == "TS") 40 else 25
#'     )
#'   }
#' )
#' sample_events_for_year_extended(sample_lib, year = 2000, n_ts = 1, n_hur = 1)
#' @seealso \code{\link{generate_daily_year_extended}}
#' @export
sample_events_for_year_extended <- function(lib, year, n_ts, n_hur, seed = NULL) {

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
    event_class <- if (sev == "HUR") "HUR" else "TS"

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
#' Creates a deterministic within-event daily wind profile from an event
#' duration and peak sustained wind.
#'
#' @param dur_days Integer scalar event duration in days.
#' @param V_peak Numeric scalar peak sustained wind in kt.
#' @param shape Character scalar pulse shape; one of \code{"cosine"} or
#'   \code{"triangle"}.
#'
#' @return Numeric vector of daily sustained wind values in kt.
#' @examples
#' event_pulse(dur_days = 5, V_peak = 60)
#' @seealso \code{\link{generate_daily_year_extended}}
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


#' Generate a daily wind + dominant-event attribute series for a single calendar year
#'
#' @description
#' Converts sampled storm events for one calendar year into a daily sustained
#' wind series while tracking the dominant event and key event attributes for
#' each day.
#'
#' @param year Integer scalar calendar year.
#' @param sampled_events Tibble from \code{sample_events_for_year_extended()}.
#' @param pulse_shape Character scalar pulse shape passed to
#'   \code{event_pulse()}.
#'
#' @return Tibble with one row per calendar day and daily hazard attributes.
#' @examples
#' sampled_events <- tibble::tibble(
#'   start_date = as.Date("2000-08-01"),
#'   dur_days = 3L,
#'   V_peak = 60,
#'   event_id = "evt_1",
#'   event_class = "TS",
#'   Pc_min_hPa = 995,
#'   dP_max_hPa = 18,
#'   RMW_mean_km = 40
#' )
#' generate_daily_year_extended(2000, sampled_events)
#' @seealso \code{\link{sample_events_for_year_extended}}, \code{\link{event_pulse}}
#' @export
generate_daily_year_extended <- function(year, sampled_events,
                                         pulse_shape = "cosine") {

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
  event_peak_map <- numeric(0)
  event_class_map <- character(0)

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
      realized_peak <- if (length(contrib) > 0L) max(contrib) else NA_real_
      event_class_k <- .classify_downscaled_event_peak(realized_peak)

      if (!is.na(id) && is.finite(realized_peak)) {
        event_peak_map[[id]] <- realized_peak
        event_class_map[[id]] <- event_class_k
      }

      wind[idx0c:idx1c] <- pmax(wind[idx0c:idx1c], contrib)

      take <- contrib > best_wind_contrib[idx0c:idx1c]
      if (any(take)) {
        ii <- (idx0c:idx1c)[take]
        best_wind_contrib[ii] <- contrib[take]
        event_id[ii] <- id
        event_class[ii] <- event_class_k
        pressure_hpa[ii]  <- sampled_events$Pc_min_hPa[k]
        pressure_deficit_hpa[ii]  <- sampled_events$dP_max_hPa[k]
        rmw_km[ii] <- sampled_events$RMW_mean_km[k]
      }
    }
  }

  if (length(event_peak_map) > 0) {
    event_check <- tibble::tibble(
      event_id = names(event_peak_map),
      realized_peak_kt = as.numeric(event_peak_map),
      event_class = unname(event_class_map)
    ) |>
      dplyr::mutate(
        implied_class = .classify_downscaled_event_peak(.data$realized_peak_kt),
        mismatch = !is.na(.data$event_class) &
          !is.na(.data$implied_class) &
          .data$event_class != .data$implied_class
      )

    if (any(event_check$mismatch, na.rm = TRUE)) {
      stop("Internal error: inconsistent `event_class` after realized-peak classification.", call. = FALSE)
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


# =============================================================================
# 5) Daily hazard-impact generation
# =============================================================================

#' Generate daily synthetic hazard and impact time series from hazard-model output
#'
#' @description
#' Converts stochastic annual hazard-model output from \code{run_hazard_model()}
#' into daily synthetic hazard and impact series for one or more target
#' locations and selected simulation years.
#'
#' For each requested location, the function:
#' \enumerate{
#'   \item builds an empirical event library from the location-specific
#'   historical track-point and event information stored in \code{out};
#'   \item reads the simulated annual storm counts from \code{out$sim}
#'   (currently tropical storms and hurricanes via \code{n_ts} and
#'   \code{n_hur});
#'   \item resamples synthetic events year by year from the empirical library;
#'   \item optionally applies event-level climate perturbations when such
#'   metadata are attached to \code{out$fit};
#'   \item downscales each sampled event to a daily sustained-wind pulse;
#'   \item combines overlapping events into a daily series, retaining the
#'   dominant event identifier and associated pressure / RMW attributes for each
#'   day; and
#'   \item derives gust wind, a simple pressure-based surge proxy, daily damage
#'   forcing, and cumulative damage.
#' }
#'
#' The function is intended as a temporal downscaling bridge between the annual
#' synthetic hazard model and downstream impact models that operate at daily
#' time step. It does not re-simulate storm occurrence from first principles;
#' instead, it reuses the synthetic annual counts already generated by
#' \code{run_hazard_model()} and combines them with empirical event structure
#' extracted from the historical event library. The return value is always a
#' named list of tibbles, one element per requested location.
#'
#' @details
#' The daily series is generated through the following helper workflow:
#'
#' \enumerate{
#'   \item \code{\link{build_event_library_from_out}()} extracts the
#'   location-specific historical event library from \code{out$trackpoints} and
#'   \code{out$events}. Depending on configuration stored in \code{out$cfg}, the
#'   underlying sampler can use either stratified historical resampling or the
#'   \code{"copula_nn"} nearest-neighbour copula workflow.
#'   \item \code{sample_events_for_year_extended()} samples the required number
#'   of tropical storms and hurricanes for each simulation year using the
#'   annual counts in \code{out$sim}. Sampled events carry start date, duration,
#'   peak wind, central pressure, pressure deficit, and mean radius of maximum
#'   wind.
#'   \item If climate perturbation metadata are attached to \code{out$fit}, the
#'   sampled event set is modified with \code{perturb_event()} before daily
#'   downscaling. This is only applied when perturbation settings are present
#'   and a finite scalar \code{delta_sst} is available.
#'   \item \code{\link{generate_daily_year_extended}()} converts each sampled
#'   event into a daily time series by applying \code{\link{event_pulse}()} to
#'   the event duration and peak wind, then merges overlapping events by taking
#'   the daily maximum wind contribution. The dominant event for each day
#'   carries the \code{event_id}, \code{event_class}, \code{pressure_hpa},
#'   \code{pressure_deficit_hpa}, and \code{rmw_km} fields into the final daily
#'   table.
#'   \item Damage forcing is then added using either
#'   \code{\link{add_damage_forcing}()} or
#'   \code{\link{damage_rate_from_wind}()}, depending on
#'   \code{damage$method}. \code{cum_damage} is computed as the cumulative sum
#'   of daily \code{damage_rate} within the generated series.
#'   \item Finally, \code{wind_gust_kt} is derived as
#'   \code{wind_kt * gust_factor}, and \code{surge_m} is derived from central
#'   pressure using the simple pressure-deficit proxy
#'   \code{0.14 * (1013 - pressure_hpa)} where pressure is available.
#' }
#'
#' Scientific parameterization:
#'
#' \itemize{
#'   \item \strong{Gust factor.}
#'   \code{gust_factor} is a simple multiplicative conversion from daily
#'   sustained wind speed to daily gust wind speed:
#'   \deqn{wind\_gust\_kt = wind\_kt \times gust\_factor}
#'   No additional boundary-layer or exposure correction is applied inside this
#'   function. A value of 1 leaves gust equal to sustained wind; values above 1
#'   represent stronger short-duration peak gusts.
#'
#'   \item \strong{Pulse shape.}
#'   \code{pulse_shape} controls the within-event temporal evolution of daily
#'   sustained wind through \code{\link{event_pulse}()}. Supported values are
#'   currently \code{"cosine"} and \code{"triangle"}. The cosine option uses
#'   \code{sin(pi * (t - 0.5) / d)} over event day \code{t = 1, ..., d}, giving
#'   a smooth rise and fall with zero at the event edges. The triangle option
#'   uses a piecewise-linear ramp up and down around the event midpoint. In both
#'   cases the profile is scaled to the sampled event peak wind \code{V_peak}
#'   and truncated at zero.
#'
#'   \item \strong{Damage intensity index.}
#'   Regardless of \code{damage$method}, the function first computes a bounded
#'   hazard-intensity index from sustained daily wind:
#'   \deqn{I = \mathrm{clip}_{[0,1]}\left(\frac{wind\_kt - V0}{V1 - V0}\right)^p}
#'   where \code{V0} is the lower damage threshold, \code{V1} is the wind speed
#'   at which the normalized intensity saturates at 1, and \code{p} controls
#'   convexity. This index is stored in \code{damage_intensity}. For
#'   \code{damage$method = "intensity"}, these values come from the intensity
#'   damage specification; for \code{damage$method = "powerlaw"}, the generic
#'   index uses \code{V0 = 34} kt and \code{V1 = 120} kt with the shared
#'   exponent \code{p}.
#'
#'   \item \strong{Damage method = "intensity".}
#'   This route calls \code{\link{add_damage_forcing}()}, which maps the bounded
#'   intensity index to daily damage rate as
#'   \deqn{damage\_rate = dmax \times I}
#'   where \code{dmax} is the maximum daily damage fraction. Defaults are
#'   \code{V0 = 34} kt, \code{V1 = 120} kt, \code{p = 3}, and
#'   \code{dmax = 0.02}. This is a capped intensity-index model rather than a
#'   direct engineering vulnerability curve.
#'
#'   \item \strong{Damage method = "powerlaw".}
#'   This route calls \code{\link{damage_rate_from_wind}()}, which applies a
#'   thresholded power law calibrated at a reference wind:
#'   \deqn{x = \max\left(0, \frac{wind\_kt - thr}{V\_ref - thr}\right)}
#'   \deqn{damage\_rate = \min(d\_max,\; d\_ref \times x^p)}
#'   Here \code{thr} is the wind threshold below which damage is zero,
#'   \code{V_ref} is a calibration wind at which damage equals \code{d_ref},
#'   \code{p} is the exponent, and \code{d_max} is the daily cap. Defaults are
#'   \code{thr = 34} kt, \code{V_ref = 80} kt, \code{d_ref = 0.03},
#'   \code{p = 3}, and \code{d_max = 0.10}. This method exposes the reference
#'   point explicitly and is easier to align with an external damage assumption
#'   at a chosen hazard level.
#'
#'   \item \strong{Cumulative damage.}
#'   \code{cum_damage} is computed as the running cumulative sum of
#'   \code{damage_rate} over the generated daily sequence for each location.
#'   It is therefore a synthetic cumulative forcing index, not necessarily a
#'   literal stock-loss fraction unless the user calibrates the daily damage
#'   formulation accordingly.
#' }
#'
#' @param out List returned by \code{run_hazard_model()}. Must contain
#'   location-specific simulation output in \code{out$sim} and the historical
#'   event / track-point information needed to build the empirical event library.
#'
#' @param location Character vector of one or more target location names for
#'   which daily series should be generated. Each requested location must be
#'   present in the hazard-model output.
#'
#' @param sim_years Integer vector of simulation-year indices to generate.
#'   These refer to the synthetic years in \code{out$sim}, not historical
#'   calendar years.
#'
#' @param year0 Integer scalar base calendar year used to map synthetic year
#'   indices to dates. The function interprets \code{sim_year == 1} as
#'   \code{year0}, \code{sim_year == 2} as \code{year0 + 1}, and so on.
#'
#' @param gust_factor Numeric scalar multiplier used to convert daily sustained
#'   wind \code{wind_kt} to daily gust wind \code{wind_gust_kt}. No internal
#'   validation against observed gust-duration conventions is performed.
#'
#' @param damage Named list defining the daily damage model. Must include a
#'   scalar character element \code{method} equal to \code{"intensity"} or
#'   \code{"powerlaw"}. Recognized method-specific fields are:
#'   \describe{
#'     \item{\code{method = "intensity"}}{\code{V0}, \code{V1}, \code{p},
#'       \code{dmax}. Defaults are \code{V0 = 34} kt, \code{V1 = 120} kt,
#'       \code{p = 3}, and \code{dmax = 0.02}.}
#'     \item{\code{method = "powerlaw"}}{\code{thr}, \code{V_ref},
#'       \code{d_ref}, \code{p}, \code{d_max}. Defaults are \code{thr = 34} kt,
#'       \code{V_ref = 80} kt, \code{d_ref = 0.03}, \code{p = 3}, and
#'       \code{d_max = 0.10}.}
#'   }
#'   Unknown fields for the selected method are rejected. The output still
#'   includes \code{damage_intensity}; for \code{method = "powerlaw"}, it is a
#'   generic bounded wind-intensity index and does not change the power-law
#'   damage-rate formula.
#'
#' @param pulse_shape Character scalar pulse shape passed to
#'   \code{\link{event_pulse}()}. Currently supported values are
#'   \code{"cosine"} and \code{"triangle"}.
#'
#' @param scenario Optional character scalar scenario label copied into the
#'   output rows without modification.
#'
#' @param seed Optional integer scalar random seed. When \code{NULL} (the
#'   default), the seed is inherited from \code{out$run_metadata$seed} if
#'   available, and falls back to \code{1L} otherwise. Pass an explicit integer
#'   to override and decouple from the hazard-model run. The wrapper offsets the
#'   resolved seed by location index so repeated calls are deterministic but
#'   distinct across requested locations.
#'
#' @return
#' A named list of tibbles, one element per requested location. Each tibble
#' contains a daily synthetic hazard-impact series with columns:
#' \describe{
#'   \item{\code{location}}{Target location name.}
#'   \item{\code{sim_year}}{Synthetic simulation-year index from \code{out$sim}.}
#'   \item{\code{scenario}}{Scenario label passed through from the function call.}
#'   \item{\code{date}}{Daily calendar date derived from \code{year0}.}
#'   \item{\code{wind_kt}}{Daily sustained wind speed (kt).}
#'   \item{\code{wind_gust_kt}}{Daily gust wind speed (kt), computed using
#'     \code{gust_factor}.}
#'   \item{\code{surge_m}}{Simple pressure-based surge proxy (m).}
#'   \item{\code{event_id}}{Identifier of the dominant event contributing the
#'     highest wind on that day, or \code{NA} outside events.}
#'   \item{\code{event_class}}{Dominant daily event class, currently \code{"TS"}
#'     or \code{"HUR"}, based on realized downscaled peak.}
#'   \item{\code{pressure_hpa}}{Minimum central pressure associated with the
#'     dominant event on that day, when available.}
#'   \item{\code{pressure_deficit_hpa}}{Pressure deficit associated with the
#'     dominant event on that day, when available.}
#'   \item{\code{rmw_km}}{Mean radius of maximum wind for the dominant event,
#'     when available.}
#'   \item{\code{damage_intensity}}{Bounded normalized hazard-intensity index.}
#'   \item{\code{damage_rate}}{Daily damage forcing fraction from the selected
#'     damage model.}
#'   \item{\code{cum_damage}}{Cumulative sum of \code{damage_rate}.}
#' }
#' Each location tibble also carries a \code{"gust_factor"} attribute.
#'
#' @seealso
#' \code{\link{build_event_library_from_out}},
#' \code{\link{generate_daily_year_extended}},
#' \code{\link{event_pulse}},
#' \code{\link{add_damage_forcing}},
#' \code{\link{damage_rate_from_wind}}
#'
#' @examples
#' \dontrun{
#' # Single location — returns a named list of length 1
#' res <- generate_daily_hazard_impact(
#'   out = out,
#'   location = "Saba",
#'   sim_years = 1:10,
#'   year0 = 2000,
#'   gust_factor = 1.3,
#'   damage = list(method = "intensity", V0 = 34, V1 = 120, p = 3, dmax = 0.02),
#'   pulse_shape = "cosine",
#'   scenario = "baseline",
#'   seed = 1
#' )
#' res$Saba
#'
#' # Multiple locations with the power-law damage option
#' res <- generate_daily_hazard_impact(
#'   out = out,
#'   location = c("Saba", "St. Eustatius"),
#'   sim_years = 1:5,
#'   damage = list(
#'     method = "powerlaw",
#'     thr = 34,
#'     V_ref = 80,
#'     d_ref = 0.03,
#'     p = 3,
#'     d_max = 0.10
#'   )
#' )
#' res$Saba
#' res$`St. Eustatius`
#' }
#' @export
generate_daily_hazard_impact <- function(
    out,
    location,
    sim_years = 1:1000,
    year0 = 2000,
    gust_factor = 1,
    damage = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario = NA_character_,
    seed = NULL) {
  stopifnot(is.character(location), length(location) >= 1L)

  if (is.null(seed)) {
    seed <- if (!is.null(out$run_metadata$seed)) out$run_metadata$seed else 1L
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be NULL or a single finite numeric value.", call. = FALSE)
  }
  seed <- as.integer(seed)

  damage <- .validate_damage_spec(damage)

  # --- Run each location via the internal worker ---
  results <- stats::setNames(
    vector("list", length(location)),
    location
  )

  for (loc in location) {
    # Per-location seed: deterministic but distinct across locations
    loc_seed <- seed + match(loc, location) - 1L

    results[[loc]] <- .generate_daily_hazard_impact_single(
      out = out,
      location = loc,
      sim_years = sim_years,
      year0 = year0,
      gust_factor = gust_factor,
      damage = damage,
      pulse_shape = pulse_shape,
      scenario = scenario,
      seed = loc_seed
    )
  }

  results
}

#' Validate and fill damage specification defaults
#'
#' @param damage Named list supplied to \code{generate_daily_hazard_impact()}.
#' @return Named list with validated method-specific defaults filled in.
#' @keywords internal
.validate_damage_spec <- function(damage) {
  defaults <- list(
    intensity = list(V0 = 34, V1 = 120, p = 3, dmax = 0.02),
    powerlaw = list(thr = 34, V_ref = 80, d_ref = 0.03, p = 3, d_max = 0.10)
  )

  if (!is.list(damage) || is.data.frame(damage)) {
    stop("`damage` must be a named list.", call. = FALSE)
  }
  if (length(damage) == 0L) {
    stop("`damage$method` must be a single non-empty character string.", call. = FALSE)
  }
  if (is.null(names(damage)) || any(names(damage) == "") || anyDuplicated(names(damage))) {
    stop("`damage` must be a named list with unique names.", call. = FALSE)
  }

  method <- damage[["method"]]
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("`damage$method` must be a single non-empty character string.", call. = FALSE)
  }
  if (!(method %in% names(defaults))) {
    stop("`damage$method` must be one of: intensity, powerlaw.", call. = FALSE)
  }

  allowed <- c("method", names(defaults[[method]]))
  unknown <- setdiff(names(damage), allowed)
  if (length(unknown) > 0L) {
    stop(
      "`damage` contains unsupported field(s) for method '", method, "': ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }

  for (nm in setdiff(names(damage), "method")) {
    value <- damage[[nm]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
      stop("`damage$", nm, "` must be a single finite numeric value.", call. = FALSE)
    }
  }

  utils::modifyList(list(method = method), defaults[[method]]) |>
    utils::modifyList(damage)
}

#' Generate daily hazard and damage for one location
#'
#' @description
#' Internal worker that builds one location's event library, samples annual
#' events from the stochastic hazard output, and assembles the final daily
#' hazard-impact table.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param location Character scalar location name.
#' @param sim_years Integer vector of simulation-year indices to generate.
#' @param year0 Integer scalar base calendar year corresponding to
#'   \code{sim_year == 1}.
#' @param gust_factor Numeric scalar gust multiplier.
#' @param damage Named list of validated damage parameters.
#' @param pulse_shape Character scalar pulse-shape identifier.
#' @param scenario Optional character scalar scenario label carried into output.
#' @param seed Integer scalar seed for deterministic per-location generation.
#' @return Tibble with the daily hazard-impact schema used by
#'   \code{generate_daily_hazard_impact()}.
#' @keywords internal
.generate_daily_hazard_impact_single <- function(
    out,
    location,
    sim_years,
    year0,
    gust_factor,
    damage,
    pulse_shape,
    scenario,
    seed) {

  # seed is forwarded to build_event_library_from_out(), which calls
  # set.seed(seed) internally; sampling after library construction inherits
  # that RNG state deterministically.
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
    dplyr::filter(.data$location == !!location, .data$sim_year %in% sim_years)

  if (nrow(sim) == 0) stop("No sim years found for location '", location, "' in out$sim.", call. = FALSE)

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

  delta_sst <- if (!is.null(out$fit)) attr(out$fit, "delta_sst") else NULL
  perturb_cfg <- if (!is.null(out$fit)) attr(out$fit, "perturb") else NULL

  perturb_enabled <- !is.null(perturb_cfg)

  # Pre-compute intensity parameters once — these are constant across all years.
  intensity_pars <- if (identical(damage$method, "intensity")) {
    damage[c("V0", "V1", "p")]
  } else {
    list(V0 = 34, V1 = 120, p = damage$p)
  }

  daily_list <- vector("list", nrow(sim))

  for (i in seq_len(nrow(sim))) {
    yr <- year0 + (sim$sim_year[i] - 1L)

    sampled <- sample_events_for_year_extended(
      lib = lib,
      year = yr,
      n_ts = n_ts_vec[i],
      n_hur = n_hur_vec[i]
    )

    # Apply storm perturbation only when explicitly enabled in the climate config.
    if (perturb_enabled && nrow(sampled) > 0) {
      if (!is.numeric(delta_sst) || length(delta_sst) != 1L || !is.finite(delta_sst)) {
        stop("delta_sst must be a single finite numeric value for perturb_event().", call. = FALSE)
      }
      sampled <- perturb_event(sampled, delta_sst = delta_sst, cc_params = perturb_cfg)
    }

    daily_list[[i]] <- generate_daily_year_extended(
      year = yr,
      sampled_events = sampled,
      pulse_shape = pulse_shape
    )
  }

  # Assemble full series, then add constant and derived columns in bulk rather
  # than repeating the same dplyr calls once per simulation year.
  result <- dplyr::bind_rows(daily_list)

  # sim_year: each year block may be 365 or 366 rows (leap years), so use
  # per-element lengths rather than assuming a fixed block size.
  n_days_per_year <- vapply(daily_list, nrow, integer(1L))
  result$sim_year  <- rep(sim$sim_year, times = n_days_per_year)
  result$location  <- location
  result$scenario  <- scenario

  # Damage intensity (vectorized over the full assembled series)
  result$damage_intensity <- pmax(
    0,
    pmin(1, (result$wind_kt - intensity_pars$V0) / (intensity_pars$V1 - intensity_pars$V0))
  )^intensity_pars$p

  # Damage rate (vectorized over the full assembled series)
  if (damage$method == "intensity") {
    result <- do.call(add_damage_forcing, c(list(daily = result), damage[c("V0", "V1", "p", "dmax")]))
  } else {
    result$damage_rate <- do.call(
      damage_rate_from_wind,
      c(list(wind_kt = result$wind_kt), damage[c("thr", "V_ref", "d_ref", "p", "d_max")])
    )
  }

  # cum_damage: per-sim_year cumulative sum using base ave() for grouping
  result$cum_damage <- ave(
    dplyr::coalesce(result$damage_rate, 0),
    result$sim_year,
    FUN = cumsum
  )

  result <- result |>
    dplyr::mutate(
      wind_gust_kt = .data$wind_kt * gust_factor,
      surge_m = ifelse(is.finite(.data$pressure_hpa), 0.14 * (1013 - .data$pressure_hpa), NA_real_)
    ) |>
    dplyr::select(
      "location", "sim_year", "scenario", "date",
      "wind_kt", "wind_gust_kt", "surge_m",
      "event_id", "event_class", "pressure_hpa",
      "pressure_deficit_hpa", "rmw_km", "damage_intensity",
      "damage_rate", "cum_damage"
    ) |>
    dplyr::relocate("location", "sim_year", "scenario", "date")
  attr(result, "gust_factor") <- gust_factor
  result
}


# =============================================================================
# 5b) Spatially coherent daily hazard-impact generation (shared event pool)
# =============================================================================

#' Build a shared event pool from multiple location event libraries
#'
#' @description
#' Unions all events from a named list of location event libraries, keeping one
#' row per SID (the row with the highest site-level peak wind across all
#' locations). Only \code{"TS"} and \code{"HUR"} class events are retained.
#'
#' @param libs Named list of event libraries as returned by
#'   \code{build_event_library_from_out()}, one element per location.
#' @return Tibble with one row per unique SID and a \code{storm_class} column
#'   alongside any shared event attributes carried from the library tables.
#' @keywords internal
.build_shared_event_pool <- function(libs) {

  .get_sid_col <- function(ev) {
    if ("SID" %in% names(ev)) return(ev$SID)
    if ("storm_id" %in% names(ev)) return(ev$storm_id)
    rep(NA_character_, nrow(ev))
  }

  .get_vpk <- function(ev) {
    if ("V_site_max_kt" %in% names(ev)) return(as.numeric(ev$V_site_max_kt))
    if ("peak_wind_kt"  %in% names(ev)) return(as.numeric(ev$peak_wind_kt))
    if ("wind_max_kt"   %in% names(ev)) return(as.numeric(ev$wind_max_kt))
    rep(NA_real_, nrow(ev))
  }

  parts <- lapply(libs, function(lib) {
    ev <- lib$events
    if (is.null(ev) || nrow(ev) == 0L) return(NULL)
    ev <- dplyr::as_tibble(ev)
    ev$pool_SID <- .get_sid_col(ev)
    ev$vpk      <- .get_vpk(ev)
    ev <- ev[!is.na(ev$pool_SID), ]
    # Drop the original SID/storm_id columns now so that rename(SID = pool_SID)
    # does not fail after bind_rows when those columns are already present.
    ev[, setdiff(names(ev), c("SID", "storm_id")), drop = FALSE]
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L)
    return(tibble::tibble(SID = character(), storm_class = character()))

  all_ev <- dplyr::bind_rows(parts)

  # Derive storm_class if absent
  if (!"storm_class" %in% names(all_ev)) {
    all_ev$storm_class <- dplyr::case_when(
      is.finite(all_ev$vpk) & all_ev$vpk >= 64 ~ "HUR",
      is.finite(all_ev$vpk) & all_ev$vpk >= 34 ~ "TS",
      TRUE ~ NA_character_
    )
  }
  all_ev <- all_ev[all_ev$storm_class %in% c("TS", "HUR"), ]
  if (nrow(all_ev) == 0L)
    return(tibble::tibble(SID = character(), storm_class = character()))

  # One row per SID — highest site-level peak wind
  all_ev |>
    dplyr::group_by(.data$pool_SID, .data$storm_class) |>
    dplyr::slice_max(order_by = .data$vpk, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(SID = "pool_SID")
}


#' Sample SIDs from the shared event pool
#'
#' @param pool Tibble from \code{.build_shared_event_pool()}.
#' @param storm_class Character scalar; one of \code{"TS"} or \code{"HUR"}.
#' @param n Integer scalar number of SIDs to draw (with replacement).
#' @return Character vector of length \code{n}, or \code{character(0)} when
#'   \code{n == 0} or the class sub-pool is empty.
#' @keywords internal
.sample_shared_sids <- function(pool, storm_class, n) {
  if (n == 0L || nrow(pool) == 0L) return(character(0))
  sub <- pool[pool$storm_class == storm_class, , drop = FALSE]
  if (nrow(sub) == 0L) return(character(0))
  sample(sub$SID, size = n, replace = TRUE)
}


#' Zero-row sampled-events tibble
#'
#' @description
#' Returns a zero-row tibble with the schema expected by
#' \code{generate_daily_year_extended()}.
#'
#' @return Zero-row tibble.
#' @keywords internal
.empty_sampled_events_tbl <- function() {
  tibble::tibble(
    severity    = character(0),
    start_date  = as.Date(character(0)),
    dur_days    = integer(0),
    V_peak      = numeric(0),
    event_id    = character(0),
    event_class = character(0),
    Pc_min_hPa  = numeric(0),
    dP_max_hPa  = numeric(0),
    RMW_mean_km = numeric(0)
  )
}


#' Build event rows from shared SIDs using a location-specific library
#'
#' @description
#' For each SID in \code{sids} that exists in \code{lib$events}, constructs one
#' event row in the schema produced by \code{sample_events_for_year_extended()}.
#' SIDs absent from the location library are silently skipped, so the number of
#' rows returned may be less than \code{length(sids)}.
#'
#' @param lib Event library from \code{build_event_library_from_out()}.
#' @param sids Character vector of SIDs to look up.
#' @param sev Character scalar severity; one of \code{"TS"} or \code{"HUR"}.
#' @param year Integer scalar calendar year for date construction and
#'   \code{event_id} encoding.
#' @param counter_start Integer scalar starting value for the within-year event
#'   counter (used to produce unique \code{event_id} values across TS and HUR
#'   draws for the same year and location).
#' @return Named list with elements:
#'   \describe{
#'     \item{\code{rows}}{Tibble of event rows (zero or more).}
#'     \item{\code{n_filled}}{Integer count of rows produced.}
#'     \item{\code{counter}}{Updated event counter.}
#'   }
#' @keywords internal
.build_event_rows_from_sids <- function(lib, sids, sev, year, counter_start) {
  year    <- as.integer(year)
  counter <- as.integer(counter_start)

  if (length(sids) == 0L || is.null(lib$events) || nrow(lib$events) == 0L)
    return(list(rows = .empty_sampled_events_tbl(), n_filled = 0L, counter = counter))

  ev <- dplyr::as_tibble(lib$events)
  if ("SID" %in% names(ev)) {
    ev$.sid <- ev$SID
  } else if ("storm_id" %in% names(ev)) {
    ev$.sid <- ev$storm_id
  } else {
    return(list(rows = .empty_sampled_events_tbl(), n_filled = 0L, counter = counter))
  }

  # Local helpers — mirrors internals of sample_events_for_year_extended()
  get_dur_days <- function(row) {
    if ("dur_days" %in% names(row) && is.finite(row$dur_days[1]) && row$dur_days[1] > 0)
      return(as.integer(row$dur_days[1]))
    if (all(c("start_time", "end_time") %in% names(row)) &&
        !is.na(row$start_time[1]) && !is.na(row$end_time[1])) {
      d <- as.numeric(difftime(row$end_time[1], row$start_time[1], units = "days"))
      return(max(1L, as.integer(floor(d) + 1L)))
    }
    if ("n_points" %in% names(row) && is.finite(row$n_points[1]) && row$n_points[1] > 0)
      return(max(1L, as.integer(ceiling(row$n_points[1] / 4))))
    1L
  }

  get_V_peak <- function(row, sev) {
    v <- NA_real_
    if ("V_site_max_kt" %in% names(row)) v <- as.numeric(row$V_site_max_kt[1])
    if (!is.finite(v) && "wind_max_kt" %in% names(row))  v <- as.numeric(row$wind_max_kt[1])
    if (!is.finite(v) && "peak_wind_kt" %in% names(row)) v <- as.numeric(row$peak_wind_kt[1])
    if (!is.finite(v) || v <= 0) v <- if (sev == "HUR") 80 else 40
    v
  }

  safe_num <- function(row, col) {
    if (col %in% names(row)) {
      val <- as.numeric(row[[col]][1])
      if (is.finite(val)) return(val)
    }
    NA_real_
  }

  get_doy <- function(row) {
    if ("doy" %in% names(row) && is.finite(row$doy[1]) &&
        row$doy[1] >= 1 && row$doy[1] <= 366)
      return(as.integer(row$doy[1]))
    if ("start_time" %in% names(row) && !is.na(row$start_time[1]))
      return(as.integer(format(as.Date(row$start_time[1]), "%j")))
    220L  # tropical-season centre fallback
  }

  row_list <- vector("list", length(sids))
  n_filled <- 0L

  for (k in seq_along(sids)) {
    sid  <- sids[k]
    midx <- which(ev$.sid == sid)
    if (length(midx) == 0L) next  # SID not in this location's library

    # Multiple library rows for the same SID: pick the one with highest V_peak
    if (length(midx) > 1L) {
      vpk  <- vapply(midx, function(i) get_V_peak(ev[i, ], sev), numeric(1L))
      midx <- midx[which.max(vpk)]
    }

    row        <- ev[midx, , drop = FALSE]
    doy0       <- get_doy(row)
    start_date <- as.Date(sprintf("%d-01-01", year)) + (doy0 - 1L)
    dur_days   <- get_dur_days(row)
    V_peak     <- get_V_peak(row, sev)

    counter    <- counter + 1L
    eid        <- paste0(sid, "_y", year, "_", counter)

    row_list[[k]] <- tibble::tibble(
      severity    = sev,
      start_date  = start_date,
      dur_days    = as.integer(dur_days),
      V_peak      = as.numeric(V_peak),
      event_id    = eid,
      event_class = if (sev == "HUR") "HUR" else "TS",
      Pc_min_hPa  = dplyr::coalesce(
        safe_num(row, "min_pressure_hpa"), safe_num(row, "Pc_min_hPa")),
      dP_max_hPa  = dplyr::coalesce(
        safe_num(row, "pressure_deficit_hpa"), safe_num(row, "dP_max_hPa")),
      RMW_mean_km = dplyr::coalesce(
        safe_num(row, "rmw_mean_km"), safe_num(row, "RMW_mean_km"))
    )
    n_filled <- n_filled + 1L
  }

  rows <- dplyr::bind_rows(Filter(Negate(is.null), row_list))
  if (nrow(rows) == 0L) rows <- .empty_sampled_events_tbl()
  list(rows = rows, n_filled = n_filled, counter = counter)
}


#' Generate spatially coherent daily synthetic hazard and impact time series
#'
#' @description
#' Drop-in replacement for \code{\link{generate_daily_hazard_impact}()} that
#' enforces spatial coherence across locations by sampling storms once at the
#' basin level per simulated year and assigning each drawn storm to every
#' location whose event library contains it (Option 1 — shared event pool).
#'
#' The key difference from the independent-sampling baseline is that the
#' annual storm draw happens \emph{before} the location loop rather than
#' inside it. Consequently, if Hurricane Irma (\code{"2017242N16333"}) is
#' drawn in simulated year 47, it will appear simultaneously at every location
#' that has an Irma entry in its event library (e.g. both St. Martin and Saba),
#' each retaining its own site-level wind profile, duration, and pressure
#' attributes. Storms whose tracks never came within the search radius of a
#' location are absent from that location's library and are therefore skipped
#' silently for that location, so per-location event counts may be lower than
#' the basin-level draw.
#'
#' Basin-level annual counts are resolved as \code{max(n_ts)} and
#' \code{max(n_hur)} across all requested locations for each simulation year,
#' using the counts already present in \code{out$sim}.
#'
#' All downstream processing (climate perturbation, pulse generation, damage
#' forcing) is identical to the independent-sampling variant.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param location Character vector of one or more target location names.
#' @param sim_years Integer vector of simulation-year indices to generate.
#' @param year0 Integer scalar base calendar year for \code{sim_year == 1}.
#' @param gust_factor Numeric scalar gust multiplier applied to \code{wind_kt}.
#' @param damage Named list defining the daily damage model; same specification
#'   as \code{\link{generate_daily_hazard_impact}()}.
#' @param pulse_shape Character scalar pulse shape; \code{"cosine"} or
#'   \code{"triangle"}.
#' @param scenario Optional character scalar scenario label carried into output.
#' @param seed Optional integer scalar seed. Defaults to
#'   \code{out$run_metadata$seed} or \code{1L}. Per-location library seeds are
#'   offset by location index for reproducibility.
#'
#' @return Named list of tibbles — one per requested location — with the same
#'   column schema as \code{\link{generate_daily_hazard_impact}()}.
#'
#' @seealso
#' \code{\link{generate_daily_hazard_impact}},
#' \code{\link{build_event_library_from_out}},
#' \code{\link{generate_daily_year_extended}}
#'
#' @examples
#' \dontrun{
#' daily_spatial <- generate_daily_hazard_impact_spatial(
#'   out         = hazard_out,
#'   location    = c("Saba", "St_Martin", "Statia"),
#'   sim_years   = 1:2000,
#'   year0       = 2025,
#'   gust_factor = 1.3,
#'   damage      = list(method = "intensity"),
#'   pulse_shape = "cosine",
#'   scenario    = "stationary"
#' )
#' # Verify spatial coherence: Irma should appear in the same sim_years
#' # at every location that has Irma in its event library.
#' lapply(daily_spatial, function(loc_tbl) {
#'   loc_tbl |>
#'     dplyr::filter(grepl("^2017242N16333", event_id)) |>
#'     dplyr::distinct(sim_year) |>
#'     dplyr::pull(sim_year)
#' })
#' }
#' @export
generate_daily_hazard_impact_spatial <- function(
    out,
    location,
    sim_years   = 1:1000,
    year0       = 2000,
    gust_factor = 1,
    damage      = list(method = "intensity"),
    pulse_shape = "cosine",
    scenario    = NA_character_,
    seed        = NULL,
    verbose     = TRUE) {

  stopifnot(is.character(location), length(location) >= 1L)

  if (is.null(seed)) {
    seed <- if (!is.null(out$run_metadata$seed)) out$run_metadata$seed else 1L
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be NULL or a single finite numeric value.", call. = FALSE)
  }
  seed <- as.integer(seed)

  damage <- .validate_damage_spec(damage)

  # ---- Build all location event libraries ------------------------------------
  method <- if (!is.null(out$cfg) && !is.null(out$cfg$resampling_method))
    out$cfg$resampling_method else "stratified"
  copula_min_n        <- if (!is.null(out$cfg$copula_min_n))        out$cfg$copula_min_n        else 30L
  copula_k            <- if (!is.null(out$cfg$copula_k))            out$cfg$copula_k            else 1L
  copula_robust_scale <- if (!is.null(out$cfg$copula_robust_scale)) out$cfg$copula_robust_scale else TRUE

  libs <- stats::setNames(vector("list", length(location)), location)
  for (i in seq_along(location)) {
    libs[[location[i]]] <- build_event_library_from_out(
      out,
      location            = location[i],
      seed                = seed + i - 1L,
      resampling_method   = method,
      copula_min_n        = copula_min_n,
      copula_k            = copula_k,
      copula_robust_scale = copula_robust_scale
    )
  }

  # ---- Build shared basin-level event pool -----------------------------------
  set.seed(seed)
  pool <- .build_shared_event_pool(libs)

  # ---- Validate and cache per-location sim tables ----------------------------
  if (is.null(out$sim)) stop("out$sim is NULL.", call. = FALSE)

  sim_by_loc <- stats::setNames(vector("list", length(location)), location)
  for (loc in location) {
    s <- out$sim |>
      dplyr::filter(.data$location == loc, .data$sim_year %in% sim_years)
    if (nrow(s) == 0L)
      stop("No sim years found for location '", loc, "' in out$sim.", call. = FALSE)
    sim_by_loc[[loc]] <- s
  }

  # Climate perturbation settings (forwarded unchanged from the independent variant)
  delta_sst       <- if (!is.null(out$fit)) attr(out$fit, "delta_sst") else NULL
  perturb_cfg     <- if (!is.null(out$fit)) attr(out$fit, "perturb")   else NULL
  perturb_enabled <- !is.null(perturb_cfg)

  # ---- Progress header -------------------------------------------------------
  n_years <- length(sim_years)
  t_start <- proc.time()[["elapsed"]]
  if (verbose) {
    scen_label <- if (is.na(scenario) || !nzchar(scenario)) "baseline" else scenario
    perturb_label <- if (perturb_enabled) {
      sprintf("perturb ON (delta_sst = %+.2f C)", delta_sst)
    } else {
      "perturb OFF"
    }
    message(sprintf(
      "[daily] %d years x %d location(s) | scenario: %s | seed: %d | %s",
      n_years, length(location), scen_label, seed, perturb_label
    ))
    # Quarterly progress ticks (only when the run is large enough to matter)
    tick_at <- if (n_years >= 100L) {
      unique(floor(seq(n_years / 4, n_years - 1, length.out = 3L)))
    } else {
      integer(0)
    }
  }

  # ---- Pre-allocate output storage -------------------------------------------
  results <- stats::setNames(
    lapply(location, function(loc) vector("list", length(sim_years))),
    location
  )

  # ---- Main loop: year first, then location -----------------------------------
  for (yr_idx in seq_along(sim_years)) {
    sim_yr <- sim_years[yr_idx]
    cal_yr <- as.integer(year0) + (sim_yr - 1L)

    if (verbose && yr_idx %in% tick_at) {
      message(sprintf("[daily]   %d / %d years ...", yr_idx, n_years))
    }

    # Basin-level counts: max across locations for this sim_year
    n_ts_vec <- vapply(location, function(loc) {
      s <- sim_by_loc[[loc]]
      r <- s[s$sim_year == sim_yr, ]
      if (nrow(r) == 0L) return(0L)
      as.integer(r$n_ts[1L])
    }, integer(1L))

    n_hur_vec <- vapply(location, function(loc) {
      s <- sim_by_loc[[loc]]
      r <- s[s$sim_year == sim_yr, ]
      if (nrow(r) == 0L) return(0L)
      as.integer(r$n_hur[1L])
    }, integer(1L))

    n_ts_basin  <- max(n_ts_vec,  0L)
    n_hur_basin <- max(n_hur_vec, 0L)

    # Draw shared SIDs once at basin level
    shared_ts_sids  <- .sample_shared_sids(pool, "TS",  n_ts_basin)
    shared_hur_sids <- .sample_shared_sids(pool, "HUR", n_hur_basin)

    # ---- Per-location: resolve shared SIDs and generate daily series ---------
    for (loc in location) {
      lib <- libs[[loc]]

      ts_res  <- .build_event_rows_from_sids(lib, shared_ts_sids,  "TS",  cal_yr, 0L)
      hur_res <- .build_event_rows_from_sids(lib, shared_hur_sids, "HUR", cal_yr, ts_res$counter)
      sampled <- dplyr::bind_rows(ts_res$rows, hur_res$rows)

      # Apply climate perturbation when enabled
      if (perturb_enabled && nrow(sampled) > 0L) {
        if (!is.numeric(delta_sst) || length(delta_sst) != 1L || !is.finite(delta_sst))
          stop("delta_sst must be a single finite numeric value for perturb_event().", call. = FALSE)
        sampled <- perturb_event(sampled, delta_sst = delta_sst, cc_params = perturb_cfg)
      }

      daily0 <- generate_daily_year_extended(
        year           = cal_yr,
        sampled_events = sampled,
        pulse_shape    = pulse_shape
      ) |>
        dplyr::mutate(sim_year = sim_yr, location = loc, scenario = scenario)

      # Damage intensity index
      intensity_pars <- if (identical(damage$method, "intensity")) {
        damage[c("V0", "V1", "p")]
      } else {
        list(V0 = 34, V1 = 120, p = damage$p)
      }
      daily0 <- daily0 |>
        dplyr::mutate(
          damage_intensity = pmax(
            0,
            pmin(1, (.data$wind_kt - intensity_pars$V0) /
                      (intensity_pars$V1 - intensity_pars$V0))
          )^intensity_pars$p
        )

      # Damage rate and cumulative damage
      if (damage$method == "intensity") {
        daily1 <- do.call(
          add_damage_forcing,
          c(list(daily = daily0), damage[c("V0", "V1", "p", "dmax")])
        )
        daily1 <- daily1 |>
          dplyr::mutate(cum_damage = cumsum(dplyr::coalesce(.data$damage_rate, 0)))
      } else {
        daily1 <- daily0 |>
          dplyr::mutate(
            damage_rate = do.call(
              damage_rate_from_wind,
              c(list(wind_kt = .data$wind_kt), damage[c("thr", "V_ref", "d_ref", "p", "d_max")])
            ),
            cum_damage = cumsum(dplyr::coalesce(.data$damage_rate, 0))
          )
      }

      results[[loc]][[yr_idx]] <- daily1
    }
  }

  if (verbose) {
    elapsed <- proc.time()[["elapsed"]] - t_start
    message(sprintf("[daily] Done in %.1fs", elapsed))
  }

  # ---- Assemble and return final output --------------------------------------
  for (loc in location) {
    tbl <- dplyr::bind_rows(results[[loc]]) |>
      dplyr::mutate(
        wind_gust_kt = .data$wind_kt * gust_factor,
        surge_m      = ifelse(
          is.finite(.data$pressure_hpa),
          0.14 * (1013 - .data$pressure_hpa),
          NA_real_
        )
      ) |>
      dplyr::select(
        "location", "sim_year", "scenario", "date",
        "wind_kt", "wind_gust_kt", "surge_m",
        "event_id", "event_class", "pressure_hpa",
        "pressure_deficit_hpa", "rmw_km",
        "damage_intensity", "damage_rate", "cum_damage"
      ) |>
      dplyr::relocate("location", "sim_year", "scenario", "date")
    attr(tbl, "gust_factor") <- gust_factor
    results[[loc]] <- tbl
  }

  results
}


# =============================================================================
# 6) Damage forcing
# =============================================================================

#' Add hazard intensity and damage forcing from daily wind
#'
#' @description
#' Maps daily sustained wind speed to a bounded hazard intensity index and a
#' bounded daily damage rate for downstream impact calculations.
#'
#' @param daily Tibble/data.frame with at least a \code{wind_kt} column.
#' @param V0 Numeric scalar threshold wind in kt below which intensity is zero.
#' @param V1 Numeric scalar wind in kt at which intensity saturates at one.
#' @param p Numeric scalar nonlinearity exponent.
#' @param dmax Numeric scalar maximum daily damage fraction.
#'
#' @return Tibble with added \code{damage_intensity} and \code{damage_rate}
#'   columns.
#'
#' @examples
#' daily <- tibble::tibble(wind_kt = c(20, 40, 80))
#' add_damage_forcing(daily)
#' @seealso \code{\link{damage_rate_from_wind}}
#' @export
add_damage_forcing <- function(daily,
                               V0 = 34, V1 = 120,
                               p = 3,
                               dmax = 0.02) {

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
#' Converts sustained wind speed to a bounded daily damage fraction using a
#' thresholded power-law response calibrated at \code{V_ref}.
#'
#' @param wind_kt Numeric vector of sustained wind speeds in kt.
#' @param thr Numeric scalar threshold wind in kt below which damage is zero.
#' @param V_ref Numeric scalar reference wind in kt where damage equals
#'   \code{d_ref}.
#' @param d_ref Numeric scalar damage fraction at \code{V_ref}.
#' @param p Numeric scalar exponent controlling nonlinearity.
#' @param d_max Numeric scalar upper cap on daily damage fraction.
#'
#' @return Numeric vector of daily damage rates.
#' @examples
#' damage_rate_from_wind(c(20, 40, 80))
#' @seealso \code{\link{add_damage_forcing}}
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
