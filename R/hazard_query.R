# =============================================================================
# Script overview: synthetic-year event-query helpers
# - .resolve_daily_tbl():          coerce daily input to a flat tibble.
# - .reference_storm_threshold():  derive impact threshold for a named storm.
# - lookup_storm_id():             find IBTrACS SIDs by year / date / intensity.
# - query_storm_track_years():     years where a specific storm's track was sampled.
# - query_impact_years():          years where site impact >= reference storm level.
# - query_aftermath_impact():      post-event window impact for candidate years.
# =============================================================================


# =============================================================================
# 0) Internal helpers
# =============================================================================

#' Coerce daily input to a single flat tibble
#'
#' @description
#' Accepts either a named list of tibbles (as returned by
#' \code{\link{generate_daily_hazard_impact}()}) or a single tibble, and
#' returns a flat tibble filtered to the requested location(s).
#'
#' @param daily Named list of tibbles or single tibble.
#' @param location Character scalar/vector or \code{NULL}. When \code{NULL}
#'   and \code{daily} is a list, all locations are bound together.
#' @return Tibble with at least the columns present in the daily schema.
#' @keywords internal
.resolve_daily_tbl <- function(daily, location = NULL) {
  if (is.data.frame(daily)) {
    tbl <- tibble::as_tibble(daily)
  } else if (is.list(daily) && !is.data.frame(daily)) {
    if (!is.null(location)) {
      missing_locs <- setdiff(location, names(daily))
      if (length(missing_locs) > 0L) {
        stop(
          "location(s) not found in daily output: ",
          paste(missing_locs, collapse = ", "),
          call. = FALSE
        )
      }
      tbl <- dplyr::bind_rows(daily[location])
    } else {
      tbl <- dplyr::bind_rows(daily)
    }
    tbl <- tibble::as_tibble(tbl)
  } else {
    stop(
      "`daily` must be a tibble or a named list from generate_daily_hazard_impact().",
      call. = FALSE
    )
  }

  # Apply location filter when daily was a plain tibble and location is specified
  if (is.data.frame(daily) && !is.null(location)) {
    if ("location" %in% names(tbl)) {
      tbl <- dplyr::filter(tbl, .data$location %in% !!location)
    }
  }

  if (nrow(tbl) == 0L) {
    warning("No rows remain after resolving daily input and location filter.", call. = FALSE)
  }
  tbl
}


#' Derive an impact threshold for a reference storm
#'
#' @description
#' Resolves a numeric threshold for a named metric from either:
#' \enumerate{
#'   \item The historical event record in \code{out$events} (for
#'         \code{metric = "peak_wind_kt"}).
#'   \item The rows of \code{daily} where that storm appears as the sampled
#'         event (for damage metrics), taking the median annual metric across
#'         those occurrences.
#' }
#'
#' @param daily Flat tibble (already resolved by \code{.resolve_daily_tbl()}).
#' @param out List returned by \code{run_hazard_model()}.
#' @param storm_id Character scalar IBTrACS SID of the reference storm.
#' @param metric Character scalar metric name; one of \code{"peak_wind_kt"},
#'   \code{"cum_damage"}, or \code{"max_damage_rate"}.
#' @param location Character vector of locations to consider, or \code{NULL}
#'   for all.
#' @return Named numeric vector of thresholds, one entry per location.
#' @keywords internal
.reference_storm_threshold <- function(daily, out, storm_id, metric, location = NULL) {

  locs <- if (!is.null(location)) location else unique(daily$location)
  if (length(locs) == 0L) stop("No locations found in daily data.", call. = FALSE)

  thresholds <- stats::setNames(vector("numeric", length(locs)), locs)

  for (loc in locs) {
    thr <- NA_real_

    if (metric == "peak_wind_kt") {
      # Prefer historical record: storm's site-level peak wind at this location
      if (!is.null(out$events) && all(c("storm_id", "location", "peak_wind_kt") %in% names(out$events))) {
        ev_row <- dplyr::filter(
          out$events,
          .data$storm_id == !!storm_id,
          .data$location == !!loc
        )
        if (nrow(ev_row) > 0L) {
          thr <- max(ev_row$peak_wind_kt, na.rm = TRUE)
        }
      }
    }

    # Fallback for all metrics (including damage metrics that are absent from
    # out$events): find storm occurrences in the simulation and take the
    # median annual metric across those years.
    if (!is.finite(thr)) {
      loc_daily <- dplyr::filter(daily, .data$location == !!loc)
      storm_days <- dplyr::filter(
        loc_daily,
        !is.na(.data$event_id),
        grepl(paste0("^", storm_id, "_y"), .data$event_id, fixed = FALSE)
      )

      if (nrow(storm_days) == 0L) {
        stop(
          "Storm '", storm_id, "' not found in daily output for location '", loc, "'. ",
          "Supply an explicit `threshold` value.",
          call. = FALSE
        )
      }

      storm_years <- unique(storm_days$sim_year)
      annual_tbl  <- .annual_metric_tbl(loc_daily, metric)
      ref_vals    <- annual_tbl$annual_metric[annual_tbl$sim_year %in% storm_years]

      thr <- if (length(ref_vals) > 0L) stats::median(ref_vals, na.rm = TRUE) else NA_real_
    }

    if (!is.finite(thr)) {
      stop(
        "Could not derive a finite threshold for storm '", storm_id,
        "' at location '", loc, "' for metric '", metric, "'. ",
        "Supply an explicit `threshold` value.",
        call. = FALSE
      )
    }

    thresholds[[loc]] <- thr
  }

  thresholds
}


#' Compute per-sim_year aggregate metric for one location's daily series
#'
#' @param loc_daily Tibble of daily rows for a single location.
#' @param metric Character scalar; one of \code{"peak_wind_kt"},
#'   \code{"cum_damage"}, or \code{"max_damage_rate"}.
#' @return Tibble with columns \code{location}, \code{sim_year},
#'   \code{annual_metric}.
#' @keywords internal
.annual_metric_tbl <- function(loc_daily, metric) {
  if (metric == "peak_wind_kt") {
    if (!("wind_kt" %in% names(loc_daily))) {
      stop("`daily` must contain a `wind_kt` column for metric 'peak_wind_kt'.", call. = FALSE)
    }
    loc_daily |>
      dplyr::group_by(.data$location, .data$sim_year) |>
      dplyr::summarise(annual_metric = max(.data$wind_kt, na.rm = TRUE), .groups = "drop")

  } else if (metric == "cum_damage") {
    if (!("cum_damage" %in% names(loc_daily))) {
      stop("`daily` must contain a `cum_damage` column for metric 'cum_damage'.", call. = FALSE)
    }
    loc_daily |>
      dplyr::group_by(.data$location, .data$sim_year) |>
      dplyr::summarise(annual_metric = max(.data$cum_damage, na.rm = TRUE), .groups = "drop")

  } else if (metric == "max_damage_rate") {
    if (!("damage_rate" %in% names(loc_daily))) {
      stop("`daily` must contain a `damage_rate` column for metric 'max_damage_rate'.", call. = FALSE)
    }
    loc_daily |>
      dplyr::group_by(.data$location, .data$sim_year) |>
      dplyr::summarise(annual_metric = max(.data$damage_rate, na.rm = TRUE), .groups = "drop")

  } else {
    stop("Unknown metric '", metric, "'. Use 'peak_wind_kt', 'cum_damage', or 'max_damage_rate'.",
         call. = FALSE)
  }
}


# =============================================================================
# 1) Storm ID lookup
# =============================================================================

#' Look up IBTrACS storm IDs from the hazard model event record
#'
#' @description
#' The hazard model stores storms under their \strong{IBTrACS native SID}
#' (format: \code{YYYYDDDLLLBBB}, e.g. \code{"2017242N16333"} for Hurricane
#' Irma), \emph{not} the ATCF identifier familiar from NHC advisories
#' (\code{"2017242N16333"}).  This function filters \code{out$events} so you can
#' find the correct \code{storm_id} to pass to \code{\link{query_storm_track_years}()}
#' and \code{\link{query_impact_years}()} using human-readable criteria.
#'
#' All filter arguments are optional and combined with AND logic.
#' Omitting all arguments returns the full event record.
#'
#' @param out List returned by \code{\link{run_hazard_model}()}.
#' @param year Integer scalar or vector: calendar season year(s) to include.
#' @param location Character scalar or vector: location name(s) to filter on.
#'   \code{NULL} aggregates across locations (returns distinct storm IDs).
#' @param date_range Length-2 vector of \code{Date} or character scalars
#'   (\code{c("2017-09-01", "2017-09-15")}) constraining \code{start_time}.
#' @param min_wind_kt Numeric scalar: minimum \code{peak_wind_kt} at the site.
#'
#' @return Tibble of matching rows from \code{out$events} with columns
#'   \code{storm_id}, \code{location}, \code{start_time}, \code{peak_wind_kt},
#'   and \code{storm_intensity_kt}, sorted by \code{start_time}.
#'   Use \code{storm_id} values directly in subsequent query functions.
#'
#' @examples
#' \dontrun{
#' # Find Irma by approximate date and intensity at Saba
#' lookup_storm_id(
#'   out          = hazard_out,
#'   year         = 2017,
#'   location     = "Saba",
#'   date_range   = c("2017-09-01", "2017-09-15"),
#'   min_wind_kt  = 50
#' )
#' # Returns storm_id "2017242N16333" — pass this to query_storm_track_years()
#'
#' # All major hurricanes (>= 96 kt) at any location
#' lookup_storm_id(out = hazard_out, min_wind_kt = 96)
#' }
#' @seealso \code{\link{query_storm_track_years}}, \code{\link{query_impact_years}}
#' @export
lookup_storm_id <- function(
    out,
    year        = NULL,
    location    = NULL,
    date_range  = NULL,
    min_wind_kt = NULL) {

  if (is.null(out$events) || !is.data.frame(out$events)) {
    stop("`out$events` not found. Pass the list returned by run_hazard_model().",
         call. = FALSE)
  }

  required <- c("storm_id", "location", "start_time", "peak_wind_kt", "storm_intensity_kt")
  missing_cols <- setdiff(required, names(out$events))
  if (length(missing_cols) > 0L) {
    stop(
      "`out$events` is missing expected columns: ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }

  ev <- out$events

  # --- Apply filters ----------------------------------------------------------
  if (!is.null(year)) {
    if (!("year" %in% names(ev))) {
      ev$year <- as.integer(format(as.Date(ev$start_time), "%Y"))
    }
    ev <- ev[ev$year %in% as.integer(year), , drop = FALSE]
  }

  if (!is.null(location)) {
    ev <- ev[ev$location %in% location, , drop = FALSE]
  }

  if (!is.null(date_range)) {
    if (length(date_range) != 2L) {
      stop("`date_range` must be a length-2 vector: c(start, end).", call. = FALSE)
    }
    d0 <- as.Date(date_range[1])
    d1 <- as.Date(date_range[2])
    ev_dates <- as.Date(ev$start_time)
    ev <- ev[!is.na(ev_dates) & ev_dates >= d0 & ev_dates <= d1, , drop = FALSE]
  }

  if (!is.null(min_wind_kt)) {
    ev <- ev[!is.na(ev$peak_wind_kt) & ev$peak_wind_kt >= min_wind_kt, , drop = FALSE]
  }

  if (nrow(ev) == 0L) {
    message("No storms match the supplied filters.")
    return(tibble::tibble(
      storm_id           = character(0),
      location           = character(0),
      start_time         = as.POSIXct(character(0)),
      peak_wind_kt       = numeric(0),
      storm_intensity_kt = numeric(0)
    ))
  }

  ev |>
    dplyr::select(
      "storm_id", "location", "start_time",
      "peak_wind_kt", "storm_intensity_kt"
    ) |>
    dplyr::arrange(.data$start_time) |>
    tibble::as_tibble()
}


# =============================================================================
# 2) Track-based query
# =============================================================================

#' Find synthetic years where a specific storm's track was sampled
#'
#' @description
#' Scans the daily hazard-impact output and returns the simulation years in
#' which a named historical storm was drawn from the event library, preserving
#' its original approach geometry, temporal wind profile, and intensity.
#'
#' The event-library sampler encodes the source storm in the \code{event_id}
#' column as \code{"<SID>_y<year>_<counter>"}. This function matches rows
#' where \code{event_id} begins with \code{"<storm_id>_y"}, so any simulation
#' year that received that track will be returned.
#'
#' @param daily Named list of tibbles returned by
#'   \code{\link{generate_daily_hazard_impact}()}, or a single tibble for one
#'   location.
#' @param storm_id Character scalar IBTrACS SID of the target storm (e.g.
#'   \code{"2017242N16333"} for Hurricane Irma — IBTrACS native SID, not ATCF "AL112017").
#' @param location Character scalar or vector of location names to query, or
#'   \code{NULL} to query all locations present in \code{daily}.
#'
#' @return Tibble with columns:
#' \describe{
#'   \item{\code{location}}{Target location name.}
#'   \item{\code{sim_year}}{Synthetic simulation-year index in which the storm
#'     was sampled.}
#' }
#' Rows are sorted by \code{location} then \code{sim_year}. An empty tibble
#' is returned (with the same schema) when the storm never appears.
#'
#' @examples
#' \dontrun{
#' # Find all synthetic years where Hurricane Irma's track was resampled
#' track_years <- query_storm_track_years(
#'   daily    = daily_out,
#'   storm_id = "2017242N16333"
#' )
#' track_years
#'
#' # Restrict to one site
#' query_storm_track_years(daily_out, storm_id = "2017242N16333", location = "Saba")
#' }
#' @seealso \code{\link{query_impact_years}},
#'   \code{\link{generate_daily_hazard_impact}}
#' @export
query_storm_track_years <- function(daily, storm_id, location = NULL) {
  if (!is.character(storm_id) || length(storm_id) != 1L || !nzchar(storm_id)) {
    stop("`storm_id` must be a single non-empty character string.", call. = FALSE)
  }

  tbl <- .resolve_daily_tbl(daily, location)

  if (!("event_id" %in% names(tbl))) {
    stop("`daily` must contain an `event_id` column.", call. = FALSE)
  }
  if (!("sim_year" %in% names(tbl))) {
    stop("`daily` must contain a `sim_year` column.", call. = FALSE)
  }

  # event_id pattern: "<SID>_y<calendar_year>_<counter>"
  pattern <- paste0("^", storm_id, "_y")

  matched <- dplyr::filter(
    tbl,
    !is.na(.data$event_id),
    grepl(pattern, .data$event_id)
  )

  if (nrow(matched) == 0L) {
    message(
      "Storm '", storm_id, "' not found in any simulated year",
      if (!is.null(location)) paste0(" for location(s): ", paste(location, collapse = ", ")) else "",
      "."
    )
    return(tibble::tibble(location = character(0), sim_year = integer(0)))
  }

  matched |>
    dplyr::select("location", "sim_year") |>
    dplyr::distinct() |>
    dplyr::arrange(.data$location, .data$sim_year)
}


# =============================================================================
# 2) Impact-based query
# =============================================================================

#' Find synthetic years with at least a reference storm's impact level
#'
#' @description
#' Returns every simulation year in which the site's annual impact metric
#' equalled or exceeded the level caused by a named reference storm, regardless
#' of which storm(s) were responsible.
#'
#' The reference impact threshold is resolved in the following order of
#' priority:
#' \enumerate{
#'   \item An explicit \code{threshold} value supplied by the caller.
#'   \item For \code{metric = "peak_wind_kt"}: the reference storm's
#'     site-level peak wind taken directly from \code{out$events}.
#'   \item For all metrics (and as a fallback): the median annual metric
#'     across the simulation years in which that storm appeared in \code{daily},
#'     as identified by \code{\link{query_storm_track_years}()}.
#' }
#'
#' @param daily Named list of tibbles returned by
#'   \code{\link{generate_daily_hazard_impact}()}, or a single tibble.
#' @param storm_id Character scalar IBTrACS SID of the reference storm (e.g.
#'   \code{"2017242N16333"} for Hurricane Irma — IBTrACS native SID, not ATCF "AL112017").
#' @param out Optional list returned by \code{\link{run_hazard_model}()}. Used
#'   to look up the reference storm's historical peak wind when
#'   \code{metric = "peak_wind_kt"} and \code{threshold = NULL}.
#' @param location Character scalar or vector of location names to query, or
#'   \code{NULL} for all locations in \code{daily}.
#' @param metric Character scalar defining how annual impact is measured.
#'   \describe{
#'     \item{\code{"peak_wind_kt"}}{Maximum daily sustained wind in the year
#'       (kt). Default.}
#'     \item{\code{"cum_damage"}}{Cumulative annual damage fraction (final
#'       value of the \code{cum_damage} column).}
#'     \item{\code{"max_damage_rate"}}{Maximum daily damage rate in the year.}
#'   }
#' @param threshold Optional numeric scalar or named numeric vector. When a
#'   scalar is supplied it is applied to all locations. When a named vector is
#'   supplied the names must match the queried locations. Overrides automatic
#'   storm-based threshold derivation; used as the lower bound when
#'   \code{percentile} is also specified.
#' @param percentile Optional numeric scalar in \code{(0, 1)}. When supplied,
#'   the empirical \code{percentile} of the annual metric distribution is
#'   computed per location and used as an additional selection gate. The
#'   effective threshold becomes
#'   \code{max(ref_threshold, percentile_threshold)}, so only years that clear
#'   \emph{both} the reference-storm level and the distributional percentile
#'   are returned. For example \code{percentile = 0.95} retains at most the top
#'   5 \% of years. Default \code{NULL} disables the percentile gate.
#' @param min_threshold Optional numeric scalar applied as an absolute floor on
#'   the effective threshold after all other computations. Useful for ensuring a
#'   minimum physical severity level regardless of the reference storm or
#'   percentile gate (e.g. \code{min_threshold = 64} for \code{metric =
#'   "peak_wind_kt"} to require at least Cat-1 strength). Default \code{NULL}.
#'
#' @return Tibble with columns:
#' \describe{
#'   \item{\code{location}}{Target location name.}
#'   \item{\code{sim_year}}{Simulation year index.}
#'   \item{\code{annual_metric}}{Annual metric value for that year (column
#'     name matches the chosen \code{metric}).}
#'   \item{\code{threshold}}{The effective threshold applied for that location,
#'     reflecting the combined result of storm-based, percentile, and
#'     \code{min_threshold} logic.}
#' }
#' Rows are sorted by \code{location} then \code{sim_year}.
#'
#' @examples
#' \dontrun{
#' # Years where Saba experienced at least Irma-level peak wind
#' impact_years <- query_impact_years(
#'   daily    = daily_out,
#'   storm_id = "2017242N16333",
#'   out      = out,
#'   location = "Saba",
#'   metric   = "peak_wind_kt"
#' )
#' impact_years
#'
#' # Hybrid: Irma wind as lower bound, 95th-percentile as selection gate
#' query_impact_years(
#'   daily      = daily_out,
#'   storm_id   = "2017242N16333",
#'   out        = out,
#'   metric     = "peak_wind_kt",
#'   percentile = 0.95
#' )
#'
#' # Hybrid with explicit Cat-2 floor (83 kt) regardless of Irma's site wind
#' query_impact_years(
#'   daily         = daily_out,
#'   storm_id      = "2017242N16333",
#'   out           = out,
#'   metric        = "peak_wind_kt",
#'   percentile    = 0.95,
#'   min_threshold = 83
#' )
#'
#' # Cumulative damage at least as large as Irma's, with explicit threshold
#' query_impact_years(
#'   daily     = daily_out,
#'   storm_id  = "2017242N16333",
#'   location  = "Saba",
#'   metric    = "cum_damage",
#'   threshold = 0.15
#' )
#' }
#' @seealso \code{\link{query_storm_track_years}},
#'   \code{\link{generate_daily_hazard_impact}}
#' @export
query_impact_years <- function(
    daily,
    storm_id,
    out           = NULL,
    location      = NULL,
    metric        = c("peak_wind_kt", "cum_damage", "max_damage_rate"),
    threshold     = NULL,
    percentile    = NULL,
    min_threshold = NULL) {

  metric <- match.arg(metric)

  if (!is.character(storm_id) || length(storm_id) != 1L || !nzchar(storm_id)) {
    stop("`storm_id` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || any(!is.finite(threshold))) {
      stop("`threshold` must be a finite numeric scalar or named numeric vector.", call. = FALSE)
    }
  }

  if (!is.null(percentile)) {
    if (!is.numeric(percentile) || length(percentile) != 1L ||
        !is.finite(percentile) || percentile <= 0 || percentile >= 1) {
      stop("`percentile` must be a single numeric value in (0, 1).", call. = FALSE)
    }
  }

  if (!is.null(min_threshold)) {
    if (!is.numeric(min_threshold) || length(min_threshold) != 1L ||
        !is.finite(min_threshold)) {
      stop("`min_threshold` must be a single finite numeric value.", call. = FALSE)
    }
  }

  tbl <- .resolve_daily_tbl(daily, location)

  locs <- unique(tbl$location)
  if (length(locs) == 0L) {
    stop("No locations found in the resolved daily data.", call. = FALSE)
  }

  # --- Resolve per-location reference thresholds ---
  if (!is.null(threshold)) {
    if (length(threshold) == 1L) {
      thr_map <- stats::setNames(rep(threshold, length(locs)), locs)
    } else {
      missing_names <- setdiff(locs, names(threshold))
      if (length(missing_names) > 0L) {
        stop(
          "Named `threshold` is missing entries for location(s): ",
          paste(missing_names, collapse = ", "),
          call. = FALSE
        )
      }
      thr_map <- threshold[locs]
    }
  } else {
    thr_map <- .reference_storm_threshold(
      daily    = tbl,
      out      = out,
      storm_id = storm_id,
      metric   = metric,
      location = locs
    )
  }

  # --- Compute annual metric across all locations ---
  annual <- .annual_metric_tbl(tbl, metric)

  # --- Build threshold table, applying percentile gate and min_threshold ------
  thr_df <- tibble::tibble(
    location  = names(thr_map),
    threshold = unname(thr_map)
  )

  # Percentile gate: compute per-location percentile and take the max with the
  # reference threshold so that both criteria must be satisfied simultaneously.
  if (!is.null(percentile)) {
    pct_df <- annual |>
      dplyr::group_by(.data$location) |>
      dplyr::summarise(
        pct_threshold = stats::quantile(.data$annual_metric, percentile, na.rm = TRUE),
        .groups = "drop"
      )
    thr_df <- thr_df |>
      dplyr::left_join(pct_df, by = "location") |>
      dplyr::mutate(
        threshold = pmax(.data$threshold, .data$pct_threshold, na.rm = TRUE)
      ) |>
      dplyr::select("location", "threshold")
  }

  # Absolute floor: applied after all other threshold logic.
  if (!is.null(min_threshold)) {
    thr_df <- thr_df |>
      dplyr::mutate(threshold = pmax(.data$threshold, min_threshold))
  }

  # --- Filter years meeting the effective threshold ---------------------------
  annual |>
    dplyr::left_join(thr_df, by = "location") |>
    dplyr::filter(.data$annual_metric >= .data$threshold) |>
    dplyr::rename(!!metric := "annual_metric") |>
    dplyr::arrange(.data$location, .data$sim_year)
}


# =============================================================================
# 3) Aftermath query
# =============================================================================

#' Measure post-event impact in the days following a primary storm
#'
#' @description
#' For each candidate \code{(location, sim_year)}, locates the end of the
#' primary event, then summarises the \code{window_days} days that follow.
#' This captures compound risk: recovery-period storms, re-intensification,
#' or successive hits that arrive while the site is already damaged.
#'
#' \strong{Primary-event anchor:}
#' \describe{
#'   \item{With \code{storm_id}}{The last calendar day on which an event
#'     whose \code{event_id} starts with \code{"<storm_id>_y"} is active.
#'     Simulation years where that storm was not sampled are silently dropped.}
#'   \item{Without \code{storm_id} (\code{NULL})}{The last day of the event
#'     that produced the year's peak wind (\code{max(wind_kt)}).  Calm years
#'     (no TC activity) anchor on the day of peak wind even when
#'     \code{event_id} is \code{NA}.}
#' }
#' The aftermath window runs from \code{anchor_date + 1} through
#' \code{anchor_date + window_days}, clipped to the end of the simulated year.
#'
#' \strong{Aftermath damage:}
#' \code{aftermath_cum_damage} is the additional cumulative damage fraction
#' accrued during the window, computed as
#' \eqn{\max(\texttt{cum\_damage}) - \texttt{cum\_damage}[\texttt{anchor\_date}]}.
#' Because \code{cum_damage} is monotonically non-decreasing within a year,
#' this delta is always \eqn{\ge 0}.
#'
#' @param daily Named list of tibbles returned by
#'   \code{\link{generate_daily_hazard_impact}()}, or a single tibble.
#' @param sim_years Optional filter applied before computing aftermath.
#'   Either:
#'   \itemize{
#'     \item An integer vector of \code{sim_year} values applied to all
#'           locations.
#'     \item A tibble with \code{location} and \code{sim_year} columns (the
#'           direct output of \code{\link{query_storm_track_years}()} or
#'           \code{\link{query_impact_years}()}).
#'     \item \code{NULL} (default) — all years are included.
#'   }
#' @param storm_id Character scalar IBTrACS SID used to anchor the aftermath
#'   window (e.g. \code{"2017242N16333"} for Irma).  When \code{NULL} the window is
#'   anchored on the year's peak-wind event.
#' @param window_days Positive integer: length of the aftermath window in
#'   calendar days (default \code{60L}).
#' @param location Character scalar or vector of location names, or \code{NULL}
#'   for all locations.
#'
#' @return Tibble with one row per \code{(location, sim_year)} and columns:
#' \describe{
#'   \item{\code{location}}{Location name.}
#'   \item{\code{sim_year}}{Simulation-year index.}
#'   \item{\code{primary_event_id}}{The \code{event_id} of the anchoring event
#'     (\code{NA} for calm years when \code{storm_id = NULL}).}
#'   \item{\code{event_end_date}}{Last calendar date of the primary event;
#'     the aftermath window starts the following day.}
#'   \item{\code{aftermath_peak_wind_kt}}{Maximum sustained wind (kt) in the
#'     window.}
#'   \item{\code{aftermath_n_events}}{Count of distinct storm events active
#'     during the window (excluding the primary event).}
#'   \item{\code{aftermath_n_hur_days}}{Days with sustained wind \eqn{\ge}
#'     64 kt during the window.}
#'   \item{\code{aftermath_cum_damage}}{Cumulative damage fraction accrued
#'     during the window (delta from primary event end).}
#'   \item{\code{aftermath_max_damage_rate}}{Maximum single-day damage rate
#'     during the window.}
#'   \item{\code{aftermath_rank}}{Integer rank within location;
#'     1 = worst aftermath by \code{aftermath_cum_damage}.}
#' }
#' Rows are sorted by \code{location} then \code{aftermath_rank}.
#'
#' @examples
#' \dontrun{
#' # Aftermath in the 60 days following an Irma-track event
#' track_years  <- query_storm_track_years(daily_out, storm_id = "2017242N16333")
#' aftermath_60 <- query_aftermath_impact(
#'   daily       = daily_out,
#'   sim_years   = track_years,
#'   storm_id    = "2017242N16333",
#'   window_days = 60L
#' )
#'
#' # 90-day aftermath anchored on each year's own peak-wind event
#' impact_years  <- query_impact_years(daily_out, storm_id = "2017242N16333",
#'                                      out = hazard_out)
#' aftermath_90  <- query_aftermath_impact(
#'   daily       = daily_out,
#'   sim_years   = impact_years,
#'   storm_id    = NULL,
#'   window_days = 90L
#' )
#' }
#' @seealso \code{\link{query_storm_track_years}}, \code{\link{query_impact_years}},
#'   \code{\link{compute_stress_year_metrics}}
#' @export
query_aftermath_impact <- function(
    daily,
    sim_years   = NULL,
    storm_id    = NULL,
    window_days = 60L,
    location    = NULL) {

  # --- Input validation -------------------------------------------------------
  window_days <- as.integer(window_days)
  if (!is.finite(window_days) || window_days < 1L) {
    stop("`window_days` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(storm_id)) {
    if (!is.character(storm_id) || length(storm_id) != 1L || !nzchar(storm_id)) {
      stop("`storm_id` must be a single non-empty character string.", call. = FALSE)
    }
  }

  # --- Resolve and filter daily data ------------------------------------------
  tbl <- .resolve_daily_tbl(daily, location)

  required <- c("date", "sim_year", "wind_kt", "event_id", "cum_damage", "damage_rate")
  missing_cols <- setdiff(required, names(tbl))
  if (length(missing_cols) > 0L) {
    stop(
      "`daily` is missing required columns: ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  tbl$date <- as.Date(tbl$date)

  # Apply sim_years filter (.resolve_sim_years_filter is defined in
  # hazard_stress_select.R; both live in the same package namespace)
  yr_filter <- .resolve_sim_years_filter(sim_years, unique(tbl$location))
  if (!is.null(yr_filter)) {
    tbl <- dplyr::bind_rows(lapply(unique(tbl$location), function(loc) {
      yrs <- yr_filter[[loc]]
      if (length(yrs) == 0L) return(tbl[0L, ])
      dplyr::filter(tbl, .data$location == loc, .data$sim_year %in% yrs)
    }))
  }

  if (nrow(tbl) == 0L) {
    warning("No rows remain after applying sim_years filter.", call. = FALSE)
    return(.empty_aftermath_tbl())
  }

  # --- Iterate over (location, sim_year) pairs --------------------------------
  yr_locs     <- dplyr::distinct(tbl, .data$location, .data$sim_year)
  result_list <- vector("list", nrow(yr_locs))

  if (!is.null(storm_id)) {
    primary_pattern <- paste0("^", storm_id, "_y")
  }

  for (i in seq_len(nrow(yr_locs))) {
    loc  <- yr_locs$location[i]
    yr   <- yr_locs$sim_year[i]
    rows <- tbl[tbl$location == loc & tbl$sim_year == yr, , drop = FALSE]
    rows <- rows[order(rows$date), , drop = FALSE]

    # --- Locate primary event and derive anchor date -------------------------
    if (!is.null(storm_id)) {
      primary_rows <- rows[
        !is.na(rows$event_id) & grepl(primary_pattern, rows$event_id),
        , drop = FALSE
      ]
      if (nrow(primary_rows) == 0L) {
        # Storm not sampled in this year; drop silently
        result_list[[i]] <- NULL
        next
      }
      event_end_date   <- max(primary_rows$date)
      # Use the event_id from the last day of the primary event for attribution
      primary_event_id <- primary_rows$event_id[which.max(primary_rows$date)]

    } else {
      # Anchor on the last day of the peak-wind event
      peak_idx         <- which.max(rows$wind_kt)
      peak_event_id    <- rows$event_id[peak_idx]
      if (!is.na(peak_event_id)) {
        event_rows     <- rows[!is.na(rows$event_id) & rows$event_id == peak_event_id, , drop = FALSE]
        event_end_date <- max(event_rows$date)
      } else {
        event_end_date <- rows$date[peak_idx]
      }
      primary_event_id <- peak_event_id   # NA for calm years
    }

    # --- Extract aftermath window --------------------------------------------
    window_start   <- event_end_date + 1L
    window_end     <- event_end_date + window_days
    aftermath_rows <- rows[rows$date >= window_start & rows$date <= window_end, , drop = FALSE]

    # Cumulative damage at anchor date (last value on or before event_end_date)
    pre_rows      <- rows[rows$date <= event_end_date, , drop = FALSE]
    damage_anchor <- if (nrow(pre_rows) > 0L) {
      max(pre_rows$cum_damage, na.rm = TRUE)
    } else {
      0
    }
    damage_anchor <- if (is.finite(damage_anchor)) damage_anchor else 0

    if (nrow(aftermath_rows) == 0L) {
      # Primary event runs to the end of the simulated year — no days remain
      result_list[[i]] <- tibble::tibble(
        location                  = loc,
        sim_year                  = as.integer(yr),
        primary_event_id          = primary_event_id,
        event_end_date            = event_end_date,
        aftermath_peak_wind_kt    = 0,
        aftermath_n_events        = 0L,
        aftermath_n_hur_days      = 0L,
        aftermath_cum_damage      = 0,
        aftermath_max_damage_rate = 0
      )
      next
    }

    # Subsequent events: distinct event_ids in the window, excluding the primary
    subsequent_event_ids <- unique(
      aftermath_rows$event_id[!is.na(aftermath_rows$event_id)]
    )
    if (!is.null(storm_id)) {
      subsequent_event_ids <- subsequent_event_ids[
        !grepl(primary_pattern, subsequent_event_ids)
      ]
    }

    damage_window <- max(aftermath_rows$cum_damage, na.rm = TRUE)
    damage_window <- if (is.finite(damage_window)) damage_window else damage_anchor

    result_list[[i]] <- tibble::tibble(
      location                  = loc,
      sim_year                  = as.integer(yr),
      primary_event_id          = primary_event_id,
      event_end_date            = event_end_date,
      aftermath_peak_wind_kt    = max(aftermath_rows$wind_kt,      na.rm = TRUE),
      aftermath_n_events        = as.integer(length(subsequent_event_ids)),
      aftermath_n_hur_days      = as.integer(sum(aftermath_rows$wind_kt >= 64, na.rm = TRUE)),
      aftermath_cum_damage      = pmax(0, damage_window - damage_anchor),
      aftermath_max_damage_rate = max(aftermath_rows$damage_rate,  na.rm = TRUE)
    )
  }

  out <- dplyr::bind_rows(result_list)

  if (nrow(out) == 0L) {
    return(.empty_aftermath_tbl())
  }

  # Replace non-finite values (calm no-TC windows produce -Inf from max())
  out <- out |>
    dplyr::mutate(
      aftermath_peak_wind_kt    = dplyr::if_else(
        is.finite(.data$aftermath_peak_wind_kt),    .data$aftermath_peak_wind_kt,    0),
      aftermath_cum_damage      = dplyr::if_else(
        is.finite(.data$aftermath_cum_damage),      .data$aftermath_cum_damage,      0),
      aftermath_max_damage_rate = dplyr::if_else(
        is.finite(.data$aftermath_max_damage_rate), .data$aftermath_max_damage_rate, 0)
    )

  # Rank within location: 1 = worst aftermath by cumulative damage
  out |>
    dplyr::group_by(.data$location) |>
    dplyr::mutate(
      aftermath_rank = as.integer(rank(-.data$aftermath_cum_damage, ties.method = "first"))
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$location, .data$aftermath_rank)
}


#' Empty aftermath tibble with correct schema
#' @keywords internal
.empty_aftermath_tbl <- function() {
  tibble::tibble(
    location                  = character(0),
    sim_year                  = integer(0),
    primary_event_id          = character(0),
    event_end_date            = as.Date(character(0)),
    aftermath_peak_wind_kt    = numeric(0),
    aftermath_n_events        = integer(0),
    aftermath_n_hur_days      = integer(0),
    aftermath_cum_damage      = numeric(0),
    aftermath_max_damage_rate = numeric(0),
    aftermath_rank            = integer(0)
  )
}
