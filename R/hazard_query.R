# =============================================================================
# Script overview: synthetic-year event-query helpers
# - .resolve_daily_tbl():          coerce daily input to a flat tibble.
# - .reference_storm_threshold():  derive impact threshold for a named storm.
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
# 1) Track-based query
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
#'   \code{"AL112017"} for Hurricane Irma 2017).
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
#'   storm_id = "AL112017"
#' )
#' track_years
#'
#' # Restrict to one site
#' query_storm_track_years(daily_out, storm_id = "AL112017", location = "Saba")
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
#'   \code{"AL112017"} for Hurricane Irma 2017).
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
#'   threshold derivation.
#'
#' @return Tibble with columns:
#' \describe{
#'   \item{\code{location}}{Target location name.}
#'   \item{\code{sim_year}}{Simulation year index.}
#'   \item{\code{annual_metric}}{Annual metric value for that year (column
#'     name matches the chosen \code{metric}).}
#'   \item{\code{threshold}}{The threshold applied for that location.}
#' }
#' Rows are sorted by \code{location} then \code{sim_year}.
#'
#' @examples
#' \dontrun{
#' # Years where Saba experienced at least Irma-level peak wind
#' impact_years <- query_impact_years(
#'   daily    = daily_out,
#'   storm_id = "AL112017",
#'   out      = out,
#'   location = "Saba",
#'   metric   = "peak_wind_kt"
#' )
#' impact_years
#'
#' # Cumulative damage at least as large as Irma's, with explicit threshold
#' query_impact_years(
#'   daily     = daily_out,
#'   storm_id  = "AL112017",
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
    out       = NULL,
    location  = NULL,
    metric    = c("peak_wind_kt", "cum_damage", "max_damage_rate"),
    threshold = NULL) {

  metric <- match.arg(metric)

  if (!is.character(storm_id) || length(storm_id) != 1L || !nzchar(storm_id)) {
    stop("`storm_id` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || any(!is.finite(threshold))) {
      stop("`threshold` must be a finite numeric scalar or named numeric vector.", call. = FALSE)
    }
  }

  tbl <- .resolve_daily_tbl(daily, location)

  locs <- unique(tbl$location)
  if (length(locs) == 0L) {
    stop("No locations found in the resolved daily data.", call. = FALSE)
  }

  # --- Resolve per-location thresholds ---
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

  # --- Filter years meeting threshold, joining per-location threshold values ---
  thr_df <- tibble::tibble(
    location  = names(thr_map),
    threshold = unname(thr_map)
  )

  annual |>
    dplyr::left_join(thr_df, by = "location") |>
    dplyr::filter(.data$annual_metric >= .data$threshold) |>
    dplyr::rename(!!metric := "annual_metric") |>
    dplyr::arrange(.data$location, .data$sim_year)
}
