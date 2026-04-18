# =============================================================================
# Script overview: stress-year characterisation and representative selection
# - compute_stress_year_metrics(): per-year multi-metric summary table.
# - aggregate_stress_metrics():    weighted composite score from chosen metrics.
# - select_stress_years():         k representative years with good sample coverage.
#
# Typical workflow — portfolio mode (single set of years across all locations)
# ─────────────────────────────────────────────────────────────────────────────
#   daily_out  <- generate_daily_hazard_impact_spatial(...)          # full ensemble
#   filtered   <- query_impact_years(daily_out, ...)         # candidate set
#
#   metrics    <- compute_stress_year_metrics(               # step 1: per (location, sim_year)
#                   daily_out,
#                   sim_years = filtered)
#
#   scored     <- aggregate_stress_metrics(                  # step 2: aggregate across locations
#                   metrics,                                 # → one row per sim_year
#                   metrics_used = c("peak_wind_kt", "cum_damage", "n_events"),
#                   weights      = c(2, 1, 0.5),
#                   location_agg = "mean")                   # "mean" | "max" | "sum"
#
#   selected   <- select_stress_years(scored, k = 10,        # step 3: single portfolio set
#                   method = "stratified")                   # no location column in output
# =============================================================================


# =============================================================================
# 0) Internal helpers
# =============================================================================

#' Resolve sim_years filter to a per-location list
#'
#' @description
#' Accepts either (a) an integer vector applied uniformly across all locations,
#' or (b) a tibble with at least \code{location} and \code{sim_year} columns
#' (the output of \code{query_storm_track_years()} or
#' \code{query_impact_years()}) and returns a named list mapping each location
#' to its allowed \code{sim_year} values.
#'
#' @param sim_years Integer vector, tibble with \code{location}/\code{sim_year}
#'   columns, or \code{NULL} (no filter).
#' @param locations Character vector of all locations present in \code{daily}.
#' @return Named list: \code{list(loc = integer_vector, ...)} or \code{NULL}.
#' @keywords internal
.resolve_sim_years_filter <- function(sim_years, locations) {
  if (is.null(sim_years)) return(NULL)

  if (is.data.frame(sim_years)) {
    if (!all(c("location", "sim_year") %in% names(sim_years))) {
      stop(
        "When `sim_years` is a tibble it must contain `location` and `sim_year` columns.",
        call. = FALSE
      )
    }
    # Build per-location lookup from the tibble
    out <- stats::setNames(vector("list", length(locations)), locations)
    for (loc in locations) {
      yrs <- sim_years$sim_year[sim_years$location == loc]
      out[[loc]] <- if (length(yrs) > 0L) as.integer(yrs) else integer(0)
    }
    return(out)
  }

  # Plain integer vector — same years for every location
  yrs <- as.integer(sim_years)
  stats::setNames(replicate(length(locations), yrs, simplify = FALSE), locations)
}


#' Compute compound-event stress metrics from a focal event and 60-day aftermath
#'
#' @param daily Named list of tibbles returned by
#'   \code{\link{generate_daily_hazard_impact_spatial}()}, or a single tibble.
#' @param sim_years Optional filter. Either an integer vector of
#'   \code{sim_year} values applied to all locations, a tibble with
#'   \code{location} and \code{sim_year} columns, or \code{NULL}.
#' @param location Character vector of location names to include.
#' @param window_days Positive integer: length of the compound window after the
#'   focal event end date.
#'
#' @return Tibble with one row per \code{sim_year}, describing the strongest
#'   focal event across the selected locations together and the cumulative
#'   damage accrued from focal-event onset through the aftermath window.
#'
#' @keywords internal
#' @noRd
compute_compound_stress_year_metrics <- function(
    daily,
    sim_years = NULL,
    location,
    window_days = 60L) {

  if (missing(location) || is.null(location) || length(location) == 0L) {
    stop("`location` must name at least one location.", call. = FALSE)
  }
  if (!is.character(location) || any(is.na(location)) || any(!nzchar(location))) {
    stop("`location` must be a character vector of non-empty names.", call. = FALSE)
  }

  window_days <- as.integer(window_days)
  if (!is.finite(window_days) || window_days < 1L) {
    stop("`window_days` must be a positive integer.", call. = FALSE)
  }

  tbl <- .resolve_daily_tbl(daily, location = location)
  required <- c("location", "sim_year", "date", "wind_kt", "event_id", "cum_damage", "damage_rate")
  missing_cols <- setdiff(required, names(tbl))
  if (length(missing_cols) > 0L) {
    stop(
      "`daily` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  tbl$date <- as.Date(tbl$date)

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
    return(tibble::tibble(
      sim_year = integer(0),
      focal_event_id = character(0),
      focal_start_date = as.Date(character(0)),
      focal_end_date = as.Date(character(0)),
      focal_peak_wind_kt = numeric(0),
      compound_window_end_date = as.Date(character(0)),
      compound_n_events = integer(0),
      compound_n_aftermath_events = integer(0),
      compound_cum_damage = numeric(0),
      compound_max_damage_rate = numeric(0)
    ))
  }

  year_ids <- sort(unique(tbl$sim_year))
  result <- vector("list", length(year_ids))

  for (i in seq_along(year_ids)) {
    yr <- year_ids[i]
    rows <- tbl[tbl$sim_year == yr, , drop = FALSE]
    rows <- rows[order(rows$date, rows$location), , drop = FALSE]

    event_rows <- rows[!is.na(rows$event_id), , drop = FALSE]
    if (nrow(event_rows) == 0L) {
      peak_idx <- which.max(rows$wind_kt)
      focal_start <- rows$date[peak_idx]
      focal_end <- rows$date[peak_idx]
      focal_event_id <- NA_character_
      focal_peak <- rows$wind_kt[peak_idx]
    } else {
      peak_row <- event_rows[which.max(event_rows$wind_kt), , drop = FALSE]
      focal_event_id <- peak_row$event_id[[1L]]
      focal_block <- event_rows[event_rows$event_id == focal_event_id, , drop = FALSE]
      focal_start <- min(focal_block$date)
      focal_end <- max(focal_block$date)
      focal_peak <- max(focal_block$wind_kt, na.rm = TRUE)
    }

    window_end <- focal_end + window_days
    window_rows <- rows[rows$date >= focal_start & rows$date <= window_end, , drop = FALSE]

    compound_ids <- unique(window_rows$event_id[!is.na(window_rows$event_id)])
    compound_n_events <- length(compound_ids)
    compound_n_aftermath_events <- sum(compound_ids != focal_event_id)

    damage_delta <- 0
    for (loc in unique(rows$location)) {
      loc_rows <- rows[rows$location == loc, , drop = FALSE]
      damage_before <- loc_rows$cum_damage[loc_rows$date < focal_start]
      damage_before <- if (length(damage_before) > 0L) {
        max(damage_before, na.rm = TRUE)
      } else {
        0
      }
      damage_after <- loc_rows$cum_damage[loc_rows$date <= window_end]
      damage_after <- if (length(damage_after) > 0L) {
        max(damage_after, na.rm = TRUE)
      } else {
        0
      }
      damage_before <- if (is.finite(damage_before)) damage_before else 0
      damage_after <- if (is.finite(damage_after)) damage_after else 0
      damage_delta <- damage_delta + max(0, damage_after - damage_before)
    }

    compound_max_damage_rate <- if (nrow(window_rows) > 0L) {
      max(window_rows$damage_rate, na.rm = TRUE)
    } else {
      0
    }
    compound_max_damage_rate <- if (is.finite(compound_max_damage_rate)) {
      compound_max_damage_rate
    } else {
      0
    }

    result[[i]] <- tibble::tibble(
      sim_year = as.integer(yr),
      focal_event_id = focal_event_id,
      focal_start_date = focal_start,
      focal_end_date = focal_end,
      focal_peak_wind_kt = focal_peak,
      compound_window_end_date = window_end,
      compound_n_events = as.integer(compound_n_events),
      compound_n_aftermath_events = as.integer(compound_n_aftermath_events),
      compound_cum_damage = damage_delta,
      compound_max_damage_rate = compound_max_damage_rate
    )
  }

  dplyr::bind_rows(result)
}


# =============================================================================
# 1) compute_stress_year_metrics()
# =============================================================================

#' Compute per-year multi-metric summary for stress-test characterisation
#'
#' @description
#' Aggregates the daily hazard-impact series to one row per
#' \code{(location, sim_year)}, capturing seven complementary properties of
#' each simulated year:
#'
#' \describe{
#'   \item{\code{peak_wind_kt}}{Maximum daily sustained wind (kt).
#'     Primary intensity indicator.}
#'   \item{\code{n_ts_days}}{Days with sustained wind \eqn{\ge} 34 kt.
#'     Captures tropical-storm-level duration exposure.}
#'   \item{\code{n_hur_days}}{Days with sustained wind \eqn{\ge} 64 kt.
#'     Captures hurricane-level duration exposure.}
#'   \item{\code{n_events}}{Count of distinct storm events in the year.
#'     Proxy for compound-event risk; years with multiple events have
#'     compounding recovery and damage burdens.}
#'   \item{\code{max_event_dur_days}}{Duration in days of the longest single
#'     event (number of days attributed to one \code{event_id}). Identifies
#'     slow-moving or stalling storms.}
#'   \item{\code{cum_damage}}{Total cumulative damage fraction over the year
#'     (\code{max(cum_damage)} across all days). Aggregate impact measure.}
#'   \item{\code{max_damage_rate}}{Maximum single-day damage rate.
#'     Captures acute shock intensity independently of duration.}
#' }
#'
#' @param daily Named list of tibbles returned by
#'   \code{\link{generate_daily_hazard_impact_spatial}()}, or a single tibble.
#' @param sim_years Optional filter. Either:
#'   \itemize{
#'     \item An integer vector of \code{sim_year} values applied to all
#'           locations.
#'     \item A tibble with \code{location} and \code{sim_year} columns (the
#'           direct output of \code{\link{query_storm_track_years}()} or
#'           \code{\link{query_impact_years}()}) for per-location filtering.
#'     \item \code{NULL} (default) to include all years.
#'   }
#' @param location Character vector of location names to include, or
#'   \code{NULL} for all.
#'
#' @return Tibble with columns \code{location}, \code{sim_year}, and the
#'   seven metrics listed above. Metric values are zero for calm years
#'   (no TC activity).
#'
#' @examples
#' \dontrun{
#' # Compute metrics for all impact-based candidate years
#' filtered   <- query_impact_years(daily_out, storm_id = "AL112017", out = hazard_out)
#' metrics    <- compute_stress_year_metrics(daily_out, sim_years = filtered)
#' }
#' @seealso \code{\link{aggregate_stress_metrics}}, \code{\link{select_stress_years}},
#'   \code{\link{query_impact_years}}, \code{\link{query_storm_track_years}}
#' @keywords internal
#' @export
compute_stress_year_metrics <- function(daily, sim_years = NULL, location = NULL) {

  tbl <- .resolve_daily_tbl(daily, location)

  # --- Apply sim_years filter (supports vector or query-output tibble) -------
  locs <- unique(tbl$location)
  yr_filter <- .resolve_sim_years_filter(sim_years, locs)

  if (!is.null(yr_filter)) {
    tbl <- dplyr::bind_rows(lapply(locs, function(loc) {
      yrs <- yr_filter[[loc]]
      if (length(yrs) == 0L) return(tbl[0L, ])
      dplyr::filter(tbl, .data$location == loc, .data$sim_year %in% yrs)
    }))
  }

  if (nrow(tbl) == 0L) {
    warning("No rows remain after applying sim_years filter.", call. = FALSE)
    return(tibble::tibble(
      location = character(0), sim_year = integer(0),
      peak_wind_kt = numeric(0), n_ts_days = integer(0),
      n_hur_days = integer(0), n_events = integer(0),
      max_event_dur_days = integer(0),
      cum_damage = numeric(0), max_damage_rate = numeric(0)
    ))
  }

  .summarise_daily_year_metrics(tbl)
}


# =============================================================================
# 2) aggregate_stress_metrics()
# =============================================================================

#' Compute a weighted composite stress score from per-year metrics
#'
#' @description
#' Aggregates metrics across locations (portfolio mode), normalises the chosen
#' metric columns to \eqn{[0, 1]} using min-max scaling, then combines them
#' into a single \code{composite_score} via a weighted mean.
#'
#' \strong{Portfolio mode (default, multiple locations):}
#' When the input contains more than one location, metric values are first
#' aggregated across locations for each \code{sim_year} using
#' \code{location_agg} (\code{"mean"}, \code{"max"}, or \code{"sum"}).  The
#' result has one row per \code{sim_year} with no \code{location} column.
#' This is the recommended path when a single representative year-set is
#' needed for the full location portfolio.
#'
#' \strong{Single-location mode:}
#' When only one location is present the input rows are scored directly and
#' the \code{location} column is retained.
#'
#' Normalisation is relative to the rows supplied, not the full ensemble.
#' Pass the full-ensemble metrics table if you need ensemble-relative scores.
#'
#' @param metrics Tibble returned by \code{\link{compute_stress_year_metrics}()}.
#' @param metrics_used Character vector of metric column names to include in
#'   the composite.  Defaults to all numeric columns except \code{sim_year}
#'   and any previously computed score/rank columns.
#' @param weights Optional numeric vector specifying relative importance.
#'   \itemize{
#'     \item \strong{Named vector} — names must match \code{metrics_used};
#'           missing names receive weight 0.
#'     \item \strong{Unnamed vector} — must have the same length as
#'           \code{metrics_used}; assigned in order.
#'     \item \code{NULL} (default) — uniform weights across all selected
#'           metrics.
#'   }
#'   Weights are automatically rescaled to sum to 1.
#' @param location_agg Character scalar controlling how metrics are aggregated
#'   across locations before scoring.  One of \code{"mean"} (default),
#'   \code{"max"}, or \code{"sum"}.  Ignored when only one location is present.
#'
#' @return Tibble with two additional columns:
#' \describe{
#'   \item{\code{composite_score}}{Weighted mean of min-max-normalised metric
#'     values; 0 = least extreme, 1 = most extreme within the provided set.}
#'   \item{\code{composite_rank}}{Integer rank (1 = highest score).}
#' }
#' When multiple locations are present (portfolio mode) the returned tibble has
#' columns \code{sim_year}, the metric columns, \code{composite_score}, and
#' \code{composite_rank} — no \code{location} column.
#'
#' @examples
#' \dontrun{
#' # Portfolio mode — single year-set across all locations (default)
#' scored <- aggregate_stress_metrics(metrics)
#'
#' # Up-weight intensity and damage, use max across locations
#' scored <- aggregate_stress_metrics(
#'   metrics,
#'   metrics_used = c("peak_wind_kt", "cum_damage", "n_events"),
#'   weights      = c(peak_wind_kt = 2, cum_damage = 2, n_events = 1),
#'   location_agg = "max"
#' )
#' }
#' @seealso \code{\link{compute_stress_year_metrics}}, \code{\link{select_stress_years}}
#' @keywords internal
#' @export
aggregate_stress_metrics <- function(
    metrics,
    metrics_used  = NULL,
    weights       = NULL,
    location_agg  = c("mean", "max", "sum")) {

  location_agg <- match.arg(location_agg)

  if (!is.data.frame(metrics)) {
    stop("`metrics` must be a data frame.", call. = FALSE)
  }

  id_cols     <- c("location", "sim_year")
  skip_cols   <- c(id_cols, "composite_score", "composite_rank")
  all_numeric <- names(metrics)[
    vapply(metrics, is.numeric, logical(1)) & !(names(metrics) %in% skip_cols)
  ]

  if (is.null(metrics_used)) {
    metrics_used <- all_numeric
  } else {
    bad <- setdiff(metrics_used, all_numeric)
    if (length(bad) > 0L) {
      stop(
        "`metrics_used` refers to columns that are not numeric metrics: ",
        paste(bad, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  if (length(metrics_used) == 0L) {
    stop("No numeric metric columns found to aggregate.", call. = FALSE)
  }

  # --- Resolve weights -------------------------------------------------------
  w <- .resolve_weights(weights, metrics_used)

  # --- Cross-location aggregation (portfolio mode) ---------------------------
  has_location <- "location" %in% names(metrics)
  n_locs       <- if (has_location) dplyr::n_distinct(metrics$location) else 1L

  if (has_location && n_locs > 1L) {
    agg_fn <- switch(location_agg, mean = mean, max = max, sum = sum)
    metrics <- metrics |>
      dplyr::group_by(.data$sim_year) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(metrics_used), ~ agg_fn(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    has_location <- FALSE
  }

  # --- Min-max normalise -----------------------------------------------------
  out      <- metrics
  norm_mat <- matrix(NA_real_, nrow = nrow(out), ncol = length(metrics_used),
                     dimnames = list(NULL, metrics_used))

  if (has_location) {
    # Single-location: normalise within that location's rows
    for (loc in unique(out$location)) {
      idx <- out$location == loc
      for (m in metrics_used) {
        x   <- out[[m]][idx]
        rng <- range(x, na.rm = TRUE)
        norm_mat[idx, m] <- if (is.finite(rng[1]) && rng[2] > rng[1]) {
          (x - rng[1]) / (rng[2] - rng[1])
        } else {
          rep(0.5, sum(idx))
        }
      }
    }
  } else {
    # Portfolio or single-location after aggregation: normalise globally
    for (m in metrics_used) {
      x   <- out[[m]]
      rng <- range(x, na.rm = TRUE)
      norm_mat[, m] <- if (is.finite(rng[1]) && rng[2] > rng[1]) {
        (x - rng[1]) / (rng[2] - rng[1])
      } else {
        rep(0.5, nrow(out))
      }
    }
  }

  # --- Weighted composite score ----------------------------------------------
  out$composite_score <- as.numeric(norm_mat %*% w[metrics_used])

  # --- Rank ------------------------------------------------------------------
  if (has_location) {
    out <- out |>
      dplyr::group_by(.data$location) |>
      dplyr::mutate(
        composite_rank = rank(-.data$composite_score, ties.method = "first")
      ) |>
      dplyr::ungroup()
  } else {
    out$composite_rank <- as.integer(
      rank(-out$composite_score, ties.method = "first")
    )
  }

  out
}


#' Resolve and validate a weights argument
#'
#' @param weights NULL, named numeric vector, or unnamed numeric vector.
#' @param metrics_used Character vector of metric names.
#' @return Named numeric vector, names = metrics_used, values summing to 1.
#' @keywords internal
.resolve_weights <- function(weights, metrics_used) {
  p <- length(metrics_used)

  if (is.null(weights)) {
    return(stats::setNames(rep(1 / p, p), metrics_used))
  }

  if (!is.numeric(weights) || any(weights < 0, na.rm = TRUE)) {
    stop("`weights` must be a non-negative numeric vector.", call. = FALSE)
  }

  if (!is.null(names(weights))) {
    unknown <- setdiff(names(weights), metrics_used)
    if (length(unknown) > 0L) {
      warning(
        "Ignoring weight(s) for unrecognised metric(s): ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    w <- stats::setNames(rep(0, p), metrics_used)
    common <- intersect(names(weights), metrics_used)
    w[common] <- weights[common]
  } else {
    if (length(weights) != p) {
      stop(
        "Unnamed `weights` must have length ", p,
        " (one per metric in `metrics_used`), got ", length(weights), ".",
        call. = FALSE
      )
    }
    w <- stats::setNames(as.numeric(weights), metrics_used)
  }

  if (sum(w) == 0) stop("At least one weight must be positive.", call. = FALSE)
  w / sum(w)
}


# =============================================================================
# 3) select_stress_years()
# =============================================================================

#' Select k representative stress-test years with good sample coverage
#'
#' @description
#' Picks a small subset of \code{k} years from the candidate set that
#' collectively represent the full distribution of stress severity.  Three
#' selection strategies are available:
#'
#' \describe{
#'   \item{\code{"stratified"} (default)}{Divides the \code{composite_score}
#'     range into \code{k} equal-count bins and returns the year whose score
#'     is closest to each bin's median.  Guarantees one representative per
#'     severity level — from near-miss to worst-case.}
#'   \item{\code{"diverse"}}{Greedy maximin selection in normalised metric
#'     space.  Starts from the most extreme year, then iteratively adds the
#'     year that is furthest from all already-selected years.  Maximises
#'     spread in the multi-dimensional metric space.}
#'   \item{\code{"top"}}{Returns the \code{k} years with the highest
#'     \code{composite_score}.  Useful when the goal is a set of purely
#'     high-severity scenarios.}
#' }
#'
#' \strong{Portfolio mode vs per-location mode:}
#' When \code{metrics} has no \code{location} column (the typical output of
#' \code{\link{aggregate_stress_metrics}()} with multiple locations), all rows
#' are treated as a single pool and \code{k} years are selected for the full
#' portfolio.  When a \code{location} column is present (single-location input
#' or a pre-aggregated table), selection is performed independently within
#' each location and the result is sorted by \code{location} then
#' \code{selection_rank}.
#'
#' If \code{composite_score} is absent from \code{metrics},
#' \code{\link{aggregate_stress_metrics}()} is called internally using the
#' supplied \code{weights} and \code{metrics_used}.
#'
#' @param metrics Tibble from \code{\link{compute_stress_year_metrics}()} or
#'   \code{\link{aggregate_stress_metrics}()}.
#' @param k Positive integer: number of years to select (portfolio-wide, or
#'   per location when a \code{location} column is present).
#'   If the candidate pool has \eqn{\le k} rows, all are returned.
#' @param method Character scalar selection strategy; one of
#'   \code{"stratified"} (default), \code{"diverse"}, or \code{"top"}.
#' @param metrics_used Forwarded to \code{\link{aggregate_stress_metrics}()}
#'   if \code{composite_score} is not already present.
#' @param weights Forwarded to \code{\link{aggregate_stress_metrics}()} if
#'   \code{composite_score} is not already present.
#'
#' @return Subset of \code{metrics} rows for the selected years with an
#'   additional integer column \code{selection_rank} (1 = first/most important
#'   pick, as determined by the method).  In portfolio mode rows are sorted by
#'   \code{selection_rank}; in per-location mode by \code{location} then
#'   \code{selection_rank}.
#'
#' @examples
#' \dontrun{
#' # Portfolio mode — k years for the full location set (typical use)
#' scored   <- aggregate_stress_metrics(metrics)       # no location column
#' selected <- select_stress_years(scored, k = 10)
#'
#' # Diversity-maximising selection with custom scoring
#' selected <- select_stress_years(
#'   scored,
#'   k            = 8,
#'   method       = "diverse",
#'   metrics_used = c("peak_wind_kt", "cum_damage", "n_events"),
#'   weights      = c(2, 1, 0.5)
#' )
#' }
#' @seealso \code{\link{compute_stress_year_metrics}},
#'   \code{\link{aggregate_stress_metrics}}
#' @keywords internal
#' @export
select_stress_years <- function(
    metrics,
    k,
    method       = c("stratified", "diverse", "top"),
    metrics_used = NULL,
    weights      = NULL) {

  method <- match.arg(method)

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  # Ensure composite_score exists
  if (!("composite_score" %in% names(metrics))) {
    metrics <- aggregate_stress_metrics(
      metrics,
      metrics_used = metrics_used,
      weights      = weights
    )
  }

  # Columns to use for the diverse-method distance computation
  id_cols   <- c("location", "sim_year", "composite_score", "composite_rank", "selection_rank")
  dist_cols <- if (!is.null(metrics_used)) {
    metrics_used
  } else {
    names(metrics)[vapply(metrics, is.numeric, logical(1)) & !(names(metrics) %in% id_cols)]
  }

  # ---- Portfolio mode: no location column — select from a single pool -------
  has_location <- "location" %in% names(metrics)

  if (!has_location) {
    if (nrow(metrics) <= k) {
      metrics$selection_rank <- as.integer(
        rank(-metrics$composite_score, ties.method = "first")
      )
      return(dplyr::arrange(metrics, .data$selection_rank))
    }
    result <- switch(
      method,
      stratified = .select_stratified(metrics, k),
      diverse    = .select_diverse(metrics, k, dist_cols),
      top        = .select_top(metrics, k)
    )
    return(dplyr::arrange(result, .data$selection_rank))
  }

  # ---- Per-location mode: location column present ---------------------------
  locs        <- unique(metrics$location)
  result_list <- vector("list", length(locs))

  for (i in seq_along(locs)) {
    loc    <- locs[i]
    loc_df <- metrics[metrics$location == loc, , drop = FALSE]

    if (nrow(loc_df) <= k) {
      loc_df$selection_rank <- rank(-loc_df$composite_score, ties.method = "first")
      result_list[[i]] <- loc_df
      next
    }

    result_list[[i]] <- switch(
      method,
      stratified = .select_stratified(loc_df, k),
      diverse    = .select_diverse(loc_df, k, dist_cols),
      top        = .select_top(loc_df, k)
    )
  }

  dplyr::bind_rows(result_list) |>
    dplyr::arrange(.data$location, .data$selection_rank)
}


# --- Selection helpers --------------------------------------------------------

#' @keywords internal
.select_top <- function(loc_df, k) {
  loc_df |>
    dplyr::arrange(dplyr::desc(.data$composite_score)) |>
    utils::head(k) |>
    dplyr::mutate(selection_rank = seq_len(k))
}

#' @keywords internal
.select_stratified <- function(loc_df, k) {
  n       <- nrow(loc_df)
  sorted  <- loc_df[order(loc_df$composite_score), , drop = FALSE]

  # Assign each row to one of k equal-count bins (1 = lowest score)
  sorted$bin <- ceiling(seq_len(n) * k / n)

  selected <- vector("list", k)
  for (b in seq_len(k)) {
    bin_rows <- sorted[sorted$bin == b, , drop = FALSE]
    med      <- stats::median(bin_rows$composite_score)
    # Pick the year within the bin whose score is closest to the bin median
    best     <- bin_rows[which.min(abs(bin_rows$composite_score - med)), , drop = FALSE]
    best$selection_rank <- b   # rank 1 = lowest severity bin, k = highest
    selected[[b]] <- best
  }

  dplyr::bind_rows(selected) |>
    dplyr::select(-"bin")
}

#' @keywords internal
.select_diverse <- function(loc_df, k, dist_cols) {
  # Normalise dist_cols to [0, 1] for Euclidean distance computation
  mat <- as.matrix(loc_df[, intersect(dist_cols, names(loc_df)), drop = FALSE])

  for (j in seq_len(ncol(mat))) {
    rng <- range(mat[, j], na.rm = TRUE)
    mat[, j] <- if (is.finite(rng[1]) && rng[2] > rng[1]) {
      (mat[, j] - rng[1]) / (rng[2] - rng[1])
    } else {
      rep(0.5, nrow(mat))
    }
  }

  n            <- nrow(loc_df)
  selected_idx <- integer(k)

  # Seed: start with the most extreme year (highest composite_score)
  selected_idx[1L] <- which.max(loc_df$composite_score)

  for (i in seq(2L, k)) {
    # For each candidate compute its minimum distance to any selected year
    min_dist <- numeric(n)
    for (j in seq_len(n)) {
      dists <- vapply(
        selected_idx[seq_len(i - 1L)],
        function(s) sqrt(sum((mat[j, ] - mat[s, ])^2, na.rm = TRUE)),
        numeric(1)
      )
      min_dist[j] <- min(dists)
    }
    # Exclude already-selected years
    min_dist[selected_idx[seq_len(i - 1L)]] <- -Inf
    selected_idx[i] <- which.max(min_dist)
  }

  loc_df[selected_idx, , drop = FALSE] |>
    dplyr::mutate(selection_rank = seq_len(k))
}
