# =============================================================================
# Visualization functions for hurricane hazard model
# =============================================================================

# -----------------------------------------------------------------------------
# Theme and palette
# -----------------------------------------------------------------------------

plot_theme <- function(base_size = 11) {

  theme_light(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(color = "grey20", face = "bold"),
      legend.position = "bottom"
    )
}

plot_colors <- list(
  event = c(TS = "#E69F00", HUR = "#d73027"),
  wind = c(mean = "grey40", max = "#a50026"),
  threshold = c(TS = "#E69F00", HUR = "#d73027"),
  quantile = c(median = "grey30", p95 = "#91bfdb", p99 = "#4575b4")
)

# -----------------------------------------------------------------------------
# Data preparation helpers
# -----------------------------------------------------------------------------

#' Add day/month/year fields to daily hazard data
#'
#' @param daily_impact A data frame or tibble with at least a `date` column coercible
#'   by `format()` (typically `Date`), one row per day.
#'
#' @return `daily_impact` with three added integer columns: `doy` (1-366), `month`
#'   (1-12), and `year` (4-digit calendar year).
#' @keywords internal
prep_daily <- function(daily_impact) {

  daily_impact %>%
    mutate(

      doy   = as.integer(format(date, "%j")),
      month = as.integer(format(date, "%m")),
      year  = as.integer(format(date, "%Y"))
    )
}

#' Summarize event-level features from daily rows
#'
#' @param daily A data frame or tibble containing daily rows with columns
#'   `location`, `sim_year`, `event_id`, `event_class`, `date`, and `wind_kt`.
#'   Rows with `NA` in `event_id` are excluded.
#'
#' @return A tibble with one row per unique (`location`, `sim_year`, `event_id`)
#'   and columns `event_class`, `start_date`, `end_date`, `dur_days`,
#'   `max_wind_kt`, `start_doy`, and `start_month`.
#' @keywords internal
prep_events <- function(daily) {
  daily %>%
    filter(!is.na(event_id)) %>%
    group_by(location, sim_year, event_id) %>%
    summarise(
      event_class = first(na.omit(event_class)),
      start_date  = min(date),
      end_date    = max(date),
      dur_days    = as.integer(end_date - start_date) + 1L,
      max_wind_kt = max(wind_kt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      event_class = ifelse(event_class == "HUR", "HUR", "TS"),
      start_doy   = as.integer(format(start_date, "%j")),
      start_month = as.integer(format(start_date, "%m"))
    )
}

# -----------------------------------------------------------------------------
# Plot 1: Wind time series with event overlay
# -----------------------------------------------------------------------------

#' Plot daily wind speed with event-duration overlays
#'
#' @description Draws a daily wind time series and overlays each event as a
#'   horizontal segment spanning event start-to-end dates at that event's peak
#'   wind speed.
#'
#' @details If `events` is `NULL`, event summaries are computed from `daily`.
#'   Segment colors collapse event classes to `TS` and `HUR` (anything not equal
#'   to `"HUR"` is treated as `TS`). `NA` values in `event_id` are ignored when
#'   computing events.
#'
#' @param daily A data frame or tibble of daily records with columns `date`
#'   (`Date`), `wind_kt` (numeric, knots), `event_id` (event identifier; `NA`
#'   means no event), and `event_class` (character/factor).
#' @param events Optional precomputed event table (data frame/tibble) with
#'   columns `start_date`, `end_date` (`Date`), `max_wind_kt` (numeric, knots),
#'   and `event_class` (`"TS"`/`"HUR"`). If `NULL`, it is derived from `daily`.
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots.
#' @param show_thresholds Logical scalar; if `TRUE`, dashed horizontal threshold
#'   lines are added.
#' @param title Character scalar used as plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-08-01") + 0:9,
#'   wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 12, 8),
#'   event_id = c(NA, NA, 1, 1, 1, NA, NA, NA, NA, NA),
#'   event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, NA, NA, NA),
#'   location = "A",
#'   sim_year = 2001
#' )
#' plot_wind_timeseries(d)
#' @export
plot_wind_timeseries <- function(daily,
                                 events = NULL,
                                 thr_tc = 34,
                                 thr_hur = 64,
                                 show_thresholds = TRUE,
                                 title = "Daily Wind with TS/Hurricane Events") {

  if (is.null(events)) events <- prep_events(daily)


  p <- ggplot(daily, aes(x = date, y = wind_kt)) +
    geom_line(color = "grey50", linewidth = 0.4) +
    geom_segment(
      data = events,
      aes(x = start_date, xend = end_date,
          y = max_wind_kt, yend = max_wind_kt,
          color = event_class),
      linewidth = 1.5, lineend = "round",
      inherit.aes = FALSE, alpha = 0.85
    ) +
    scale_color_manual(values = plot_colors$event, name = "Event") +
    scale_x_date(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Wind speed (kt)", title = title,
         subtitle = "Segments show event duration at peak intensity") +
    plot_theme()

  if (show_thresholds) {
    p <- p +
      geom_hline(yintercept = thr_tc,  linetype = "dashed",
                 color = plot_colors$threshold["TS"], alpha = 0.7) +
      geom_hline(yintercept = thr_hur, linetype = "dashed",
                 color = plot_colors$threshold["HUR"], alpha = 0.7)
  }

  p
}

# -----------------------------------------------------------------------------
# Plot 2: Seasonal distribution (DOY histogram)
# -----------------------------------------------------------------------------

#' Plot day-of-year distribution of event timing
#'
#' @description Plots a day-of-year histogram for either event-active days or
#'   event starts, optionally faceted by event class.
#'
#' @details For `metric = "event_days"`, each day with non-`NA` `event_id`
#'   contributes one count. For `metric = "starts"`, counts are based on event
#'   start dates derived from grouped events. Event classes are collapsed to
#'   `TS`/`HUR` using the same rule as `prep_events()`.
#'
#' @param daily A data frame or tibble with at least `date` (`Date`),
#'   `event_id`, `event_class`, plus fields required by `prep_events()` when
#'   `metric = "starts"` (`location`, `sim_year`, `wind_kt`).
#' @param metric Character scalar; one of `"event_days"` or `"starts"`.
#' @param facet_class Logical scalar; if `TRUE`, create one panel per class.
#' @param binwidth Numeric scalar > 0; histogram bin width in day-of-year units
#'   (days).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-08-01") + 0:9,
#'   wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 12, 8),
#'   event_id = c(NA, NA, 1, 1, 1, NA, NA, 2, 2, NA),
#'   event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, "HUR", "HUR", NA),
#'   location = "A",
#'   sim_year = 2001
#' )
#' plot_seasonality_doy(d, metric = "event_days")
#' @export
plot_seasonality_doy <- function(daily,
                                 metric = c("event_days", "starts"),
                                 facet_class = TRUE,
                                 binwidth = 7) {


  metric <- match.arg(metric)

  if (metric == "event_days") {
    # Duration-weighted: count each day an event is active
    plot_data <- daily %>%
      filter(!is.na(event_id), !is.na(event_class)) %>%
      mutate(
        doy = as.integer(format(date, "%j")),
        event_class = ifelse(event_class == "HUR", "HUR", "TS")
      )
    ylab <- "Event-days"
    subtitle <- "Duration-weighted: each day of event exposure counted"

  } else {
    # Starts only: count event initiations
    events <- prep_events(daily)
    plot_data <- events %>%
      transmute(doy = start_doy, event_class)
    ylab <- "Event starts"
    subtitle <- "Count of event initiations by day of year"
  }

  p <- ggplot(plot_data, aes(x = doy, fill = event_class)) +
    geom_histogram(position = "dodge", binwidth = binwidth, boundary = 0, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = plot_colors$event, name = "Class") +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ylab,
         title = "Seasonality of TS/HUR Activity",
         subtitle = subtitle) +
    plot_theme()

  if (facet_class) {
    p <- p + facet_wrap(~ event_class, ncol = 1, scales = "free_y")
  }

  p
}

# -----------------------------------------------------------------------------
# Plot 3: Monthly event counts (bar chart)
# -----------------------------------------------------------------------------

#' Plot monthly counts of event starts
#'
#' @description Summarizes event starts by calendar month and plots grouped bars
#'   by event class, optionally normalized to an annual rate.
#'
#' @details Event starts are computed from `prep_events()`. Missing month/class
#'   combinations are filled with zero. When `normalize = TRUE`, counts are
#'   divided by `n_distinct(daily$sim_year)`.
#'
#' @param daily A data frame or tibble with columns required by `prep_events()`
#'   (`location`, `sim_year`, `event_id`, `event_class`, `date`, `wind_kt`).
#'   Rows with `NA` `event_id` are excluded during event extraction.
#' @param normalize Logical scalar; if `TRUE`, plot events per year.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-08-01") + 0:9,
#'   wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 12, 8),
#'   event_id = c(NA, NA, 1, 1, 1, NA, NA, 2, 2, NA),
#'   event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, "HUR", "HUR", NA),
#'   location = "A",
#'   sim_year = 2001
#' )
#' plot_monthly_events(d)
#' @export
plot_monthly_events <- function(daily, normalize = FALSE) {

  events <- prep_events(daily)
  n_years <- n_distinct(daily$sim_year)

  # Count by month, fill zeros
  monthly <- events %>%
    count(event_class, start_month, name = "n") %>%
    complete(event_class = c("TS", "HUR"), start_month = 1:12, fill = list(n = 0))

  if (normalize) {
    monthly <- monthly %>% mutate(n = n / n_years)
    ylab <- "Events per year"
  } else {
    ylab <- "Total event starts"
  }

  ggplot(monthly, aes(x = factor(start_month), y = n, fill = event_class)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_x_discrete(labels = month.abb) +
    scale_fill_manual(values = plot_colors$event, name = "Class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ylab, title = "Monthly Distribution of Event Starts") +
    plot_theme()
}

# -----------------------------------------------------------------------------
# Plot 4: DOY wind summary (mean + max)
# -----------------------------------------------------------------------------

#' Plot mean and maximum wind by day of year
#'
#' @description Aggregates daily wind by day-of-year and plots mean and maximum
#'   wind speed curves with optional loess smoothing.
#'
#' @details Summary statistics use `na.rm = TRUE`; all-`NA` day groups propagate
#'   non-finite summaries. Threshold lines are always drawn at `thr_tc` and
#'   `thr_hur`.
#'
#' @param daily A data frame or tibble with columns `date` (`Date`) and
#'   `wind_kt` (numeric, knots).
#' @param smooth Logical scalar; if `TRUE`, use `geom_smooth(method = "loess")`,
#'   otherwise plot unsmoothed lines.
#' @param span Numeric scalar in `(0, 1]`; loess span used when `smooth = TRUE`.
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-01-01") + 0:29,
#'   wind_kt = 15 + sin((1:30) / 5) * 10
#' )
#' plot_doy_wind(d, smooth = FALSE)
#' @export
plot_doy_wind <- function(daily,
                          smooth = TRUE,
                          span = 0.15,
                          thr_tc = 34,
                          thr_hur = 64) {

  doy_summary <- daily %>%
    mutate(doy = as.integer(format(date, "%j"))) %>%
    group_by(doy) %>%
    summarise(
      mean_wind = mean(wind_kt, na.rm = TRUE),
      max_wind  = max(wind_kt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean_wind, max_wind),
                 names_to = "stat", values_to = "wind_kt") %>%
    mutate(stat = factor(stat, levels = c("mean_wind", "max_wind"),
                         labels = c("Mean", "Maximum")))

  p <- ggplot(doy_summary, aes(x = doy, y = wind_kt, color = stat))

  if (smooth) {
    p <- p + geom_smooth(method = "loess", span = span, se = FALSE, linewidth = 1)
  } else {
    p <- p + geom_line(linewidth = 0.6)
  }

  p +
    geom_hline(yintercept = thr_tc,  linetype = "dashed",
               color = plot_colors$threshold["TS"], alpha = 0.7) +
    geom_hline(yintercept = thr_hur, linetype = "dashed",
               color = plot_colors$threshold["HUR"], alpha = 0.7) +
    scale_color_manual(values = c(Mean = "grey40", Maximum = "steelblue"), name = NULL) +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Wind speed (kt)",
         title = "Daily Wind Speed by Day of Year") +
    plot_theme()
}

# -----------------------------------------------------------------------------
# Plot 5: Monthly wind quantiles
# -----------------------------------------------------------------------------

#' Plot monthly wind quantile summary
#'
#' @description Computes monthly wind quantiles and plots a median-to-upper
#'   ribbon with upper and extreme quantile lines.
#'
#' @details `probs` is indexed as lower/median display (`probs[1]`), ribbon
#'   upper bound (`probs[2]`), and extreme line (`probs[3]`). If `log_scale` is
#'   `TRUE`, quantile values are floored at 0.5 before `scale_y_log10()` to
#'   avoid non-positive values.
#'
#' @param daily A data frame or tibble with columns `date` (`Date`) and
#'   `wind_kt` (numeric, knots; `NA` allowed and ignored in quantiles).
#' @param probs Numeric vector of length 3 with probabilities in `[0, 1]`, used
#'   in
#'   order `probs[1]`, `probs[2]`, `probs[3]`.
#' @param log_scale Logical scalar; if `TRUE`, use a log10 y-axis.
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-01-01") + 0:59,
#'   wind_kt = pmax(1, 20 + rnorm(60, sd = 5))
#' )
#' plot_monthly_quantiles(d)
#' @export
plot_monthly_quantiles <- function(daily,
                                   probs = c(0.50, 0.95, 0.99),
                                   log_scale = FALSE,
                                   thr_tc = 34,
                                   thr_hur = 64) {

  monthly <- daily %>%
    mutate(month = as.integer(format(date, "%m"))) %>%
    group_by(month) %>%
    summarise(
      median = quantile(wind_kt, probs[1], na.rm = TRUE),
      upper  = quantile(wind_kt, probs[2], na.rm = TRUE),
      extreme = quantile(wind_kt, probs[3], na.rm = TRUE),
      .groups = "drop"
    )

  # For log scale, floor small values
  if (log_scale) {
    monthly <- monthly %>%
      mutate(across(c(median, upper, extreme), ~ pmax(.x, 0.5)))
  }

  p <- ggplot(monthly, aes(x = month)) +
    geom_ribbon(aes(ymin = median, ymax = upper), fill = "grey80", alpha = 0.8) +
    geom_line(aes(y = median), linewidth = 1, color = plot_colors$quantile["median"]) +
    geom_line(aes(y = upper), linewidth = 0.7, linetype = "dashed",
              color = plot_colors$quantile["p95"]) +
    geom_line(aes(y = extreme), linewidth = 0.7, linetype = "dotted",
              color = plot_colors$quantile["p99"]) +
    geom_hline(yintercept = thr_tc,  linetype = "dashed",
               color = plot_colors$threshold["TS"], alpha = 0.6) +
    geom_hline(yintercept = thr_hur, linetype = "dashed",
               color = plot_colors$threshold["HUR"], alpha = 0.6) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(
      x = NULL,
      y = "Wind speed (kt)",
      title = "Monthly Wind Distribution",
      subtitle = sprintf("Ribbon: median-P%d | Dotted: P%d",
                         probs[2]*100, probs[3]*100)
    ) +
    plot_theme()

  if (log_scale) p <- p + scale_y_log10()

  p
}

# -----------------------------------------------------------------------------
# Plot 6: Annual event count distribution
# -----------------------------------------------------------------------------
#' Plot distribution of annual event totals
#'
#' @description Computes annual totals from daily event flags and plots a
#'   histogram of either distinct events or event-days.
#'
#' @details Annual `n_events` counts distinct non-`NA` `event_id` values per
#'   `sim_year`; annual `n_days` counts daily rows with non-`NA` `event_id`. A
#'   Poisson expected-count overlay is added only when `metric = "events"` and
#'   `show_poisson = TRUE`.
#'
#' @param daily A data frame or tibble with columns `sim_year` and `event_id`
#'   (`NA` means no event that day).
#' @param metric Character scalar; one of `"events"` or `"days"`.
#' @param show_poisson Logical scalar; if `TRUE`, add Poisson expected points for
#'   the `"events"` metric.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   sim_year = rep(2001:2003, each = 5),
#'   event_id = c(NA, 1, 1, NA, 2, NA, NA, 3, 3, 3, NA, 4, NA, 5, NA)
#' )
#' plot_annual_counts(d, metric = "events")
#' @export
plot_annual_counts <- function(daily,
                               metric = c("events", "days"),
                               show_poisson = TRUE) {

  metric <- match.arg(metric)

  yearly <- daily %>%
    group_by(sim_year) %>%
    summarise(
      n_events = n_distinct(event_id[!is.na(event_id)]),
      n_days   = sum(!is.na(event_id)),
      .groups = "drop"
    )

  if (metric == "events") {
    yearly <- yearly %>% mutate(count = n_events)
    xlab <- "Events per year"
    title <- "Distribution of Annual Event Counts"
  } else {
    yearly <- yearly %>% mutate(count = n_days)
    xlab <- "Event-days per year"
    title <- "Distribution of Annual Event-Days"
  }

  lambda <- mean(yearly$count)
  n_years <- nrow(yearly)

  p <- ggplot(yearly, aes(x = count)) +
    geom_bar(fill = "steelblue", color = "white", width = 0.8) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = xlab, y = "Number of years", title = title,
         subtitle = sprintf("Mean: %.2f | n = %d years", lambda, n_years)) +
    plot_theme()

  if (show_poisson && metric == "events") {
    # Overlay Poisson PMF
    x_range <- 0:max(yearly$count + 2)
    pois_df <- data.frame(
      count = x_range,
      expected = dpois(x_range, lambda) * n_years
    )
    p <- p + geom_point(data = pois_df, aes(x = count, y = expected),
                        color = "red", size = 2, shape = 1)
  }

  p
}

# -----------------------------------------------------------------------------
# Plot 7: Event intensity vs duration scatter
# -----------------------------------------------------------------------------

#' Plot event intensity versus duration
#'
#' @description Creates a scatter plot of event duration (days) against event
#'   peak wind speed, colored by event class.
#'
#' @details If `events` is `NULL`, `daily` must be supplied and is summarized via
#'   `prep_events()`. The function errors if both are `NULL`.
#'
#' @param daily Optional data frame/tibble of daily rows used only when
#'   `events = NULL`; must contain columns required by `prep_events()`.
#' @param events Optional data frame/tibble with columns `dur_days` (integer
#'   days), `max_wind_kt` (numeric, knots), and `event_class` (`"TS"`/`"HUR"`).
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   date = as.Date("2001-08-01") + 0:9,
#'   wind_kt = c(10, 20, 40, 55, 30, 15, 5, 10, 70, 65),
#'   event_id = c(NA, NA, 1, 1, 1, NA, NA, NA, 2, 2),
#'   event_class = c(NA, NA, "TS", "TS", "TS", NA, NA, NA, "HUR", "HUR"),
#'   location = "A",
#'   sim_year = 2001
#' )
#' plot_intensity_duration(daily = d)
#' @export
plot_intensity_duration <- function(daily = NULL,
                                    events = NULL,
                                    thr_tc = 34,
                                    thr_hur = 64) {

  if (is.null(events)) {
    stopifnot(!is.null(daily))
    events <- prep_events(daily)
  }

  ggplot(events, aes(x = dur_days, y = max_wind_kt, color = event_class)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = thr_tc,  linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = thr_hur, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = plot_colors$event, name = "Class") +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    labs(x = "Event duration (days)", y = "Peak wind (kt)",
         title = "Event Intensity vs Duration") +
    plot_theme()
}

# -----------------------------------------------------------------------------
# Plot 8: Wind histogram / density
# -----------------------------------------------------------------------------

#' Plot distribution of daily wind speed
#'
#' @description Plots daily wind speed as either a histogram or kernel density
#'   estimate with TS/HUR threshold reference lines.
#'
#' @details `type` is matched with `match.arg()`. If `log_scale = TRUE`,
#'   `scale_x_log10()` is applied; non-positive `wind_kt` values are therefore
#'   not displayable on the transformed axis.
#'
#' @param daily A data frame or tibble with column `wind_kt` (numeric, knots;
#'   `NA` values are ignored by ggplot stat layers).
#' @param type Character scalar; one of `"histogram"` or `"density"`.
#' @param log_scale Logical scalar; if `TRUE`, use a log10 x-axis.
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(wind_kt = c(5, 10, 15, 20, 35, 40, 65, 70))
#' plot_wind_distribution(d, type = "histogram")
#' @export
plot_wind_distribution <- function(daily,
                                   type = c("histogram", "density"),
                                   log_scale = FALSE,
                                   thr_tc = 34,
                                   thr_hur = 64) {

  type <- match.arg(type)

  p <- ggplot(daily, aes(x = wind_kt))

  if (type == "histogram") {
    p <- p + geom_histogram(binwidth = 5, fill = "steelblue",
                            color = "white", boundary = 0)
  } else {
    p <- p + geom_density(fill = "steelblue", alpha = 0.6, color = NA)
  }

  p <- p +
    geom_vline(xintercept = thr_tc,  linetype = "dashed",
               color = plot_colors$threshold["TS"]) +
    geom_vline(xintercept = thr_hur, linetype = "dashed",
               color = plot_colors$threshold["HUR"]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Wind speed (kt)", y = if(type == "histogram") "Count" else "Density",
         title = "Wind Speed Distribution") +
    plot_theme()

  if (log_scale) p <- p + scale_x_log10()

  p
}

# -----------------------------------------------------------------------------
# Plot 9: Return level plot
# -----------------------------------------------------------------------------

#' Plot empirical wind return levels
#'
#' @description Computes empirical return-period points using either annual block
#'   maxima or peaks-over-threshold exceedances and plots them on a log-scaled
#'   return-period axis.
#'
#' @details For `block_maxima = TRUE`, one annual maximum is taken per
#'   `sim_year`. For `block_maxima = FALSE`, `threshold` must be non-`NULL`, and
#'   only rows with `wind_kt > threshold` are used. Wind summaries use
#'   `na.rm = TRUE`.
#'
#' @param daily A data frame or tibble with columns `sim_year` and `wind_kt`
#'   (numeric, knots).
#' @param block_maxima Logical scalar; if `TRUE`, use annual maxima, otherwise
#'   use exceedances over `threshold`.
#' @param threshold Optional numeric scalar (knots); required when
#'   `block_maxima = FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   sim_year = rep(2001:2005, each = 3),
#'   wind_kt = c(20, 35, 50, 15, 30, 70, 10, 25, 45, 18, 40, 55, 22, 60, 65)
#' )
#' plot_return_levels(d, block_maxima = TRUE)
#' @export
plot_return_levels <- function(daily,
                               block_maxima = TRUE,
                               threshold = NULL) {

  if (block_maxima) {
    # Annual maxima approach
    maxima <- daily %>%
      group_by(sim_year) %>%
      summarise(max_wind = max(wind_kt, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_wind)) %>%
      mutate(
        rank = row_number(),
        n = n(),
        exceedance_prob = rank / (n + 1),
        return_period = 1 / exceedance_prob
      )

    subtitle <- "Block maxima (annual)"
  } else {
    # POT approach
    stopifnot(!is.null(threshold))
    n_years <- n_distinct(daily$sim_year)

    maxima <- daily %>%
      filter(wind_kt > threshold) %>%
      arrange(desc(wind_kt)) %>%
      mutate(
        rank = row_number(),
        n = n(),
        exceedance_prob = rank / (n + 1),
        annual_rate = n / n_years,
        return_period = 1 / (exceedance_prob * annual_rate / n_years)
      ) %>%
      rename(max_wind = wind_kt)

    subtitle <- sprintf("Peaks over threshold (%.0f kt)", threshold)
  }

  ggplot(maxima, aes(x = return_period, y = max_wind)) +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    scale_x_log10(
      breaks = c(1, 2, 5, 10, 25, 50, 100, 250, 500),
      labels = c("1", "2", "5", "10", "25", "50", "100", "250", "500")
    ) +
    labs(x = "Return period (years)", y = "Wind speed (kt)",
         title = "Empirical Return Levels",
         subtitle = subtitle) +
    plot_theme()
}

# -----------------------------------------------------------------------------
# Quick summary panel
# -----------------------------------------------------------------------------

#' Build a four-panel hurricane hazard summary figure
#'
#' @description Combines seasonality, monthly quantiles, annual event counts,
#'   and intensity-duration plots into a single 2x2 patchwork layout.
#'
#' @details Requires the `patchwork` package at runtime (`requireNamespace()` is
#'   called). Input requirements combine those of the component plotting
#'   functions.
#'
#' @param daily A data frame or tibble with columns needed by
#'   `plot_seasonality_doy()`, `plot_monthly_quantiles()`,
#'   `plot_annual_counts()`, and `plot_intensity_duration()`.
#' @param thr_tc Numeric scalar; tropical-cyclone threshold in knots passed to
#'   component plots.
#' @param thr_hur Numeric scalar; hurricane threshold in knots passed to
#'   component plots.
#'
#' @return A `patchwork`/`ggplot` composed plot object.
#'
#' @examples
#' \dontrun{
#' d <- data.frame(
#'   date = as.Date("2001-01-01") + 0:364,
#'   wind_kt = pmax(1, 20 + rnorm(365, sd = 8)),
#'   event_id = sample(c(NA, 1:8), 365, replace = TRUE),
#'   event_class = sample(c(NA, "TS", "HUR"), 365, replace = TRUE),
#'   location = "A",
#'   sim_year = 2001
#' )
#' plot_summary_panel(d)
#' }
#' @export
plot_summary_panel <- function(daily, thr_tc = 34, thr_hur = 64) {

  requireNamespace("patchwork", quietly = TRUE)

  p1 <- plot_seasonality_doy(daily, metric = "event_days", facet_class = FALSE)
  p2 <- plot_monthly_quantiles(daily, thr_tc = thr_tc, thr_hur = thr_hur)
  p3 <- plot_annual_counts(daily, metric = "events")
  p4 <- plot_intensity_duration(daily, thr_tc = thr_tc, thr_hur = thr_hur)

  (p1 + p2) / (p3 + p4) +
    patchwork::plot_annotation(
      title = "Hurricane Hazard Summary",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}
