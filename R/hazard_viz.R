# =============================================================================
# Visualization functions for hurricane hazard model
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Theme and palette
# -----------------------------------------------------------------------------

passaat_theme <- function(base_size = 11) {
  
  theme_light(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(color = "grey20", face = "bold"),
      legend.position = "bottom"
    )
}

passaat_colors <- list(
  event = c(TC = "#E69F00", HUR = "#D55E00"),
  wind = c(mean = "grey40", max = "steelblue"),
  threshold = c(tc = "#E69F00", hur = "#D55E00"),
  quantile = c(median = "grey30", p95 = "steelblue", p99 = "#D55E00")
)

# -----------------------------------------------------------------------------
# Data preparation helpers
# -----------------------------------------------------------------------------

#' Prepare daily data with time keys
#' @param daily_impact Daily impact data frame (requires: date, wind_kt)
prep_daily <- function(daily_impact) {
  
  daily_impact %>%
    mutate(
      
      doy   = as.integer(format(date, "%j")),
      month = as.integer(format(date, "%m")),
      year  = as.integer(format(date, "%Y"))
    )
}

#' Extract event-level summary from daily data
#' @param daily Daily data with event_id column
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
      event_class = ifelse(event_class == "HUR", "HUR", "TC"),
      start_doy   = as.integer(format(start_date, "%j")),
      start_month = as.integer(format(start_date, "%m"))
    )
}

# -----------------------------------------------------------------------------
# Plot 1: Wind time series with event overlay
# -----------------------------------------------------------------------------

#' Plot daily wind time series with event segments
#'
#' @param daily Daily data frame (date, wind_kt, event_id, event_class)
#' @param events Optional pre-computed events table; computed if NULL
#' @param thr_tc Tropical cyclone threshold (kt)
#' @param thr_hur Hurricane threshold (kt)
#' @param show_thresholds Show threshold lines
#' @param title Plot title
plot_wind_timeseries <- function(daily,
                                 events = NULL,
                                 thr_tc = 34,
                                 thr_hur = 64,
                                 show_thresholds = TRUE,
                                 title = "Daily Wind with TC/Hurricane Events") {
  
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
    scale_color_manual(values = passaat_colors$event, name = "Event") +
    scale_x_date(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Wind speed (kt)", title = title,
         subtitle = "Segments show event duration at peak intensity") +
    passaat_theme()
  
  if (show_thresholds) {
    p <- p +
      geom_hline(yintercept = thr_tc,  linetype = "dashed", 
                 color = passaat_colors$threshold["tc"], alpha = 0.7) +
      geom_hline(yintercept = thr_hur, linetype = "dashed", 
                 color = passaat_colors$threshold["hur"], alpha = 0.7)
  }
  
  p
}

# -----------------------------------------------------------------------------
# Plot 2: Seasonal distribution (DOY histogram)
# -----------------------------------------------------------------------------

#' Plot day-of-year distribution of event activity
#'
#' @param daily Daily data with event_id, event_class
#' @param metric One of "event_days" (duration-weighted) or "starts" (initiation)
#' @param facet_class Facet by event class (TC/HUR)
#' @param binwidth DOY bin width in days
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
        event_class = ifelse(event_class == "HUR", "HUR", "TC")
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
    geom_histogram(binwidth = binwidth, boundary = 0, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = passaat_colors$event, name = "Class") +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ylab, 
         title = "Seasonality of TC Activity",
         subtitle = subtitle) +
    passaat_theme()
  
  if (facet_class) {
    p <- p + facet_wrap(~ event_class, ncol = 1, scales = "free_y")
  }
  
  p
}

# -----------------------------------------------------------------------------
# Plot 3: Monthly event counts (bar chart)
# -----------------------------------------------------------------------------

#' Plot monthly summary of event starts
#'
#' @param daily Daily data frame
#' @param normalize Divide by number of years to show annual rate
plot_monthly_events <- function(daily, normalize = FALSE) {
  
  events <- prep_events(daily)
  n_years <- n_distinct(daily$sim_year)
  
  # Count by month, fill zeros
  monthly <- events %>%
    count(event_class, start_month, name = "n") %>%
    complete(event_class = c("TC", "HUR"), start_month = 1:12, fill = list(n = 0))
  
  if (normalize) {
    monthly <- monthly %>% mutate(n = n / n_years)
    ylab <- "Events per year"
  } else {
    ylab <- "Total event starts"
  }
  
  ggplot(monthly, aes(x = factor(start_month), y = n, fill = event_class)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_x_discrete(labels = month.abb) +
    scale_fill_manual(values = passaat_colors$event, name = "Class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ylab, title = "Monthly Distribution of Event Starts") +
    passaat_theme()
}

# -----------------------------------------------------------------------------
# Plot 4: DOY wind summary (mean + max)
# -----------------------------------------------------------------------------

#' Plot mean and maximum wind speed by day of year
#'
#' @param daily Daily data frame
#' @param smooth Apply loess smoothing
#' @param span Loess span (if smooth = TRUE)
#' @param thr_tc TC threshold
#' @param thr_hur Hurricane threshold
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
               color = passaat_colors$threshold["tc"], alpha = 0.7) +
    geom_hline(yintercept = thr_hur, linetype = "dashed", 
               color = passaat_colors$threshold["hur"], alpha = 0.7) +
    scale_color_manual(values = c(Mean = "grey40", Maximum = "steelblue"), name = NULL) +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = month.abb,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Wind speed (kt)", 
         title = "Daily Wind Speed by Day of Year") +
    passaat_theme()
}

# -----------------------------------------------------------------------------
# Plot 5: Monthly wind quantiles
# -----------------------------------------------------------------------------

#' Plot monthly wind quantile ribbon
#'
#' @param daily Daily data frame
#' @param probs Quantile probabilities (lower, median, upper)
#' @param log_scale Use log y-axis
#' @param thr_tc Tropical cyclone threshold (kt)
#' @param thr_hur Hurricane threshold (kt)
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
    geom_line(aes(y = median), linewidth = 1, color = passaat_colors$quantile["median"]) +
    geom_line(aes(y = upper), linewidth = 0.7, linetype = "dashed",
              color = passaat_colors$quantile["p95"]) +
    geom_line(aes(y = extreme), linewidth = 0.7, linetype = "dotted",
              color = passaat_colors$quantile["p99"]) +
    geom_hline(yintercept = thr_tc,  linetype = "dashed", 
               color = passaat_colors$threshold["tc"], alpha = 0.6) +
    geom_hline(yintercept = thr_hur, linetype = "dashed", 
               color = passaat_colors$threshold["hur"], alpha = 0.6) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(
      x = NULL, 
      y = "Wind speed (kt)",
      title = "Monthly Wind Distribution",
      subtitle = sprintf("Ribbon: median-P%d | Dotted: P%d", 
                         probs[2]*100, probs[3]*100)
    ) +
    passaat_theme()
  
  if (log_scale) p <- p + scale_y_log10()
  
  p
}

# -----------------------------------------------------------------------------
# Plot 6: Annual event count distribution
# -----------------------------------------------------------------------------
#' Plot histogram of annual event counts
#'
#' @param daily Daily data frame
#' @param metric "events" (distinct events) or "days" (event-days)
#' @param show_poisson Overlay Poisson distribution
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
    passaat_theme()
  
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

#' Scatter plot of event intensity vs duration
#'
#' @param daily Daily data frame (or events table directly)
#' @param events Pre-computed events table (optional)
#' @param thr_tc TC threshold
#' @param thr_hur Hurricane threshold
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
    scale_color_manual(values = passaat_colors$event, name = "Class") +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    labs(x = "Event duration (days)", y = "Peak wind (kt)",
         title = "Event Intensity vs Duration") +
    passaat_theme()
}

# -----------------------------------------------------------------------------
# Plot 8: Wind histogram / density
# -----------------------------------------------------------------------------

#' Plot wind speed distribution
#'
#' @param daily Daily data frame
#' @param type "histogram" or "density"
#' @param log_scale Use log x-axis
#' @param thr_tc Tropical cyclone threshold (kt)
#' @param thr_hur Hurricane threshold (kt)
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
               color = passaat_colors$threshold["tc"]) +
    geom_vline(xintercept = thr_hur, linetype = "dashed", 
               color = passaat_colors$threshold["hur"]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Wind speed (kt)", y = if(type == "histogram") "Count" else "Density",
         title = "Wind Speed Distribution") +
    passaat_theme()
  
  if (log_scale) p <- p + scale_x_log10()
  
  p
}

# -----------------------------------------------------------------------------
# Plot 9: Return level plot
# -----------------------------------------------------------------------------

#' Plot empirical return levels
#'
#' @param daily Daily data frame
#' @param block_maxima Use annual maxima (TRUE) or all exceedances (FALSE
#' @param threshold Threshold for POT approach (if block_maxima = FALSE)
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
    passaat_theme()
}

# -----------------------------------------------------------------------------
# Quick summary panel
# -----------------------------------------------------------------------------

#' Generate multi-panel summary figure
#'
#' @param daily Daily data frame
#' @param thr_tc TC threshold
#' @param thr_hur Hurricane threshold
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
