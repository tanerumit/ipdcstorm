
# =============================================================================
# Script overview: visualization and derived summaries.
# - No function definitions in this script.
# - Builds daily/weekly/yearly/event summaries from `daily_saba_impact` and
#   produces ggplot2 figures for seasonality, intensity, and duration.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

ggtheme <- theme_light()
event_colors <- c("TC" = "orange", "HUR" = "red")

# Process outputs --------------------------------------------------------------


# -----------------------------------------------------------------------------#
# 0) Base daily table + standard derived time keys
# -----------------------------------------------------------------------------#
daily <- daily_saba_impact %>%
  mutate(
    doy   = as.integer(format(date, "%j")),
    month = as.integer(format(date, "%m")),
    iso_weekday = as.integer(format(date, "%u")),                 # 1=Mon,...7=Sun
    week_start  = as.Date(date) - (iso_weekday - 1L),             # ISO week Monday
    week_end    = week_start + 6L,
    week_num    = as.integer(format(date, "%V")),                 # ISO week number
    week_year   = as.integer(format(date, "%G"))                  # ISO week-year
  )

# -----------------------------------------------------------------------------#
# 1) Event-level table (one row per event_id per sim_year)
# -----------------------------------------------------------------------------#
events <- daily %>%
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
  mutate(event_class = ifelse(event_class == "HUR", "HUR", "TC"))

# Event starts (explicit)
event_starts <- events %>%
  transmute(
    location, sim_year, event_id,
    event_class,
    start_date,
    start_doy   = as.integer(format(start_date, "%j")),
    start_month = as.integer(format(start_date, "%m"))
  )

# -----------------------------------------------------------------------------#
# 2) Weekly summary (one row per week)
# -----------------------------------------------------------------------------#
weekly <- daily %>%
  group_by(location, sim_year, week_year, week_num, week_start, week_end) %>%
  summarise(
    # hazard (weekly intensity)
    wind_max_kt = max(wind_kt, na.rm = TRUE),
    wind_p95_kt = quantile(wind_kt, 0.95, na.rm = TRUE),
    
    # hazard occurrence / persistence
    n_event_days  = sum(!is.na(event_id)),
    n_events_week = n_distinct(event_id[!is.na(event_id)]),
    any_tc_week   = any(is_tc_day %in% TRUE),
    any_hur_week  = any(is_hur_day %in% TRUE),
    
    # disruptions
    port_disrupt_any  = if ("port_disrupt" %in% names(cur_data_all())) any(port_disrupt %in% TRUE) else NA,
    infra_disrupt_any = if ("infra_disrupt" %in% names(cur_data_all())) any(infra_disrupt %in% TRUE) else NA,
    
    # damage (assumes these columns exist in intensity method)
    damage_rate_sum = if ("damage_rate" %in% names(cur_data_all())) sum(damage_rate, na.rm = TRUE) else NA_real_,
    damage_inc_sum  = if ("damage_rate" %in% names(cur_data_all())) sum(damage_rate, na.rm = TRUE) else NA_real_,
    cum_damage_end  = if ("cum_damage" %in% names(cur_data_all())) last(cum_damage) else NA_real_,
    
    .groups = "drop"
  )

# NOTE: cur_data_all() is current (dplyr >=1.1). If you want strict schema,
# remove the conditionals and force columns to exist.

# -----------------------------------------------------------------------------#
# 3) Yearly summary (one row per sim_year)
# -----------------------------------------------------------------------------#
yearly <- daily %>%
  group_by(location, sim_year) %>%
  summarise(
    # hazard (annual intensity)
    wind_max_kt = max(wind_kt, na.rm = TRUE),
    wind_p99_kt = quantile(wind_kt, 0.99, na.rm = TRUE),
    
    # hazard occurrence
    n_events_year = n_distinct(event_id[!is.na(event_id)]),
    n_event_days  = sum(!is.na(event_id)),
    any_tc_year   = any(is_tc_day %in% TRUE),
    any_hur_year  = any(is_hur_day %in% TRUE),
    
    # disruptions
    port_disrupt_any  = if ("port_disrupt" %in% names(cur_data_all())) any(port_disrupt %in% TRUE) else NA,
    infra_disrupt_any = if ("infra_disrupt" %in% names(cur_data_all())) any(infra_disrupt %in% TRUE) else NA,
    
    # damage
    damage_rate_sum = if ("damage_rate" %in% names(cur_data_all())) sum(damage_rate, na.rm = TRUE) else NA_real_,
    damage_inc_sum  = if ("damage_rate" %in% names(cur_data_all())) sum(damage_rate, na.rm = TRUE) else NA_real_,
    cum_damage_end  = if ("cum_damage" %in% names(cur_data_all())) last(cum_damage) else NA_real_,
    
    .groups = "drop"
  )

# -----------------------------------------------------------------------------#
# 4) Seasonality summaries (use event_starts / daily once)
# -----------------------------------------------------------------------------#

# Day-of-year: event-days (duration-weighted)
doy_event_days <- daily %>%
  filter(!is.na(event_id), !is.na(event_class)) %>%
  mutate(event_class = ifelse(event_class == "HUR", "HUR", "TC")) %>%
  count(event_class, doy, name = "n_event_days")

# Month-of-year: event starts (initiation only; keeps all months)
monthly_event_starts <- event_starts %>%
  count(event_class, start_month, name = "n_event_starts")

monthly_event_starts_full <- tidyr::expand_grid(
  event_class = c("TC", "HUR"),
  start_month = 1:12
) %>%
  left_join(monthly_event_starts, by = c("event_class", "start_month")) %>%
  mutate(n_event_starts = tidyr::replace_na(n_event_starts, 0L))

# DOY wind summary (mean + max)
doy_wind_both <- daily %>%
  group_by(doy) %>%
  summarise(
    mean_wind_kt = mean(wind_kt, na.rm = TRUE),
    max_wind_kt  = max(wind_kt, na.rm = TRUE),
    .groups = "drop"
  )

# Monthly wind quantiles (mean, median, p99)
monthly_wind_quantile <- daily %>%
  group_by(month) %>%
  summarise(
    mean   = mean(wind_kt, na.rm = TRUE),
    median = median(wind_kt, na.rm = TRUE),
    p99    = quantile(wind_kt, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

# Annual event-count distribution input (for compounding)
year_event_counts <- yearly %>%
  transmute(sim_year, n_events = n_events_year, n_days = n_event_days)


daily %>% summarise(n = n(), date_min = min(date), date_max = max(date))
weekly %>% summarise(n = n(), week_min = min(week_start), week_max = max(week_start))
yearly %>% summarise(n = n(), sim_year_min = min(sim_year), sim_year_max = max(sim_year))

yearly %>%
  summarise(
    mean_events = mean(n_events_year),
    p95_events  = quantile(n_events_year, 0.95),
    mean_wmax   = mean(wind_max_kt),
    p99_wmax    = quantile(wind_max_kt, 0.99, na.rm = TRUE)
  )



################################################################################

# -----------------------------------------------------------------------------
#  VISUALIZATIONS 
# -----------------------------------------------------------------------------


# 1) Daily wind time series with event segment

p1 <- ggplot(daily, aes(x = date, y = wind_kt)) +
  geom_line(color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = thr_tc,  linetype = "dashed", color = "orange") +
  geom_hline(yintercept = hurricane_threshold_kt, linetype = "dashed", color = "red") +
  geom_segment(
    data = events,
    aes(x = start_date, xend = end_date, y = max_wind_kt, yend = max_wind_kt, 
        color = event_class),
    linewidth = 1.2, lineend = "round",
    inherit.aes = FALSE, alpha = 0.8
  ) +
  scale_color_manual(values = event_colors) +
  scale_x_date(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Date", y = "Wind speed (kt)",
    title = "Daily Wind with Tropical Cyclone and Hurricane Events",
    subtitle = "Segments show event duration; segment height shows realised peak intensity"
  ) +
  ggtheme

  ggsave("output/plots/hurricane_event_ts.png", p1, height = 5, width = 8)

# 2) DOY event frequency (event-days) (uses doy_event_days)

p2 <- ggplot(doy_event_days, aes(x = doy, y = n_event_days, fill = event_class)) +
  geom_col() +
  facet_wrap(~ event_class, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = event_colors) +
  scale_x_continuous(breaks = seq(0, 366, by = 30), limits = c(1, 366), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Day of year (1–366)",
    y = "Count of event-days across synthetic record",
    title = "Seasonality of Synthetic TC and Hurricane Activity"
  ) +
  ggtheme

ggsave("output/plots/hurricane_event_doy.png", height = 5, width = 8)

# 3) Monthly event starts (initiation) (uses monthly_event_starts_full)
p3 <- ggplot(monthly_event_starts_full,
             aes(x = factor(start_month, levels = 1:12), y = n_event_starts, fill = event_class)) +
  geom_col() +
  facet_wrap(~ event_class, ncol = 1, scales = "free_y") +
  scale_x_discrete(drop = FALSE, labels = month.abb) +
  scale_fill_manual(values = event_colors) +
  labs(
    x = "Month",
    y = "Count of event starts across synthetic record",
    title = "Monthly Seasonality of Synthetic TC and Hurricane Event Starts"
  ) +
  ggtheme


ggsave("output/plots/hurricane_event_seasonality.png", height = 5, width = 8)


#4) Event intensity vs duration (uses events)

p4 <- ggplot(events, aes(x = dur_days, y = max_wind_kt, color = event_class)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = thr_tc,  linetype = "dashed") +
  geom_hline(yintercept = hurricane_threshold_kt, linetype = "dashed") +
  scale_color_manual(values = event_colors) +
  labs(
    x = "Event duration (days)",
    y = "Maximum wind (kt)",
    color = "Class",
    title = "Event Intensity vs Duration (TC vs Hurricane)"
  ) +
  ggtheme


#5) Mean and max wind by day-of-year (uses doy_wind_both)
p5 <- ggplot(doy_wind_both, aes(x = doy)) +
  geom_line(aes(y = mean_wind_kt), linewidth = 0.7) +
  geom_line(aes(y = max_wind_kt),  linewidth = 0.7) +
  geom_hline(yintercept = thr_tc,  linetype = "dashed", color = "orange") +
  geom_hline(yintercept = hurricane_threshold_kt, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0)) +
  labs(
    x = "Day of year",
    y = "Wind speed (kt)",
    title = "Mean and Maximum Daily Wind Speed by Day of Year"
  ) +
  ggtheme
p5


#6) Monthly wind quantiles with log y (uses monthly_wind_quantile)
eps <- 0.1

p7 <- ggplot(monthly_wind_quantile, aes(x = month)) +
  geom_ribbon(
    aes(ymin = pmax(mean, eps), ymax = pmax(p99, eps)),
    fill = "grey80"
  ) +
  geom_line(aes(y = pmax(median, eps)), linewidth = 0.9) +
  geom_line(aes(y = pmax(p99,    eps)), linewidth = 0.7) +
  geom_hline(yintercept = thr_tc,  linetype = "dashed", color = "orange") +
  geom_hline(yintercept = hurricane_threshold_kt, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(1, 12)) +
  scale_y_log10() +
  labs(
    x = "Month",
    y = "Wind speed (kt, log scale)",
    title = "Monthly Wind Distribution (Median and P99)"
  ) +
  ggtheme
p7

#7) Distribution of annual event counts (compounding) (uses year_event_counts)
p8 <- ggplot(year_event_counts, aes(x = n_events)) +
  geom_bar() +
  scale_x_continuous(breaks = 0:max(year_event_counts$n_events)) +
  labs(
    x = "Number of events per year",
    y = "Number of synthetic years",
    title = "Distribution of Annual Tropical Cyclone Events"
  ) +
  ggtheme

ggsave("output/plots/hurricane_event_histogram.png", height = 5, width = 8)

#7) Distribution of annual event counts (compounding) (uses year_event_counts)
p9 <- ggplot(year_event_counts, aes(x = n_days)) +
  geom_bar() +
  scale_x_continuous(breaks = 0:max(year_event_counts$n_days)) +
  labs(
    x = "Number of days per year",
    y = "Number of synthetic years",
    title = "Distribution of Days with Annual Tropical Cyclones"
  ) +
  ggtheme

ggsave("output/plots/hurricane_duration_histogram.png", height = 5, width = 8)

