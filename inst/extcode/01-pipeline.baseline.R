
# Load packages & parameters
library(ipdcstorm)

# root for all saved outputs
out_dir <- "output/baseline"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1) Configuration -------------------------------------------------------

cfg <- make_hazard_cfg(
  data_path       = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  search_radius_km = 800,
  start_year       = 1970L,
  n_sim_years      = 1000L
)

# Stationary baseline — no climate modifications
sst_cfg <- make_sst_cfg(enabled = FALSE)

targets <- tibble::tibble(
  name = c("Saba", "St. Eustatius", "St. Martin"),
  lat  = c(17.63, 17.49, 18.04),
  lon  = c(-63.24, -62.97, -63.07)
)

# --- 2) Run hazard model ----------------------------------------------------

out <- run_hazard_model(
  cfg        = cfg,
  targets    = targets,
  sst_cfg    = sst_cfg
)

# --- 3) Generate daily hazard + impact series --------------------------------

daily_all <- generate_daily_hazard_impact(
  out       = out,
  location  = targets$name,
  sim_years = seq_len(cfg$n_sim_years),
  year0     = 2000L,                        # calendar anchor for synthetic years
  gust_factor    = 1.3,                     # sustained → 3-sec gust conversion
  damage_method  = "powerlaw",
  damage_params  = list(thr = 34, V_ref = 80, d_ref = 0.03, p = 3, d_max = 0.10),
  pulse_shape    = "cosine",
  scenario       = "stationary",
  seed           = 42
)

# --- 4) Per-location visualizations ------------------------------------------

loc <- "Saba"


for (loc in names(daily_all)) {

  daily <- daily_all[[loc]]                 # extract this location's daily tibble
  loc_dir <- file.path(out_dir, gsub("[^A-Za-z0-9_]", "_", loc))
  if (!dir.exists(loc_dir)) dir.create(loc_dir, recursive = TRUE)

  # -- P1: Wind time series (5-year window to keep readable) --
  yr_window <- 1:5                          # first 5 simulated years
  daily_window <- daily |> dplyr::filter(.data$sim_year %in% yr_window)
  p1 <- plot_wind_timeseries(
    daily  = daily_window,
    title  = paste(loc, "— Daily Wind (Years 1–5)")
  )
  ggsave(file.path(loc_dir, "01_wind_timeseries.png"), p1,
         width = 12, height = 4, dpi = 150)

  # -- P2: Seasonality (DOY histogram, event-days) --
  p2 <- plot_seasonality_doy(daily, metric = "event_days")
  ggsave(file.path(loc_dir, "02_seasonality_doy.png"), p2,
         width = 8, height = 6, dpi = 150)

  # -- P3: Monthly event counts (normalized to annual rate) --
  p3 <- plot_monthly_events(daily, normalize = TRUE)
  ggsave(file.path(loc_dir, "03_monthly_events.png"), p3,
         width = 8, height = 5, dpi = 150)

  # -- P4: DOY mean/max wind with loess smooth --
  p4 <- plot_doy_wind(daily, smooth = TRUE, span = 0.15)
  ggsave(file.path(loc_dir, "04_doy_wind.png"), p4,
         width = 8, height = 5, dpi = 150)

  # -- P5: Monthly wind quantile ribbons (P50, P95, P99) --
  p5 <- plot_monthly_quantiles(daily, probs = c(0.50, 0.95, 0.99))
  ggsave(file.path(loc_dir, "05_monthly_quantiles.png"), p5,
         width = 8, height = 5, dpi = 150)

  # -- P6: Annual event count distribution with Poisson overlay --
  p6 <- plot_annual_counts(daily, metric = "events", show_poisson = TRUE)
  ggsave(file.path(loc_dir, "06_annual_counts.png"), p6,
         width = 7, height = 5, dpi = 150)

  # -- P7: Intensity vs duration scatter --
  p7 <- plot_intensity_duration(daily = daily)
  ggsave(file.path(loc_dir, "07_intensity_duration.png"), p7,
         width = 7, height = 5, dpi = 150)

  # -- P8: Wind speed distribution (density) --
  p8_data <- daily |> dplyr::filter(.data$wind_kt > 0)  # exclude zero-wind days for clarity
  p8 <- plot_wind_distribution(p8_data, type = "density")
  ggsave(file.path(loc_dir, "08_wind_distribution.png"), p8,
         width = 7, height = 5, dpi = 150)

  # -- P9: Empirical return levels (annual maxima) --
  p9 <- plot_return_levels(daily, block_maxima = TRUE)
  ggsave(file.path(loc_dir, "09_return_levels.png"), p9,
         width = 7, height = 5, dpi = 150)

  # -- P10: Four-panel summary --
  p10 <- plot_summary_panel(daily)
  ggsave(file.path(loc_dir, "10_summary_panel.png"), p10,
         width = 14, height = 10, dpi = 150)

  message("[VIZ] Saved 10 plots for ", loc, " → ", loc_dir)
}

# --- 5) Cross-location comparison plots --------------------------------------

comp_dir <- file.path(out_dir, "comparison")
if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

# 5A: Rate table — observed annual Poisson rates by island
rates <- out$rates                          # columns: location, storm_class, lambda, …

p_rates <- ggplot(rates, aes(x = location, y = lambda, fill = storm_class)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", lambda)),
            position = position_dodge(0.6), vjust = -0.4, size = 3.2) +
  scale_fill_manual(
    values = c(TS = "#E69F00", HUR64plus = "#D55E00"),
    labels = c(TS = "Tropical Storm", HUR64plus = "Hurricane"),
    name   = NULL
  ) +
  labs(x = NULL, y = expression(lambda ~ "(events yr"^-1*")"),
       title = "Historical Annual Storm Rates by Island") +
  passaat_theme()

ggsave(file.path(comp_dir, "rate_comparison.png"), p_rates,
       width = 8, height = 5, dpi = 150)

# 5B: Historical event peak winds across islands
events_hist <- out$events                   # columns: location, storm_id, peak_wind_kt, …

p_peaks <- ggplot(events_hist, aes(x = location, y = peak_wind_kt, fill = storm_class)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.size = 1.2) +
  scale_fill_manual(
    values = c(TD = "grey70", TS = "#E69F00", HUR64plus = "#D55E00"),
    name   = NULL
  ) +
  labs(x = NULL, y = "Peak site wind (kt)",
       title = "Distribution of Historical Event Intensities") +
  passaat_theme()

ggsave(file.path(comp_dir, "peak_wind_boxplot.png"), p_peaks,
       width = 8, height = 5, dpi = 150)

# 5C: Return level comparison across islands (overlay)
rl_data <- purrr::imap_dfr(daily_all, function(daily, loc) {
  daily |>
    dplyr::group_by(.data$sim_year) |>
    dplyr::summarise(max_wind = max(.data$wind_kt, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$max_wind)) |>
    dplyr::mutate(
      rank          = dplyr::row_number(),
      n             = dplyr::n(),
      return_period = (n + 1) / rank,       # Weibull plotting position
      location      = loc
    )
})

p_rl <- ggplot(rl_data, aes(x = return_period, y = max_wind, color = location)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.4) +
  scale_x_log10(
    breaks = c(2, 5, 10, 25, 50, 100, 250, 500, 1000),
    labels = c("2", "5", "10", "25", "50", "100", "250", "500", "1000")
  ) +
  geom_hline(yintercept = 64, linetype = "dashed", alpha = 0.5) + # hurricane threshold
  annotate("text", x = 2, y = 67, label = "Hurricane (64 kt)",
           hjust = 0, size = 3, color = "grey40") +
  labs(x = "Return period (years)", y = "Annual max wind (kt)",
       title = "Empirical Return Levels — All Islands",
       color = NULL) +
  passaat_theme()

ggsave(file.path(comp_dir, "return_level_comparison.png"), p_rl,
       width = 9, height = 5.5, dpi = 150)

# 5D: Simulated annual count distributions side-by-side
sim_summary <- out$sim |>                   # columns: location, sim_year, n_total, n_ts, n_hur
  dplyr::select("location", "sim_year", "n_total") |>
  dplyr::group_by(.data$location) |>
  dplyr::mutate(mean_n = mean(.data$n_total, na.rm = TRUE))

p_sim <- ggplot(sim_summary, aes(x = n_total)) +
  geom_bar(fill = "steelblue", color = "white", width = 0.8) +
  geom_vline(aes(xintercept = mean_n), linetype = "dashed", color = "red") +
  facet_wrap(~ location, ncol = 1, scales = "free_y") +
  labs(x = "Total events per year", y = "Count (simulated years)",
       title = "Simulated Annual Event Totals",
       subtitle = sprintf("n = %s synthetic years | dashed = mean",
                          format(cfg$n_sim_years, big.mark = ","))) +
  passaat_theme()

ggsave(file.path(comp_dir, "sim_annual_totals.png"), p_sim,
       width = 7, height = 8, dpi = 150)

message("\n[DONE] All baseline outputs saved to: ", out_dir)
