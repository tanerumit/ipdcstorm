# =============================================================================
# Stress-test visualization: wind time-series across baseline and KNMI scenarios
#
# Purpose:
#   Visualize the 5 selected stress-test years as daily wind time-series.
#   Each panel shows one synthetic stress-test trace (a selected sim year).
#   Within each panel, the stationary baseline is overlaid against four
#   KNMI'23 climate scenarios at the same year index and seed, so that any
#   difference is attributable to climate forcing rather than sampling noise.
#
# Panel layout:
#   - 5 facets (one per selected stress-test year), stacked vertically
#   - Each facet: baseline (black) + 4 KNMI'23 scenario lines (coloured)
#   - x-axis: day of year; y-axis: sustained wind speed (kt)
#   - Reference lines at TS (34 kt) and hurricane (64 kt) thresholds
#
# Prerequisites:
#   Run 05-stress-test-experiment.R first. The following objects must be
#   present in the workspace:
#     selected_ids    — integer[5]: stress-test year indices
#     QUANTILE_PROBS  — numeric[5]: quantile levels used for selection
#     KNMI_SCENARIOS  — character[4]: "knmi_Ld", "knmi_Ln", "knmi_Hd", "knmi_Hn"
#     SEED            — integer: random seed used for the baseline run
#     N_SIM           — integer: number of simulated years
#     baseline_stress — named list[location]: daily tibbles for 5 baseline years
#     cc_stress       — named list[scenario][location]: daily tibbles for 5 years
# =============================================================================

library(dplyr)
library(ggplot2)
library(ipdcstorm)

# =============================================================================
# 1) Plot settings
# =============================================================================

# Focal location for the panel plot — must be a name present in baseline_stress
focal_location <- "Saba"

# Wind speed thresholds (kt) shown as reference lines
thr_ts  <- 34L   # tropical storm
thr_hur <- 64L   # hurricane (Category 1)

# Scenario display labels for the legend
scenario_labels <- c(
  "Baseline" = "Baseline (stationary)",
  "knmi_Ld"  = "KNMI L-d (low warming, dry)",
  "knmi_Hd"  = "KNMI H-d (high warming, dry)"
)

# Color palette: baseline in black, low warming in blue, high warming in red
scenario_colors <- c(
  "Baseline" = "black",
  "knmi_Ld"  = "#4575b4",   # blue — low, dry
  "knmi_Hd"  = "#d73027"    # red  — high, dry
)

# Hurricane season window: June 1 (DOY 152) – Dec 31 (DOY 365), non-leap year
season_doy_start <- 152L   # June 1
season_doy_end   <- 365L   # Dec 31

# Month-start DOY values for the hurricane season months (Jun–Dec)
month_doy_all    <- cumsum(c(1L, 31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L))
season_months    <- 6:12                        # June through December
month_doy_starts <- month_doy_all[season_months]
month_labels     <- month.abb[season_months]

# =============================================================================
# 2) Assemble combined plotting data
# =============================================================================

# Baseline: already filtered to selected_ids inside baseline_stress
base_df <- baseline_stress[[focal_location]] |>
  mutate(
    scenario = "Baseline",
    doy      = as.integer(format(date, "%j"))
  )

# KNMI scenarios: generated with sim_years = selected_ids, so no further filter
knmi_list <- lapply(KNMI_SCENARIOS, function(scen) {
  cc_stress[[scen]][[focal_location]] |>
    mutate(
      scenario = scen,
      doy      = as.integer(format(date, "%j"))
    )
})
knmi_df <- bind_rows(knmi_list)

# Merge, restrict to hurricane season, and factorise scenario (baseline on top)
plot_df <- bind_rows(knmi_df, base_df) |>
  filter(doy >= season_doy_start, doy <= season_doy_end) |>
  mutate(
    scenario = factor(scenario, levels = c(KNMI_SCENARIOS, "Baseline"))
  )

# --- Panel labels ---
# One label per selected stress-test year, ordered by severity quantile.
# Format: "Trace 1 — Q50%" … "Trace 5 — Q99%"
trace_labels <- setNames(
  paste0(
    "Trace ", seq_along(selected_ids),
    "  \u2014  Q", round(QUANTILE_PROBS * 100), "%"
  ),
  as.character(selected_ids)
)

plot_df <- plot_df |>
  mutate(
    trace_label = trace_labels[as.character(sim_year)],
    trace_label = factor(trace_label, levels = trace_labels)
  )

# =============================================================================
# 3) Threshold annotation data (appears in all facets)
# =============================================================================

thr_ann <- data.frame(
  yintercept  = c(thr_ts, thr_hur),
  label       = c("TS (34 kt)", "HUR (64 kt)"),
  label_color = c("#E69F00", "#d73027")
)

# =============================================================================
# 4) Panel plot
# =============================================================================

p_wind <- ggplot(plot_df, aes(x = doy, y = wind_kt)) +

  # --- Threshold reference lines ---
  geom_hline(
    yintercept = thr_ts,
    linetype   = "dashed",
    colour     = "#E69F00",
    linewidth  = 0.35,
    alpha      = 0.8
  ) +
  geom_hline(
    yintercept = thr_hur,
    linetype   = "dashed",
    colour     = "#d73027",
    linewidth  = 0.35,
    alpha      = 0.8
  ) +

  # --- KNMI scenario traces (plotted first, below baseline) ---
  geom_line(
    data      = ~ filter(.x, scenario != "Baseline"),
    aes(colour = scenario),
    linewidth = 0.9,
    alpha     = 0.9
  ) +

  # --- Baseline trace (plotted on top, slightly thicker, fully opaque) ---
  geom_line(
    data      = ~ filter(.x, scenario == "Baseline"),
    aes(colour = scenario),
    linewidth = 1.2,
    alpha     = 1.0
  ) +

  # --- Color scale ---
  scale_colour_manual(
    values = scenario_colors,
    labels = scenario_labels,
    name   = NULL,
    # Ensure factor order in legend: Baseline first, then KNMI
    breaks = c("Baseline", KNMI_SCENARIOS)
  ) +

  # --- Axes ---
  scale_x_continuous(
    limits = c(season_doy_start, season_doy_end),
    breaks = month_doy_starts,
    labels = month_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +

  # --- Facets: one panel per selected stress-test year ---
  facet_wrap(
    ~ trace_label,
    nrow   = 5,
    ncol   = 1,
    scales = "fixed"   # shared y-axis so intensity is comparable across traces
  ) +

  # --- Labels ---
  labs(
    x = NULL,
    y = "Sustained wind speed (kt)",
    title = sprintf(
      "%s \u2014 Stress-test wind traces: baseline vs KNMI\u201923 scenarios",
      focal_location
    ),
    subtitle = sprintf(
      "5 high-impact years from %d synthetic years (seed = %d, target year %d) \u2014 hurricane season (Jun\u2013Dec)",
      N_SIM, SEED, TARGET_YEAR
    ),
    caption = "Dashed lines: TS threshold (34 kt, orange) and hurricane threshold (64 kt, red)"
  ) +

  # --- Theme ---
  theme_light(base_size = 10) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(1.4, "cm"),
    legend.text      = element_text(size = 8),
    strip.text       = element_text(size = 8.5, face = "bold"),
    axis.text.x      = element_text(size = 7.5),
    axis.text.y      = element_text(size = 7.5),
    panel.spacing.y  = unit(0.4, "lines"),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    plot.caption     = element_text(size = 7, colour = "grey50", hjust = 0)
  ) +
  guides(
    colour = guide_legend(
      nrow          = 1,
      override.aes  = list(linewidth = 1.1, alpha = 1)
    )
  )

print(p_wind)

# =============================================================================
# 5) Optional: save
# =============================================================================

# Adjust output path as needed
# ggsave(
#   "output/stress_test_wind_traces.png",
#   plot   = p_wind,
#   width  = 10,
#   height = 14,
#   dpi    = 300,
#   units  = "in"
# )
