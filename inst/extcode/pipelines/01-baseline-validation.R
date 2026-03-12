# ==============================================================================
# Baseline hazard validation pipeline
# ==============================================================================

library(ipdcstorm)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# Paths and output folders
# ------------------------------------------------------------------------------

ibtracs_file_path <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
baseline_output_dir <- "output/baseline"
validation_output_dir <- "output/validation"

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

ensure_dir(baseline_output_dir)
ensure_dir(validation_output_dir)


# ------------------------------------------------------------------------------
# Target locations
# ------------------------------------------------------------------------------

targets <- tibble::tribble(
  ~location,      ~lat,     ~lon,
  "St_Martin",    18.0708, -63.0501,
  "Saba",         17.6350, -63.2300,
  "Statia",       17.4890, -62.9740,
  "Puerto_Rico",  18.2208, -66.5901,
  "Miami",        25.7617, -80.1918
)


# ------------------------------------------------------------------------------
# Model configuration
# ------------------------------------------------------------------------------

hazard_cfg <- make_hazard_cfg(
  data_path = ibtracs_file_path,
  search_radius_km = 800,
  historical_start_year = 1970L,
  simulation_years = 500L
)

validation_cfg <- make_validation_cfg(
  holdout_years = 10L,
  n_sim = 500L,
  return_periods = c(5, 10, 25, 50),
  conf_level = 0.90,
  seed = 42L,
  out_dir = validation_output_dir,
  save_plots = TRUE,
  save_tables = TRUE,
  advanced = NULL
)


# ------------------------------------------------------------------------------
# Run baseline hazard model and validation
# ------------------------------------------------------------------------------

hazard_out <- run_hazard_model(
  cfg = hazard_cfg,
  targets = targets,
  seed = 42L,
  climate = make_climate_cfg(scenario = "stationary")
)

validation_out <- run_validation_suite(
  out = hazard_out,
  cfg = validation_cfg
)


# ------------------------------------------------------------------------------
# Generate daily hazard-impact series
# ------------------------------------------------------------------------------

daily_by_location <- generate_daily_hazard_impact(
  out = hazard_out,
  location = targets$location,
  sim_years = seq_len(hazard_cfg$n_sim),
  year0 = hazard_cfg$start_year,
  gust_factor = 1.3,
  damage_method = "powerlaw",
  damage_params = list(
    thr = 34,
    V_ref = 80,
    d_ref = 0.03,
    p = 3,
    d_max = 0.10
  ),
  pulse_shape = "cosine",
  scenario = "stationary",
  seed = 42L
)


# ------------------------------------------------------------------------------
# Save per-location visualizations
# ------------------------------------------------------------------------------

for (location in names(daily_by_location)) {
  location_output_dir <- file.path(
    baseline_output_dir,
    gsub("[^A-Za-z0-9_]", "_", location)
  )

  save_hazard_viz_plots(
    daily = daily_by_location[[location]],
    output_dir = location_output_dir,
    location_name = location,
    width = 9,
    height = 6,
    dpi = 300,
    base_size = 11,
    thr_tc = 34,
    thr_hur = 64
  )
}


# ------------------------------------------------------------------------------
# Save cross-location comparison plot
# ------------------------------------------------------------------------------

comparison_output_dir <- file.path(baseline_output_dir, "comparison")
ensure_dir(comparison_output_dir)

rate_comparison_plot <- ggplot(
  hazard_out$rates,
  aes(x = location, y = lambda, fill = storm_class)
) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", lambda)),
    position = position_dodge(0.6),
    vjust = -0.4,
    size = 3.2
  ) +
  scale_fill_manual(
    values = c(TS = "#E69F00", HUR = "#D55E00"),
    labels = c(TS = "Tropical Storm", HUR = "Hurricane"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = expression(lambda ~ "(events yr"^-1 * ")"),
    title = "Historical Annual Storm Rates by Island"
  )

ggsave(
  filename = file.path(comparison_output_dir, "rate_comparison.png"),
  plot = rate_comparison_plot,
  width = 8,
  height = 4
)


# Keep key results available for interactive use
baseline_validation <- list(
  hazard = hazard_out,
  validation = validation_out,
  daily = daily_by_location
)
