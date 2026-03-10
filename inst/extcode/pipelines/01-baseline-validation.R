
# ==============================================================================
# Model configuration ----------------------------------------------------------


# Required packages
library(ipdcstorm)
library(ggplot2)
library(dplyr)

# Files and folders
ibtracs_file_path <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
output_dir <- "output/baseline"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Target locations for storm hazard analysis
targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

# Model configuration parameters
model_cfg <- make_hazard_cfg(
  data_path       = ibtracs_file_path,
  search_radius_km = 800,
  start_year       = 1970L,
  n_sim_years      = 500L
)


# ==============================================================================


# 1) Run hazard model ----------------------------------------------------------

# Run hazard model
out <- run_hazard_model(
  cfg        = model_cfg,
  targets    = targets,
  sst_cfg    = NULL
)

# Validation
validation_cfg <- make_validation_cfg(
  holdout_years = 10L,
  n_sim = 500L,
  return_periods = c(5, 10, 25, 50),
  conf_level = 0.90,
  seed = 42L,
  out_dir = "output/validation",
  save_plots = TRUE,
  save_tables = TRUE,
  advanced = NULL
)

res <- run_validation_suite(
  out = out,
  cfg = validation_cfg
)

# Simulate daily hazard series and impacts
daily_all <- generate_daily_hazard_impact(
  out       = out,
  location  = targets$name,
  sim_years = seq_len(model_cfg$n_sim_years),
  year0     = model_cfg$start_year,
  gust_factor    = 1.3,   # sustained → 3-sec gust conversion
  damage_method  = "powerlaw",
  damage_params  = list(thr = 34, V_ref = 80, d_ref = 0.03, p = 3, d_max = 0.10),
  pulse_shape    = "cosine",
  scenario       = "stationary",
  seed           = 42
)

# ==============================================================================

# 2) Analyze results -----------------------------------------------------------

# Per-location
for (loc in names(daily_all)) {

  save_hazard_viz_plots(daily = daily_all[[loc]],
    output_dir = file.path(output_dir, gsub("[^A-Za-z0-9_]", "_", loc)),
    location_name = loc,
    width = 9, height = 6, dpi = 300, base_size = 11,
    thr_tc = 34, thr_hur = 64)

}

# Cross-location comparison
comp_dir <- file.path(output_dir, "comparison")
if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

# Observed annual Poisson rates by island
p_rates <- ggplot(out$rates , aes(x = location, y = lambda, fill = storm_class)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", lambda)),
            position = position_dodge(0.6), vjust = -0.4, size = 3.2) +
  scale_fill_manual(
    values = c(TS = "#E69F00", HUR = "#D55E00"),
    labels = c(TS = "Tropical Storm", HUR = "Hurricane"),
    name   = NULL) +
  labs(x = NULL, y = expression(lambda ~ "(events yr"^-1*")"),
       title = "Historical Annual Storm Rates by Island")

ggsave(file.path(comp_dir, "rate_comparison.png"), p_rates,
       width = 8, height = 4)



