
# Load packages & parameters
library(ipdcstorm)

# Parameters
cfg <- make_hazard_cfg(
  data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
  search_radius_km = 800,
  start_year = 1970,
  n_sim_years = 1000
)

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918)

per_target_cfg <- list(
  Saba    = list(thr_port = 40, thr_infra = 55),
  Statia  = list(thr_port = 38, thr_infra = 52),
  St_Martin = list(thr_port = 45, thr_infra = 60))

# Run baseline model
hazard_out <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = NULL
)

#print(out_stat$rates)

# Create daily results
daily <- generate_daily_hazard_impact(
  out           = hazard_out,
  location      = c("St_Martin", "Saba", "Statia", "Puerto_Rico", "Miami"),
  sim_years     = 1:200,         # use first 200 simulated years
  year0         = 2025,           # calendar year for sim_year = 1
  gust_factor   = 1.25,           # sustained-to-gust conversion
  damage_method = "powerlaw",
  seed          = 42,
  scenario      = "baseline"   # label (no climate mods applied)
)


# Validation
hazard_validate <- validate_hazard_model(
  cfg = cfg,
  targets = targets,
  per_target_cfg = per_target_cfg,
  severities = c("TS", "HUR64plus"),
  sst_cfg = sst_cfg,
  holdout_years = 10,
  n_sim = 5000,
  return_periods = c(5, 10, 25, 50),
  seed = 42,
  out_dir = "output/validation",
  save_plots = TRUE,
  save_tables = TRUE
)

################################################################################
################################################################################


