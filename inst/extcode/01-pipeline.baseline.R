
# Load packages & parameters
library(ipdcstorm)

# Parameters
cfg <- make_hazard_cfg(
  data_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
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
out_stat <- run_hazard_model(
  cfg = cfg, targets = targets, per_target_cfg = per_target_cfg,
  sst_cfg = NULL
)

cat("\nStationary lambdas:\n")
print(out_stat$rates)

