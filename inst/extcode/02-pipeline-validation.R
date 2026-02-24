
### RUN VALIDATION SUITE
### Convenience wrapper around `validate_hazard_model()`.
###
### - Sources the hazard model code
### - Defines configuration + targets
### - Runs end-to-end validation
### - Saves figures + tables to a single output directory

rm(list = ls())

library(tibble)

# Source model code (assumes repository root is the working directory)
source("R/hazard_core.R")
source("R/hazard_climate.R")
source("R/hazard_downscale.R")
source("R/hazard_ibtracs.R")
source("R/hazard_utils.R")
source("R/hazard_run.R")
source("R/hazard_validation.R")

set.seed(123)

# =============================================================================
# Configuration (same as run_hazard_model.R)
# =============================================================================

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
  "Miami",       25.7617,  -80.1918
)

per_target_cfg <- list(
  Saba      = list(thr_port = 40, thr_infra = 55),
  Statia    = list(thr_port = 38, thr_infra = 52),
  St_Martin = list(thr_port = 45, thr_infra = 60)
)

# =============================================================================
# SST Configuration (Level 1 Climate Modification)
# =============================================================================

sst_cfg <- make_sst_cfg(
  enabled = TRUE,
  sst_source = "builtin",            # uses built-in ERSST v5 MDR data
  baseline_years = 1991L:2020L,
  scenario = "stationary",            # stationary for validation
  scenario_start_year = 2025L,
  advanced = list(
    beta_sst = NULL,                  # NULL = estimate from data
    beta_prior = 0.6,                 # literature prior for shrinkage
    gamma_intensity = NULL,           # NULL = estimate from data
    gamma_prior = 0.065               # literature prior (6.5% per C)
  )
)

# =============================================================================
# Run end-to-end validation (model + suite + artifact saving)
# =============================================================================

res <- validate_hazard_model(
  cfg = cfg,
  targets = targets,
  severities = c("TS", "HUR"),
  sst_cfg = sst_cfg,
  holdout_years = 10,
  sim_years = 5000,
  return_periods = c(5, 10, 25, 50),
  seed = 42,
  out_dir = "output/validation",
  save_plots = TRUE,
  save_tables = TRUE
)

# Convenience handles for interactive use
out <- res$out
val <- res$val

message("\nSaved artifacts:")
print(res$artifacts)

