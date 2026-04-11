

# %%
library(ipdcstorm)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# Paths and output folders
# ------------------------------------------------------------------------------

ibtracs_file_path <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
baseline_output_dir <- "output/baseline"
validation_output_dir <- "output/validation"

seed <- 2026

# %%

# ------------------------------------------------------------------------------
# Target locations
# ------------------------------------------------------------------------------

targets <- tibble::tribble(
  ~name         , ~lat    , ~lon     ,
  "St_Martin"   , 18.0708 , -63.0501 ,
  "Saba"        , 17.6350 , -63.2300 ,
  "Statia"      , 17.4890 , -62.9740 ,
  "Puerto_Rico" , 18.2208 , -66.5901 ,
  "Miami"       , 25.7617 , -80.1918
)

# ------------------------------------------------------------------------------
# Configure and run the model
# ------------------------------------------------------------------------------

hazard_cfg <- make_hazard_cfg(
  data_path = ibtracs_file_path,
  search_radius_km = 600,
  historical_start_year = 1970L,
  simulation_years = 2000L,
  climate = make_climate_cfg(
    scenario = "stationary",
    delta_sst = 0,
    target_year = 2025
  )
)

output_hazard <- run_hazard_model(
  cfg = hazard_cfg,
  targets = targets,
  seed = seed,
  verbose = TRUE
)


output_daily <- generate_daily_hazard_impact(
  out = output_hazard,
  location = targets$name,
  sim_years = seq_len(hazard_cfg$n_sim),
  year0 = hazard_cfg$start_year,
  gust_factor = 1.3,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = "stationary",
  seed = seed
)
