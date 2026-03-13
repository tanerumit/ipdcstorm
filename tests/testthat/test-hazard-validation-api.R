test_that("hazard config maps canonical constructor names to runtime fields", {
  cfg_default <- make_hazard_cfg()
  expect_identical(cfg_default$start_year, 1970L)
  expect_identical(cfg_default$n_sim, 1000L)
  expect_s3_class(cfg_default$climate, "climate_cfg")
  expect_identical(cfg_default$climate$scenario, "stationary")
  expect_identical(cfg_default$climate$sensitivity_mode, "fixed")
  expect_no_error(ipdcstorm:::.normalize_climate_cfg(cfg_default$climate))

  cfg_custom <- make_hazard_cfg(
   historical_start_year = 1985L,
    simulation_years = 250L,
    climate = make_climate_cfg(scenario = "ssp245", target_year = 2050)
  )
  expect_identical(cfg_custom$start_year, 1985L)
  expect_identical(cfg_custom$n_sim, 250L)
  expect_identical(cfg_custom$climate$scenario, "ssp245")
  expect_false("n_sim_years" %in% names(cfg_custom))
})

test_that("validation suite inherits n_sim only from model output config", {
  out_stub <- list(
    config = list(n_sim = 250L),
    cfg = list(),
    events = tibble::tibble(
      location = character(0),
      storm_class = character(0),
      storm_id = character(0),
      year = integer(0),
      peak_wind_kt = numeric(0)
    ),
    rates = tibble::tibble(
      location = character(0),
      storm_class = character(0),
      lambda = numeric(0),
      n_years = integer(0)
    ),
    fit = tibble::tibble(),
    trackpoints = list()
  )

  cfg_inherit <- make_validation_cfg(n_sim = NULL, save_plots = FALSE, save_tables = FALSE)
  expect_message(
    run_validation_suite(out_stub, cfg = cfg_inherit),
    "Sim: 250 yr \\(model output config\\$n_sim\\)"
  )
})

test_that("run_validation_suite honors explicit validation n_sim overrides", {
  out_stub <- list(
    config = list(n_sim = 250L),
    cfg = list(),
    events = tibble::tibble(
      location = character(0),
      storm_class = character(0),
      storm_id = character(0),
      year = integer(0),
      peak_wind_kt = numeric(0)
    ),
    rates = tibble::tibble(
      location = character(0),
      storm_class = character(0),
      lambda = numeric(0),
      n_years = integer(0)
    ),
    fit = tibble::tibble(),
    trackpoints = list()
  )

  cfg_inherit <- make_validation_cfg(n_sim = NULL, save_plots = FALSE, save_tables = FALSE)
  expect_message(
    run_validation_suite(out_stub, cfg = cfg_inherit),
    "Sim: 250 yr \\(model output config\\$n_sim\\)"
  )

  cfg_override <- make_validation_cfg(n_sim = 180L, save_plots = FALSE, save_tables = FALSE)
  expect_message(
    run_validation_suite(out_stub, cfg = cfg_override),
    "Sim: 180 yr \\(validation_cfg\\$n_sim\\)"
  )
})

test_that("storm-class helpers no longer normalize legacy aliases", {
  expect_identical(ipdcstorm:::.normalize_storm_classes("TS34plus"), "TS34plus")

  rate_tbl <- tibble::tibble(
    location = c("Saba", "Saba"),
    storm_class = c("TS", "HUR"),
    lambda_model_raw = c(1.2, 0.3),
    lambda_ref = c(2.0, 0.55),
    expected_ratio = c(0.55, 0.30)
  )
  scalers <- ipdcstorm:::.lambda_scalers_from_rate_check(rate_tbl)
  expect_false(any(scalers$storm_class == "TS34plus"))

  rate_check <- ipdcstorm:::.build_rate_check_table(list(
    rates = tibble::tibble(
      location = "Saba",
      storm_class = c("TS", "HUR"),
      lambda = c(1.0, 0.2),
      n_years = c(100L, 100L)
    )
  ))
  expect_setequal(unique(rate_check$storm_class), c("TS", "HUR"))
  expect_false(any(rate_check$storm_class == "TS34plus"))
  expect_false(any(get_reference_rates()$storm_class == "TS34plus"))
})
