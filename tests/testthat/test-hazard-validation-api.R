test_that("hazard config normalizes n_sim with legacy compatibility", {
  cfg_default <- make_hazard_cfg()
  expect_identical(cfg_default$n_sim, cfg_default$n_sim_years)

  cfg_legacy <- make_hazard_cfg(n_sim_years = 250L)
  expect_identical(cfg_legacy$n_sim, 250L)
  expect_identical(cfg_legacy$n_sim_years, 250L)

  cfg_both <- make_hazard_cfg(n_sim = 180L, n_sim_years = 250L)
  expect_identical(cfg_both$n_sim, 180L)
  expect_identical(cfg_both$n_sim_years, 180L)
})

test_that("run_validation_suite inherits n_sim from model output config and honors override", {
  out_stub <- list(
    config = list(n_sim = 250L),
    cfg = list(n_sim_years = 999L),
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

  out_legacy <- out_stub
  out_legacy$config <- list(n_sim_years = 220L)
  expect_message(
    run_validation_suite(out_legacy, cfg = cfg_inherit),
    "Sim: 220 yr"
  )
})

test_that("legacy severity alias still works and TS34plus is normalized away", {
  events <- tibble::tibble(
    year = c(2000L, 2000L, 2001L),
    storm_class = c("TS", "HUR", "TS"),
    storm_id = c("a", "b", "c")
  )

  expect_warning(
    counts <- compute_annual_counts(events, severities = "TS34plus"),
    "deprecated"
  )
  expect_true(all(counts$storm_class == "TS"))

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
