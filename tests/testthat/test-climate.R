test_that("make_climate_shift validates inputs and prints fields", {
  shift <- make_climate_shift(delta_sst = 0.8)

  expect_s3_class(shift, "climate_shift")
  expect_equal(shift$delta_sst, 0.8)
  expect_equal(shift$baseline_years, 1991:2020)

  txt <- paste(capture.output(print(shift)), collapse = "\n")
  expect_match(txt, "delta_sst")
  expect_match(txt, "1991-2020")

  expect_error(make_climate_shift(delta_sst = "abc"), "single finite numeric")
  expect_error(make_climate_shift(delta_sst = c(0.5, 0.8)), "single finite numeric")
})

test_that("make_climate_response validates and fills perturb defaults", {
  response <- make_climate_response()
  expect_s3_class(response, "climate_response")
  expect_equal(response$beta_sst, 0)
  expect_equal(response$gamma, 0)
  expect_null(response$perturb)

  expect_error(make_climate_response(gamma = -0.1), "gamma must be >= 0")

  response2 <- make_climate_response(perturb = list(v_scale = 0.05))
  expect_equal(response2$perturb$v_scale, 0.05)
  expect_equal(response2$perturb$r_scale, default_cc_params()$r_scale)
  expect_equal(response2$perturb$speed_scale, default_cc_params()$speed_scale)
  expect_equal(response2$perturb$precip_scale, default_cc_params()$precip_scale)
})

test_that("make_climate_input validates classes and warns on no-op shifts", {
  shift <- make_climate_shift(delta_sst = 0.8)
  zero_response <- make_climate_response()

  expect_error(make_climate_input(list(), zero_response), "make_climate_shift")
  expect_error(make_climate_input(shift, list()), "make_climate_response")
  expect_warning(
    make_climate_input(shift, zero_response),
    "delta_sst is nonzero but all response parameters are zero"
  )
  expect_no_warning(
    make_climate_input(make_climate_shift(delta_sst = 0), zero_response)
  )

  climate <- make_climate_input(
    shift = shift,
    response = make_climate_response(beta_sst = 0.4, gamma = 0.05)
  )
  txt <- paste(capture.output(print(climate)), collapse = "\n")
  expect_match(txt, "delta_sst")
  expect_match(txt, "beta_sst")
})

test_that("get_scenario_delta interpolates and clamps by future-period midpoint", {
  info <- sst_scenario_info("all")
  ssp585 <- info[info$scenario == "ssp585", , drop = FALSE]

  expect_equal(
    get_scenario_delta("ssp585", future_period = 2035:2065),
    ssp585$delta_sst_2050[[1]]
  )

  expected_2090 <- ssp585$delta_sst_2050[[1]] +
    (2090 - 2050) / (2100 - 2050) *
    (ssp585$delta_sst_2100[[1]] - ssp585$delta_sst_2050[[1]])
  expect_equal(
    get_scenario_delta("ssp585", future_period = 2075:2105),
    expected_2090
  )

  expect_equal(
    get_scenario_delta("ssp585", future_period = 2100:2130),
    ssp585$delta_sst_2100[[1]]
  )

  expect_error(get_scenario_delta("nonexistent", future_period = 2035:2065), "Unknown SST scenario")
})

test_that("simulate_twolevel_counts uses scalar delta_sst", {
  lambda_table <- tibble::tibble(
    storm_class = c("TS", "HUR"),
    lambda = c(0.8, 0.2)
  )

  set.seed(1)
  sim0 <- simulate_twolevel_counts(
    lambda_table = lambda_table,
    k_hat = 1e6,
    n_years_sim = 5000,
    delta_sst = 0
  )

  set.seed(1)
  sim1 <- simulate_twolevel_counts(
    lambda_table = lambda_table,
    k_hat = 1e6,
    n_years_sim = 5000,
    delta_sst = 0.8,
    beta_sst = 0.5,
    gamma_intensity = 0.1
  )

  expect_true(all(c("sim_year", "activity_factor", "climate_scale", "activity_combined",
                    "p_hurricane", "n_total", "n_ts", "n_hur") %in% names(sim0)))
  expect_false("sst_anomaly" %in% names(sim0))
  expect_equal(attr(sim0, "delta_sst"), 0)
  expect_equal(attr(sim1, "delta_sst"), 0.8)
  expect_gt(mean(sim1$n_total), mean(sim0$n_total))
  expect_gt(mean(sim1$p_hurricane), mean(sim0$p_hurricane))
  expect_error(
    simulate_twolevel_counts(lambda_table, k_hat = 1e6, n_years_sim = 10, delta_sst = c(0.5, 0.8)),
    "single finite numeric"
  )
})

test_that("run_hazard_model consumes climate_input and records deterministic seeds", {
  cfg <- make_hazard_cfg(simulation_years = 120L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  climate <- make_climate_input(
    shift = make_climate_shift(delta_sst = 0.8),
    response = make_climate_response(beta_sst = 0.5, gamma = 0.1)
  )

  base1 <- run_hazard_model(cfg, targets, climate = NULL, seed = 321L, verbose = FALSE)
  base2 <- run_hazard_model(cfg, targets, climate = NULL, seed = 321L, verbose = FALSE)
  fut1 <- run_hazard_model(cfg, targets, climate = climate, seed = 321L, verbose = FALSE)
  fut2 <- run_hazard_model(cfg, targets, climate = climate, seed = 321L, verbose = FALSE)
  fut3 <- run_hazard_model(cfg, targets, climate = climate, seed = 654L, verbose = FALSE)
  auto <- run_hazard_model(cfg, targets, climate = climate, seed = NULL, verbose = FALSE)

  expect_identical(base1$sim, base2$sim)
  expect_identical(fut1$sim, fut2$sim)
  expect_false(identical(fut1$sim, fut3$sim))
  expect_true(is.integer(auto$run_metadata$seed))
  expect_true(length(auto$run_metadata$seed) == 1L)
  expect_true(is.finite(auto$run_metadata$seed))

  expect_equal(attr(base1$fit, "delta_sst"), 0)
  expect_equal(attr(fut1$fit, "delta_sst"), 0.8)
  expect_gt(mean(fut1$sim$n_total), mean(base1$sim$n_total))
  expect_false(all(diff(fut1$sim$n_total) == 0))
})

test_that("transient SST workflow functions are deleted from the public surface", {
  expect_false(exists("make_sst_cfg", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false(exists("generate_sst_scenario", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false("sst_cfg" %in% names(formals(run_hazard_model)))
})
