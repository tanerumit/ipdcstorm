.make_climate_counts_fixture <- function() {
  years <- 1991:2020
  tibble::tibble(
    year = rep(years, each = 2),
    storm_class = rep(c("TS", "HUR"), times = length(years)),
    n_events = c(
      as.integer(6 + ((years - min(years)) %% 5)),
      as.integer(2 + ((years - min(years)) %% 3))
    )
  )
}

test_that("make_climate_cfg validates inputs and prints fields", {
  cfg <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "linear_shifted",
    k_beta = 0.2,
    k_gamma = 0.1,
    perturb = list(v_scale = 0.05)
  )

  expect_s3_class(cfg, "climate_cfg")
  expect_true(cfg$enabled)
  expect_equal(cfg$scenario, "ssp245")
  expect_equal(cfg$sensitivity_mode, "linear_shifted")
  expect_equal(cfg$k_beta, 0.2)
  expect_equal(cfg$k_gamma, 0.1)
  expect_equal(cfg[["perturb"]]$v_scale, 0.05)
  expect_equal(names(cfg[["perturb"]]), "v_scale")

  resolved <- resolve_climate_inputs(
    cfg,
    annual_counts = .make_climate_counts_fixture(),
    future_period = 2035:2065,
    verbose = FALSE
  )
  expect_equal(resolved$perturb$v_scale, 0.05)
  expect_equal(resolved$perturb$r_scale, default_cc_params()$r_scale)
  expect_equal(resolved$perturb$speed_scale, default_cc_params()$speed_scale)
  expect_equal(resolved$perturb$precip_scale, default_cc_params()$precip_scale)

  txt <- paste(capture.output(print(cfg)), collapse = "\n")
  expect_match(txt, "Climate configuration")
  expect_match(txt, "ssp245")
  expect_match(txt, "Sensitivity mode")
  expect_match(txt, "Storm perturbation")

  stationary <- make_climate_cfg(scenario = "stationary")
  expect_true(stationary$enabled)
  expect_equal(stationary$scenario, "stationary")
  expect_null(stationary[["perturb"]])

  expect_warning(
    expect_error(make_climate_cfg(k_beta = "abc"), "single finite numeric"),
    "NAs introduced by coercion"
  )
  expect_error(make_climate_cfg(sensitivity_mode = "bad"))
  expect_error(make_climate_cfg(k_gamma = NA_real_), "single finite numeric")
})

test_that("stationary scenario resolves baseline climate scalars through the common path", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(scenario = "stationary", sensitivity_mode = "fixed")

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    lambda_table = tibble::tibble(
      storm_class = c("TS", "HUR"),
      lambda = c(0.8, 0.2)
    ),
    future_period = 2035:2065,
    verbose = FALSE
  )

  expect_identical(resolved$climate_mode, "baseline")
  expect_equal(resolved$delta_sst, 0)
  expect_identical(resolved$beta_sst, resolved$beta_0)
  expect_identical(resolved$gamma, resolved$gamma_0)
  expect_equal(exp(resolved$beta_sst * resolved$delta_sst), 1)
  expect_equal(1 + resolved$gamma * resolved$delta_sst, 1)
  expect_true(is.finite(resolved$p_hur_base))
})

test_that("fixed mode returns historical baseline sensitivities unchanged", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(scenario = "ssp245", sensitivity_mode = "fixed")

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    future_period = 2035:2065,
    verbose = FALSE
  )

  expect_true(resolved$delta_sst > 0)
  expect_equal(resolved$beta_sst, resolved$beta_0)
  expect_equal(resolved$gamma, resolved$gamma_0)
})

test_that("linear-shifted mode applies expected sensitivity transforms", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "linear_shifted",
    k_beta = 0.3,
    k_gamma = -0.4
  )

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    future_period = 2035:2065,
    verbose = FALSE
  )

  expect_equal(
    resolved$beta_sst,
    max(0, resolved$beta_0 * (1 + climate$k_beta * resolved$delta_sst))
  )
  expect_equal(
    resolved$gamma,
    max(0, resolved$gamma_0 * (1 + climate$k_gamma * resolved$delta_sst))
  )
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

test_that("resolver uses scenario table delta_sst path", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(scenario = "ssp585")

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    future_period = 2075:2105,
    verbose = FALSE
  )

  expect_equal(
    resolved$delta_sst,
    get_scenario_delta("ssp585", future_period = 2075:2105)
  )
  expect_identical(resolved$climate_mode, "future")
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

test_that("run_hazard_model accepts only NULL or climate_cfg and climate runs stay deterministic", {
  cfg <- make_hazard_cfg(simulation_years = 120L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  baseline <- make_climate_cfg(scenario = "stationary")
  climate <- make_climate_cfg(scenario = "ssp245")

  base1 <- run_hazard_model(cfg, targets, climate = baseline, seed = 321L, verbose = FALSE)
  base2 <- run_hazard_model(cfg, targets, climate = baseline, seed = 321L, verbose = FALSE)
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
  expect_identical(attr(base1$fit, "climate_mode"), "baseline")
  expect_gt(attr(fut1$fit, "delta_sst"), 0)
  expect_identical(attr(fut1$fit, "climate_mode"), "future")
  expect_gt(mean(fut1$sim$n_total), mean(base1$sim$n_total))
  expect_false(all(diff(fut1$sim$n_total) == 0))
})

test_that("baseline climate run matches climate = NULL numerically under fixed seed", {
  cfg <- make_hazard_cfg(simulation_years = 120L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  out_null <- run_hazard_model(cfg, targets, climate = NULL, seed = 456L, verbose = FALSE)
  out_stationary <- run_hazard_model(
    cfg,
    targets,
    climate = make_climate_cfg(scenario = "stationary"),
    seed = 456L,
    verbose = FALSE
  )

  expect_identical(out_null$sim, out_stationary$sim)
  expect_identical(out_null$events, out_stationary$events)
  expect_equal(attr(out_stationary$fit, "delta_sst"), 0)
  expect_identical(attr(out_stationary$fit, "climate_mode"), "baseline")
  expect_identical(attr(out_stationary$fit, "climate_scenario"), "stationary")
  expect_identical(attr(out_null$fit, "climate_mode"), "off")
  expect_identical(attr(out_null$fit, "climate_scenario"), "off")
})

test_that("legacy climate constructors and print methods are removed", {
  ns <- asNamespace("ipdcstorm")

  expect_false(exists("make_climate_shift", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("make_climate_response", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("make_climate_input", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("prepare_climate", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("print.climate_shift", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("print.climate_response", envir = ns, mode = "function", inherits = FALSE))
  expect_false(exists("print.climate_input", envir = ns, mode = "function", inherits = FALSE))

  expect_error(make_climate_shift(delta_sst = 0.8), "could not find function")
  expect_error(make_climate_response(beta_sst = 0.4), "could not find function")
  expect_error(make_climate_input(list(), list()), "could not find function")
  expect_error(prepare_climate(make_climate_cfg()), "could not find function")
})

test_that("run_hazard_model rejects legacy climate object classes", {
  cfg <- make_hazard_cfg(simulation_years = 40L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  legacy_shift <- structure(list(delta_sst = 0.8), class = c("climate_shift", "list"))
  legacy_response <- structure(list(beta_sst = 0.5, gamma = 0.1), class = c("climate_response", "list"))
  legacy_input <- structure(
    list(shift = legacy_shift, response = legacy_response),
    class = c("climate_input", "list")
  )

  expect_error(
    run_hazard_model(cfg, targets, climate = legacy_shift, seed = 11L, verbose = FALSE),
    "climate must be NULL or inherit from \"climate_cfg\""
  )
  expect_error(
    run_hazard_model(cfg, targets, climate = legacy_response, seed = 11L, verbose = FALSE),
    "climate must be NULL or inherit from \"climate_cfg\""
  )
  expect_error(
    run_hazard_model(cfg, targets, climate = legacy_input, seed = 11L, verbose = FALSE),
    "climate must be NULL or inherit from \"climate_cfg\""
  )
})

test_that("returned climate metadata uses only simplified resolved schema", {
  cfg <- make_hazard_cfg(simulation_years = 40L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  climate <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "linear_shifted",
    k_beta = 0.2,
    k_gamma = 0.1,
    perturb = list()
  )

  out <- run_hazard_model(cfg, targets, climate = climate, seed = 11L, verbose = FALSE)

  expect_type(out$cfg$climate, "list")
  expect_false(inherits(out$cfg$climate, "climate_input"))
  expect_false("shift" %in% names(out$cfg$climate))
  expect_false("response" %in% names(out$cfg$climate))
  expect_true(all(c(
    "enabled", "scenario", "source", "delta_sst", "baseline_years",
    "future_period", "sensitivity_mode", "k_beta", "k_gamma",
    "beta_0", "gamma_0", "beta_sst", "gamma", "p_hur_base",
    "climate_mode",
    "perturb", "perturb_state", "config"
  ) %in% names(out$cfg$climate)))
  expect_s3_class(out$cfg$climate$config, "climate_cfg")
  expect_identical(attr(out$fit, "climate_cfg"), out$cfg$climate$config)
})

test_that("baseline and future climate runs share resolved historical sensitivities", {
  cfg <- make_hazard_cfg(simulation_years = 120L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  climate_baseline <- make_climate_cfg(
    scenario = "stationary",
    sensitivity_mode = "fixed"
  )
  climate_future <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "fixed"
  )

  out_baseline <- run_hazard_model(cfg, targets, climate = climate_baseline, seed = 222L, verbose = FALSE)
  out_future <- run_hazard_model(cfg, targets, climate = climate_future, seed = 222L, verbose = FALSE)

  expect_identical(out_baseline$cfg$climate$climate_mode, "baseline")
  expect_identical(out_future$cfg$climate$climate_mode, "future")
  expect_identical(out_baseline$cfg$climate$beta_0, out_future$cfg$climate$beta_0)
  expect_identical(out_baseline$cfg$climate$gamma_0, out_future$cfg$climate$gamma_0)
  expect_identical(out_baseline$cfg$climate$beta_sst, out_future$cfg$climate$beta_sst)
  expect_identical(out_baseline$cfg$climate$gamma, out_future$cfg$climate$gamma)
  expect_equal(out_baseline$cfg$climate$delta_sst, 0)
  expect_gt(out_future$cfg$climate$delta_sst, 0)
  expect_false(identical(out_baseline$sim, out_future$sim))
})

test_that("fixed mode stays identical to neutral shifted mode under fixed seed", {
  cfg <- make_hazard_cfg(simulation_years = 120L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )
  climate_fixed <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "fixed"
  )
  climate_shift_neutral <- make_climate_cfg(
    scenario = "ssp245",
    sensitivity_mode = "linear_shifted",
    k_beta = 0,
    k_gamma = 0
  )

  out_fixed <- run_hazard_model(cfg, targets, climate = climate_fixed, seed = 222L, verbose = FALSE)
  out_shift_neutral <- run_hazard_model(cfg, targets, climate = climate_shift_neutral, seed = 222L, verbose = FALSE)

  expect_identical(out_fixed$sim, out_shift_neutral$sim)
  expect_identical(out_fixed$cfg$climate$beta_sst, out_shift_neutral$cfg$climate$beta_sst)
  expect_identical(out_fixed$cfg$climate$gamma, out_shift_neutral$cfg$climate$gamma)
})

test_that("transient SST workflow functions are deleted from the public surface", {
  expect_false(exists("make_sst_cfg", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false(exists("generate_sst_scenario", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false("sst_cfg" %in% names(formals(run_hazard_model)))
})
