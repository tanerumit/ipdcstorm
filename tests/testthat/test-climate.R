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

.hazard_test_targets <- tibble::tibble(
  name = "Saba",
  lat = 17.63,
  lon = -63.23
)

.hazard_test_targets_dup <- tibble::tibble(
  name = c("Saba_A", "Saba_B"),
  lat = c(17.63, 17.63),
  lon = c(-63.23, -63.23)
)

.run_hazard_reference <- local({
  cache <- new.env(parent = emptyenv())

  function(scenario = "stationary",
           delta_sst = NULL,
           target_year = NULL,
           sensitivity_mode = "fixed",
           k_beta = 0,
           k_gamma = 0,
           perturb = NULL,
           simulation_years = 8L,
           seed = 222L) {
    perturb_key <- if (is.null(perturb)) {
      "NULL"
    } else {
      paste(names(perturb), unlist(perturb, use.names = FALSE), sep = "=", collapse = ",")
    }
    key <- paste(
      if (is.null(scenario)) "NULL" else scenario,
      if (is.null(delta_sst)) "NULL" else format(delta_sst, digits = 8),
      if (is.null(target_year)) "NULL" else format(target_year, digits = 8),
      sensitivity_mode,
      k_beta,
      k_gamma,
      perturb_key,
      simulation_years,
      seed,
      sep = "|"
    )

    if (!exists(key, envir = cache, inherits = FALSE)) {
      climate_args <- list(
        sensitivity_mode = sensitivity_mode,
        k_beta = k_beta,
        k_gamma = k_gamma,
        perturb = perturb
      )
      if (!is.null(scenario)) climate_args$scenario <- scenario
      if (!is.null(delta_sst)) climate_args$delta_sst <- delta_sst
      if (!is.null(target_year)) climate_args$target_year <- target_year
      cfg <- make_hazard_cfg(
        simulation_years = simulation_years,
        climate = do.call(make_climate_cfg, climate_args)
      )
      out <- suppressWarnings(
        run_hazard_model(cfg, .hazard_test_targets, seed = seed, verbose = FALSE)
      )
      assign(key, out, envir = cache)
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

.summarize_hazard_reference <- function(out) {
  list(
    sim_summary = out$sim |>
      dplyr::group_by(location) |>
      dplyr::summarise(
        mean_total = mean(n_total),
        mean_ts = mean(n_ts),
        mean_hur = mean(n_hur),
        total_events = sum(n_total),
        .groups = "drop"
      ),
    rate_summary = out$rates |>
      dplyr::select(location, storm_class, lambda) |>
      dplyr::arrange(location, storm_class),
    metadata = out$run_metadata
  )
}

test_that("make_climate_cfg validates inputs and prints fields", {
  cfg <- make_climate_cfg(
    scenario = "ssp245",
    target_year = 2050,
    sensitivity_mode = "linear_shifted",
    k_beta = 0.2,
    k_gamma = 0.1,
    perturb = list(v_scale = 0.05)
  )

  expect_s3_class(cfg, "climate_cfg")
  expect_equal(cfg$scenario, "ssp245")
  expect_identical(cfg$input_mode, "scenario_helper")
  expect_equal(cfg$target_year, 2050)
  expect_equal(cfg$sensitivity_mode, "linear_shifted")
  expect_equal(cfg$k_beta, 0.2)
  expect_equal(cfg$k_gamma, 0.1)
  expect_equal(cfg[["perturb"]]$v_scale, 0.05)
  expect_equal(names(cfg[["perturb"]]), "v_scale")

  resolved <- resolve_climate_inputs(
    cfg,
    annual_counts = .make_climate_counts_fixture(),
    verbose = FALSE
  )
  expect_equal(resolved$perturb$v_scale, 0.05)
  expect_equal(resolved$perturb$r_scale, default_cc_params()$r_scale)
  expect_equal(resolved$perturb$speed_scale, default_cc_params()$speed_scale)
  expect_equal(resolved$perturb$precip_scale, default_cc_params()$precip_scale)

  txt <- paste(capture.output(print(cfg)), collapse = "\n")
  expect_match(txt, "Climate configuration")
  expect_match(txt, "Input mode")
  expect_match(txt, "ssp245")
  expect_match(txt, "Sensitivity mode")
  expect_match(txt, "Storm perturbation")

  stationary <- make_climate_cfg(scenario = "stationary")
  expect_equal(stationary$scenario, "stationary")
  expect_null(stationary[["perturb"]])
  expect_identical(stationary$input_mode, "scenario_helper")

  direct <- make_climate_cfg(delta_sst = 0.45)
  expect_identical(direct$input_mode, "direct_delta_sst")
  expect_true(is.na(direct$scenario))
  expect_equal(direct$delta_sst, 0.45)

  expect_warning(
    expect_error(make_climate_cfg(k_beta = "abc"), "single finite numeric"),
    "NAs introduced by coercion"
  )
  expect_error(make_climate_cfg(sensitivity_mode = "bad"))
  expect_error(make_climate_cfg(k_gamma = NA_real_), "single finite numeric")
  expect_error(
    make_climate_cfg(scenario = "ssp245", target_year = 2050, delta_sst = 0.1),
    "Conflicting climate specification"
  )
  expect_error(make_climate_cfg(scenario = "ssp245"), "Scenario helper mode requires target_year")
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
  climate <- make_climate_cfg(scenario = "ssp245", target_year = 2050, sensitivity_mode = "fixed")

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    verbose = FALSE
  )

  expect_true(resolved$delta_sst > 0)
  expect_equal(resolved$beta_sst_raw, resolved$beta_0)
  expect_equal(resolved$gamma, resolved$gamma_0)
})

test_that("linear-shifted mode applies expected sensitivity transforms", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(
    scenario = "ssp245",
    target_year = 2050,
    sensitivity_mode = "linear_shifted",
    k_beta = 0.3,
    k_gamma = -0.4
  )

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
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

test_that("get_scenario_delta interpolates and clamps by target year", {
  info <- sst_scenario_info("all")
  ssp585 <- info[info$scenario == "ssp585", , drop = FALSE]

  expect_equal(
    get_scenario_delta("ssp585", target_year = 2050),
    ssp585$delta_sst_2050[[1]]
  )

  expected_2090 <- ssp585$delta_sst_2050[[1]] +
    (2090 - 2050) / (2100 - 2050) *
    (ssp585$delta_sst_2100[[1]] - ssp585$delta_sst_2050[[1]])
  expect_equal(
    get_scenario_delta("ssp585", target_year = 2090),
    expected_2090
  )

  expect_equal(
    get_scenario_delta("ssp585", target_year = 2115),
    ssp585$delta_sst_2100[[1]]
  )

  expect_error(get_scenario_delta("nonexistent", target_year = 2050), "Unknown SST scenario")
  expect_error(get_scenario_delta("ssp245"), "argument \"target_year\" is missing")
})

test_that("resolver uses scenario table delta_sst path", {
  annual_counts <- .make_climate_counts_fixture()
  climate <- make_climate_cfg(scenario = "ssp585", target_year = 2090)

  resolved <- resolve_climate_inputs(
    climate,
    annual_counts = annual_counts,
    verbose = FALSE
  )

  expect_equal(
    resolved$delta_sst,
    get_scenario_delta("ssp585", target_year = 2090)
  )
  expect_identical(resolved$input_mode, "scenario_helper")
  expect_identical(resolved$climate_mode, "future")
})

test_that("scenario helper and direct delta_sst resolve identical downstream climate adjustments", {
  annual_counts <- .make_climate_counts_fixture()
  lambda_table <- tibble::tibble(
    storm_class = c("TS", "HUR"),
    lambda = c(0.8, 0.2)
  )
  helper <- resolve_climate_inputs(
    make_climate_cfg(scenario = "ssp245", target_year = 2050, sensitivity_mode = "linear_shifted", k_beta = 0.2, k_gamma = 0.1),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )
  direct <- resolve_climate_inputs(
    make_climate_cfg(delta_sst = helper$delta_sst, sensitivity_mode = "linear_shifted", k_beta = 0.2, k_gamma = 0.1),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )

  expect_identical(direct$input_mode, "direct_delta_sst")
  expect_equal(direct$delta_sst, helper$delta_sst)
  expect_equal(direct$beta_sst, helper$beta_sst)
  expect_equal(direct$beta_sst_raw, helper$beta_sst_raw)
  expect_equal(direct$gamma, helper$gamma)
  expect_equal(direct$rate_scale, helper$rate_scale)
  expect_equal(direct$raw_rate_scale, helper$raw_rate_scale)
  expect_equal(direct$response_regime$redistribution_strength, helper$response_regime$redistribution_strength)
  expect_equal(direct$response_regime$redistribution_bounds, helper$response_regime$redistribution_bounds)
})

test_that("estimate_beta_sst rejects target-conditioned annual counts", {
  annual_counts <- .make_climate_counts_fixture() |>
    dplyr::mutate(location = rep(c("A", "B"), length.out = dplyr::n()))
  sst_df <- get_mdr_sst_builtin() |>
    compute_sst_anomaly()

  expect_error(
    estimate_beta_sst(annual_counts, sst_df = sst_df, verbose = FALSE),
    "must not contain a location column"
  )
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
    n_years_sim = 2000,
    delta_sst = 0
  )

  set.seed(1)
  sim1 <- simulate_twolevel_counts(
    lambda_table = lambda_table,
    k_hat = 1e6,
    n_years_sim = 2000,
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

test_that("run_hazard_model uses embedded cfg climate and records seed metadata", {
  cfg_future <- make_hazard_cfg(
    simulation_years = 8L,
    climate = make_climate_cfg(scenario = "ssp245", target_year = 2050)
  )
  base <- .run_hazard_reference(
    scenario = "stationary",
    simulation_years = 8L,
    seed = 321L
  )

  fut1 <- .run_hazard_reference(
    scenario = "ssp245",
    target_year = 2050,
    simulation_years = 8L,
    seed = 321L
  )
  auto <- suppressWarnings(run_hazard_model(cfg_future, .hazard_test_targets, seed = NULL, verbose = FALSE))

  expect_identical(fut1$run_metadata$seed, 321L)
  expect_true(is.integer(auto$run_metadata$seed))
  expect_true(length(auto$run_metadata$seed) == 1L)
  expect_true(is.finite(auto$run_metadata$seed))
  expect_true(is.list(fut1$run_metadata$climate))
  expect_true(all(c(
    "input_mode", "target_year", "delta_sst", "beta_0", "beta_sst_raw", "beta_sst",
    "gamma_0", "gamma", "p_hur_base", "rate_scale"
  ) %in% names(fut1$run_metadata$climate)))

  expect_equal(attr(base$fit, "delta_sst"), 0)
  expect_equal(attr(base$fit, "rate_scale"), 1)
  expect_identical(attr(base$fit, "climate_mode"), "baseline")
  expect_identical(attr(base$fit, "climate_input_mode"), "scenario_helper")
  expect_gt(attr(fut1$fit, "delta_sst"), 0)
  expect_gt(attr(fut1$fit, "rate_scale"), 1)
  expect_identical(attr(fut1$fit, "climate_mode"), "future")
  expect_identical(attr(fut1$fit, "climate_input_mode"), "scenario_helper")
  expect_gt(mean(fut1$sim$n_total), mean(base$sim$n_total))
  expect_false(all(diff(fut1$sim$n_total) == 0))
})

test_that("stationary and default future helper runs remain deterministic", {
  expected_stationary <- list(
    sim_summary = tibble::tibble(
      location = "Saba",
      mean_total = 1,
      mean_ts = 0.875,
      mean_hur = 0.125,
      total_events = 8L
    ),
    rate_summary = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("HUR", "TS"),
      lambda = c(0.06382979, 0.48936170)
    ),
    metadata = list(
      seed = 222L,
      parameter_id = "params-00097b8d",
      lambda_scaling_mode = "target"
    )
  )
  expected_future <- list(
    sim_summary = tibble::tibble(
      location = "Saba",
      mean_total = 1,
      mean_ts = 0.875,
      mean_hur = 0.125,
      total_events = 8L
    ),
    rate_summary = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("HUR", "TS"),
      lambda = c(0.06382979, 0.48936170)
    ),
    metadata = list(
      seed = 222L,
      parameter_id = "params-00097b8d",
      lambda_scaling_mode = "target"
    )
  )

  out_stationary <- .run_hazard_reference(scenario = "stationary")
  out_future <- .run_hazard_reference(scenario = "ssp245", target_year = 2050)
  stationary_summary <- .summarize_hazard_reference(out_stationary)
  future_summary <- .summarize_hazard_reference(out_future)

  expect_equal(stationary_summary$sim_summary, expected_stationary$sim_summary, tolerance = 1e-8)
  expect_equal(stationary_summary$rate_summary, expected_stationary$rate_summary, tolerance = 1e-8)
  expect_identical(stationary_summary$metadata$seed, expected_stationary$metadata$seed)
  expect_identical(stationary_summary$metadata$parameter_id, expected_stationary$metadata$parameter_id)
  expect_identical(stationary_summary$metadata$lambda_scaling_mode, expected_stationary$metadata$lambda_scaling_mode)
  expect_equal(attr(out_stationary$fit, "delta_sst"), 0)
  expect_identical(attr(out_stationary$fit, "climate_mode"), "baseline")
  expect_identical(attr(out_stationary$fit, "climate_scenario"), "stationary")

  expect_equal(future_summary$sim_summary, expected_future$sim_summary, tolerance = 1e-8)
  expect_equal(future_summary$rate_summary, expected_future$rate_summary, tolerance = 1e-8)
  expect_identical(future_summary$metadata$seed, expected_future$metadata$seed)
  expect_identical(future_summary$metadata$parameter_id, expected_future$metadata$parameter_id)
  expect_identical(future_summary$metadata$lambda_scaling_mode, expected_future$metadata$lambda_scaling_mode)
  expect_equal(attr(out_future$fit, "delta_sst"), 0.50, tolerance = 1e-8)
  expect_identical(attr(out_future$fit, "climate_mode"), "future")
  expect_identical(attr(out_future$fit, "climate_scenario"), "ssp245")
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

test_that("run_hazard_model rejects missing or legacy climate config state", {
  cfg <- make_hazard_cfg(
    data_path = "definitely-missing.csv",
    simulation_years = 12L,
    climate = make_climate_cfg(scenario = "stationary")
  )
  legacy_shift <- structure(list(delta_sst = 0.8), class = c("climate_shift", "list"))
  legacy_response <- structure(list(beta_sst = 0.5, gamma = 0.1), class = c("climate_response", "list"))
  legacy_input <- structure(
    list(shift = legacy_shift, response = legacy_response),
    class = c("climate_input", "list")
  )

  expect_error(
    make_hazard_cfg(simulation_years = 40L, climate = legacy_shift),
    "climate must be created by make_climate_cfg"
  )
  expect_error(
    make_hazard_cfg(simulation_years = 40L, climate = legacy_response),
    "climate must be created by make_climate_cfg"
  )
  expect_error(
    make_hazard_cfg(simulation_years = 40L, climate = legacy_input),
    "climate must be created by make_climate_cfg"
  )

  cfg_missing <- unclass(cfg)
  class(cfg_missing) <- c("hazard_cfg", "list")
  cfg_missing$climate <- NULL
  expect_error(
    suppressWarnings(run_hazard_model(cfg_missing, .hazard_test_targets, seed = 11L, verbose = FALSE)),
    "cfg\\$climate is required"
  )

  cfg_legacy <- unclass(cfg)
  class(cfg_legacy) <- c("hazard_cfg", "list")
  cfg_legacy$climate <- structure(
    list(enabled = FALSE, scenario = "stationary", sst_source = "builtin", baseline_years = 1991L:2020L, start_year = 2025L, sensitivity_mode = "fixed", k_beta = 0, k_gamma = 0, perturb = NULL),
    class = c("climate_cfg", "list")
  )
  expect_error(
    suppressWarnings(run_hazard_model(cfg_legacy, .hazard_test_targets, seed = 11L, verbose = FALSE)),
    "Climate-off mode has been removed"
  )
})

test_that("returned climate metadata uses only simplified resolved schema", {
  out <- .run_hazard_reference(
    scenario = "ssp245",
    target_year = 2050,
    sensitivity_mode = "linear_shifted",
    k_beta = 0.2,
    k_gamma = 0.1,
    perturb = list(),
    simulation_years = 8L,
    seed = 11L
  )

  expect_type(out$cfg$climate, "list")
  expect_false(inherits(out$cfg$climate, "climate_input"))
  expect_false("shift" %in% names(out$cfg$climate))
  expect_false("response" %in% names(out$cfg$climate))
  expect_true(all(c(
    "scenario", "source", "input_mode", "delta_sst", "baseline_years",
    "target_year", "sensitivity_mode", "k_beta", "k_gamma",
    "beta_0", "gamma_0", "beta_sst_raw", "beta_sst", "gamma", "p_hur_base", "rate_scale",
    "annual_count_series", "annual_count_source",
    "beta_guardrail", "sst_scale_guardrail",
    "climate_mode",
    "perturb", "perturb_state", "config"
  ) %in% names(out$cfg$climate)))
  expect_s3_class(out$cfg$climate$config, "climate_cfg")
  expect_identical(attr(out$fit, "climate_cfg"), out$cfg$climate$config)
})

test_that("run_hazard_model produces equal climate-adjusted outputs for equal delta_sst", {
  helper <- .run_hazard_reference(
    scenario = "ssp245",
    target_year = 2050,
    simulation_years = 12L,
    seed = 444L
  )
  direct <- .run_hazard_reference(
    scenario = NULL,
    delta_sst = helper$cfg$climate$delta_sst,
    simulation_years = 12L,
    seed = 444L
  )

  expect_identical(helper$run_metadata$climate$input_mode, "scenario_helper")
  expect_identical(direct$run_metadata$climate$input_mode, "direct_delta_sst")
  expect_equal(helper$cfg$climate$delta_sst, direct$cfg$climate$delta_sst)
  expect_equal(helper$cfg$climate$beta_sst, direct$cfg$climate$beta_sst)
  expect_equal(helper$cfg$climate$gamma, direct$cfg$climate$gamma)
  expect_equal(helper$cfg$climate$rate_scale, direct$cfg$climate$rate_scale)
  expect_equal(helper$run_metadata$parameter_hash, direct$run_metadata$parameter_hash)
  expect_equal(helper$sim$n_total, direct$sim$n_total)
  expect_equal(helper$sim$n_ts, direct$sim$n_ts)
  expect_equal(helper$sim$n_hur, direct$sim$n_hur)
})

test_that("climate calibration is independent of duplicated target geometry", {
  cfg_future <- make_hazard_cfg(
    simulation_years = 12L,
    climate = make_climate_cfg(scenario = "ssp585", target_year = 2050)
  )

  out_single <- suppressWarnings(
    run_hazard_model(cfg_future, .hazard_test_targets, seed = 777L, verbose = FALSE)
  )
  out_dup <- suppressWarnings(
    run_hazard_model(cfg_future, .hazard_test_targets_dup, seed = 777L, verbose = FALSE)
  )

  expect_equal(out_single$cfg$climate$beta_0, out_dup$cfg$climate$beta_0, tolerance = 1e-8)
  expect_equal(out_single$cfg$climate$beta_sst, out_dup$cfg$climate$beta_sst, tolerance = 1e-8)
  expect_equal(out_single$cfg$climate$rate_scale, out_dup$cfg$climate$rate_scale, tolerance = 1e-8)
  expect_identical(out_dup$cfg$climate$annual_count_source, "basin_unique_storm_year_counts")
})

test_that("future climate multiplier stays within documented sanity guardrail", {
  out <- .run_hazard_reference(
    scenario = "ssp585",
    target_year = 2050,
    simulation_years = 12L,
    seed = 999L
  )

  expect_gte(out$cfg$climate$rate_scale, 1)
  expect_lte(out$cfg$climate$rate_scale, 4)
  expect_true(is.data.frame(out$cfg$climate$annual_count_series))
  expect_true(all(c("year", "N") %in% names(out$cfg$climate$annual_count_series)))
})

test_that("baseline and future climate runs share resolved historical sensitivities", {
  annual_counts <- .make_climate_counts_fixture()
  lambda_table <- tibble::tibble(
    storm_class = c("TS", "HUR"),
    lambda = c(0.8, 0.2)
  )
  baseline <- resolve_climate_inputs(
    make_climate_cfg(scenario = "stationary", sensitivity_mode = "fixed"),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )
  future <- resolve_climate_inputs(
    make_climate_cfg(scenario = "ssp245", target_year = 2050, sensitivity_mode = "fixed"),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )

  expect_identical(baseline$climate_mode, "baseline")
  expect_identical(future$climate_mode, "future")
  expect_identical(baseline$beta_0, future$beta_0)
  expect_identical(baseline$gamma_0, future$gamma_0)
  expect_identical(baseline$beta_sst_raw, future$beta_sst_raw)
  expect_identical(baseline$gamma, future$gamma)
  expect_equal(baseline$delta_sst, 0)
  expect_gt(future$delta_sst, 0)
})

test_that("neutral shifted mode clamps fixed-mode negative sensitivities to zero", {
  annual_counts <- .make_climate_counts_fixture()
  lambda_table <- tibble::tibble(
    storm_class = c("TS", "HUR"),
    lambda = c(0.8, 0.2)
  )
  fixed <- resolve_climate_inputs(
    make_climate_cfg(scenario = "ssp245", target_year = 2050, sensitivity_mode = "fixed"),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )
  shift_neutral <- resolve_climate_inputs(
    make_climate_cfg(
      scenario = "ssp245",
      target_year = 2050,
      sensitivity_mode = "linear_shifted",
      k_beta = 0,
      k_gamma = 0
    ),
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    verbose = FALSE
  )

  set.seed(222)
  sim_fixed <- simulate_twolevel_counts(
    lambda_table = lambda_table,
    k_hat = 1e6,
    n_years_sim = 64,
    delta_sst = fixed$delta_sst,
    beta_sst = fixed$beta_sst,
    gamma_intensity = fixed$gamma,
    p_hur_base = fixed$p_hur_base
  )
  set.seed(222)
  sim_shift_neutral <- simulate_twolevel_counts(
    lambda_table = lambda_table,
    k_hat = 1e6,
    n_years_sim = 64,
    delta_sst = shift_neutral$delta_sst,
    beta_sst = shift_neutral$beta_sst,
    gamma_intensity = shift_neutral$gamma,
    p_hur_base = shift_neutral$p_hur_base
  )

  expect_lt(fixed$beta_sst, 0)
  expect_equal(shift_neutral$beta_sst, 0)
  expect_equal(shift_neutral$gamma, 0)
  expect_true(all(sim_fixed$climate_scale < sim_shift_neutral$climate_scale))
})

test_that("transient SST workflow functions are deleted from the public surface", {
  expect_false(exists("make_sst_cfg", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false(exists("generate_sst_scenario", envir = asNamespace("ipdcstorm"), inherits = FALSE))
  expect_false("climate" %in% names(formals(run_hazard_model)))
  expect_false("sst_cfg" %in% names(formals(run_hazard_model)))
})
