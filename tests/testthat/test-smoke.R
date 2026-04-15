test_that("package loads and basic config constructors return lists", {
  cfg <- make_hazard_cfg()
  expect_true(is.list(cfg))
})

test_that("gamma intensity estimator is constrained to non-negative values", {
  years <- 1970:2030
  annual_counts <- dplyr::bind_rows(
    tibble::tibble(year = years, storm_class = "TS", n_events = rep(30L, length(years))),
    tibble::tibble(
      year = years,
      storm_class = "HUR",
      n_events = pmax(1L, as.integer(round(seq(30, 1, length.out = length(years)))))
    )
  )
  sst_df <- tibble::tibble(
    year = years,
    sst_anomaly = seq(-0.5, 1.5, length.out = length(years))
  )

  gamma_fit <- estimate_gamma_intensity(
    annual_counts = annual_counts,
    sst_df = sst_df,
    min_year = min(years),
    gamma_prior = 0.065,
    verbose = FALSE
  )

  expect_gte(gamma_fit$gamma, 0)
})

test_that("run_hazard_model resolves packaged IBTrACS data and accepts location targets", {
  cfg <- make_hazard_cfg(simulation_years = 1L)
  targets <- tibble::tibble(
    location = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  expect_no_error(
    out <- suppressWarnings(run_hazard_model(cfg, targets, verbose = FALSE))
  )
  expect_true(all(c("sim", "events", "rates", "fit", "cfg") %in% names(out)))
  expect_true(nrow(out$sim) > 0)
  expect_true(all(out$sim$location == "Saba"))
  expect_match(out$cfg$data_path, "ibtracs\\.NA\\.list\\.v04r01\\.csv$")
})

test_that("run_hazard_model console output is structured and user-facing", {
  cfg <- make_hazard_cfg(simulation_years = 1L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  msg <- paste(
    capture.output(
      out <- suppressWarnings(run_hazard_model(cfg, targets, seed = 99L, verbose = TRUE)),
      type = "message"
    ),
    collapse = "\n"
  )

  expect_match(msg, "Loading data")
  expect_match(msg, "Input file\\s+:")
  expect_match(msg, "Raw NA input\\s+:")
  expect_match(msg, "Model start year\\s+:")
  expect_match(msg, "Target/location filtering")
  expect_match(msg, "Model sample\\s+:")
  expect_match(msg, "Rate calibration")
  expect_match(msg, "location/class adjustments")
  expect_match(msg, "missing_ref")
  expect_match(msg, "Climate")
  expect_match(msg, "delta_sst")
  expect_match(msg, "Simulation")
  expect_match(msg, "Run metadata")
  expect_match(msg, "seed=99")
  expect_match(msg, "lambda_mode=")
  expect_false(grepl("seed=NA", msg, fixed = TRUE))
  expect_false(grepl("params=", msg, fixed = TRUE))
  expect_false(grepl("data=", msg, fixed = TRUE))
  expect_true(is.list(out$run_metadata))
})

test_that("run_hazard_model is reproducible with a fixed seed", {
  cfg <- make_hazard_cfg(simulation_years = 5L)
  targets <- tibble::tibble(location = "Saba", lat = 17.63, lon = -63.23)

  run1 <- suppressWarnings(run_hazard_model(cfg, targets, seed = 42L, verbose = FALSE))
  run2 <- suppressWarnings(run_hazard_model(cfg, targets, seed = 42L, verbose = FALSE))

  expect_identical(run1$sim$n_ts,  run2$sim$n_ts)
  expect_identical(run1$sim$n_hur, run2$sim$n_hur)
  expect_identical(run1$run_metadata$seed, 42L)
})

test_that("run_hazard_model seed flows into generate_daily_hazard_impact_spatial via run_metadata", {
  cfg <- make_hazard_cfg(simulation_years = 3L)
  targets <- tibble::tibble(location = "Saba", lat = 17.63, lon = -63.23)

  out <- suppressWarnings(run_hazard_model(cfg, targets, seed = 55L, verbose = FALSE))

  # NULL seed inherits from out$run_metadata$seed (55L)
  daily_inherited <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3
  )
  # Explicit seed = 55L must give identical results
  daily_explicit <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, seed = 55L
  )

  expect_identical(daily_inherited$Saba$wind_kt,   daily_explicit$Saba$wind_kt)
  expect_identical(daily_inherited$Saba$event_id,  daily_explicit$Saba$event_id)
})
