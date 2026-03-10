test_that("package loads and basic config constructors return lists", {
  cfg <- make_hazard_cfg()
  expect_true(is.list(cfg))

  sst_cfg <- make_sst_cfg(enabled = FALSE)
  expect_true(is.list(sst_cfg))
})

test_that("gamma intensity is constrained to non-negative values", {
  expect_error(
    make_sst_cfg(
      enabled = TRUE,
      scenario = "ssp245",
      advanced = list(gamma_intensity = -0.1)
    ),
    "advanced\\$gamma_intensity must be >= 0"
  )

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
  cfg <- make_hazard_cfg(simulation_years = 3L)
  targets <- tibble::tibble(
    location = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  expect_no_error(
    out <- run_hazard_model(cfg, targets, verbose = FALSE)
  )
  expect_true(all(c("sim", "events", "rates", "fit", "cfg") %in% names(out)))
  expect_true(nrow(out$sim) > 0)
  expect_true(all(out$sim$location == "Saba"))
  expect_match(out$cfg$data_path, "ibtracs\\.NA\\.list\\.v04r01\\.csv$")
})

test_that("run_hazard_model is deterministic for an explicit seed", {
  cfg <- make_hazard_cfg(simulation_years = 20L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  out1 <- run_hazard_model(cfg, targets, seed = 77L, verbose = FALSE)
  out2 <- run_hazard_model(cfg, targets, seed = 77L, verbose = FALSE)
  out3 <- run_hazard_model(cfg, targets, seed = 78L, verbose = FALSE)
  out_null <- run_hazard_model(cfg, targets, seed = NULL, verbose = FALSE)

  expect_identical(out1$sim, out2$sim)
  expect_identical(out1$run_metadata$seed, 77L)
  expect_false(identical(out1$sim, out3$sim))
  expect_true(is.null(out_null$run_metadata$seed))
  expect_true(nrow(out_null$sim) > 0)
})

test_that("run_hazard_model console output is structured and user-facing", {
  cfg <- make_hazard_cfg(simulation_years = 3L)
  targets <- tibble::tibble(
    name = "Saba",
    lat = 17.63,
    lon = -63.23
  )

  msg <- paste(
    capture.output(
      out <- run_hazard_model(cfg, targets, seed = 99L, verbose = TRUE),
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
  expect_match(msg, "Simulation")
  expect_match(msg, "Run metadata")
  expect_match(msg, "seed=99")
  expect_match(msg, "lambda_mode=")
  expect_false(grepl("seed=NA", msg, fixed = TRUE))
  expect_false(grepl("params=", msg, fixed = TRUE))
  expect_false(grepl("data=", msg, fixed = TRUE))
  expect_true(is.list(out$run_metadata))
})
