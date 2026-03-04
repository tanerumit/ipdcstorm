test_that("lambda scalers handle clamping and missing references deterministically", {
  rate_tbl <- tibble::tibble(
    location = c("Saba", "Saba", "Miami"),
    storm_class = c("TS34plus", "HUR", "TS34plus"),
    lambda_model_raw = c(2.0, 0.5, 1.2),
    lambda_ref = c(0.4, NA_real_, 10.0),
    expected_ratio = c(0.55, 0.30, 0.75)
  )

  scalers <- ipdcstorm:::.lambda_scalers_from_rate_check(rate_tbl)

  saba_ts <- dplyr::filter(scalers, location == "Saba", storm_class == "TS34plus")
  saba_hur <- dplyr::filter(scalers, location == "Saba", storm_class == "HUR")
  miami_ts <- dplyr::filter(scalers, location == "Miami", storm_class == "TS34plus")

  expect_equal(saba_ts$lambda_target, 0.22, tolerance = 1e-12)
  expect_equal(saba_ts$lambda_scale, 0.25, tolerance = 1e-12)
  expect_identical(saba_ts$scale_status, "clamped_low")
  expect_true(saba_ts$scale_clamped)

  expect_equal(saba_hur$lambda_scale, 1, tolerance = 1e-12)
  expect_identical(saba_hur$scale_status, "no_ref")
  expect_equal(saba_hur$lambda_adj, saba_hur$lambda_model_raw, tolerance = 1e-12)

  expect_equal(miami_ts$lambda_target, 7.5, tolerance = 1e-12)
  expect_equal(miami_ts$lambda_scale, 4, tolerance = 1e-12)
  expect_identical(miami_ts$scale_status, "clamped_high")
  expect_true(miami_ts$scale_clamped)
})

test_that("adjusted lambdas are used by the simulation count pipeline", {
  lambda_table <- tibble::tibble(
    storm_class = c("TS", "HUR"),
    lambda = c(0.8, 0.2)
  )
  lambda_scalers <- tibble::tibble(
    location = c("Saba", "Saba"),
    storm_class = c("TS34plus", "HUR"),
    lambda_model_raw = c(1.0, 0.2),
    lambda_ref = c(1.5, 0.1),
    expected_ratio = c(1, 1),
    lambda_target = c(1.5, 0.1),
    lambda_scale = c(1.5, 0.5),
    lambda_adj = c(1.5, 0.1),
    scale_status = c("ok", "ok"),
    scale_clamped = c(FALSE, FALSE)
  )

  lambda_adj <- ipdcstorm:::.apply_lambda_scalers_to_lambda_table(
    lambda_table = lambda_table,
    location = "Saba",
    lambda_scalers = lambda_scalers
  )

  expect_equal(
    dplyr::pull(dplyr::filter(lambda_adj, storm_class == "TS"), lambda),
    1.4,
    tolerance = 1e-12
  )
  expect_equal(
    dplyr::pull(dplyr::filter(lambda_adj, storm_class == "HUR"), lambda),
    0.1,
    tolerance = 1e-12
  )

  set.seed(1)
  sim_raw <- simulate_twolevel_counts(lambda_table, k_hat = 1e6, n_years_sim = 20000)
  set.seed(1)
  sim_adj <- simulate_twolevel_counts(lambda_adj, k_hat = 1e6, n_years_sim = 20000)

  expect_gt(mean(sim_adj$n_total), mean(sim_raw$n_total) * 1.35)
  expect_lt(mean(sim_adj$p_hurricane), mean(sim_raw$p_hurricane))
})
