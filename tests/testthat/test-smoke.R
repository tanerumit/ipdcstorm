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
