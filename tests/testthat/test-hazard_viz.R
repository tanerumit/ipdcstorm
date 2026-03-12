hazard_viz_daily_fixture <- function() {
  tibble::tibble(
    date = as.Date(c(
      "2001-08-01", "2001-08-02", "2001-08-03", "2001-08-04",
      "2002-08-01", "2002-08-02", "2002-08-03", "2002-08-04",
      "2001-09-01", "2001-09-02", "2001-09-03", "2001-09-04",
      "2002-09-01", "2002-09-02", "2002-09-03", "2002-09-04"
    )),
    wind_kt = c(
      20, 40, 55, 25,
      18, 45, 65, 22,
      15, 38, 60, 20,
      12, 42, 70, 18
    ),
    event_id = c(
      NA, 1, 1, NA,
      NA, 2, 2, NA,
      NA, 1, 1, NA,
      NA, 2, 2, NA
    ),
    event_class = c(
      NA, "TS", "TS", NA,
      NA, "HUR", "HUR", NA,
      NA, "TS", "TS", NA,
      NA, "HUR", "HUR", NA
    ),
    location = c(
      rep("St. Eustatius", 8),
      rep("Saba", 8)
    ),
    sim_year = c(
      2001L, 2001L, 2001L, 2001L,
      2002L, 2002L, 2002L, 2002L,
      2001L, 2001L, 2001L, 2001L,
      2002L, 2002L, 2002L, 2002L
    )
  )
}

test_that("save_hazard_viz_plots writes validation plots into output/validation", {
  tmp_dir <- withr::local_tempdir()
  withr::local_dir(tmp_dir)

  result <- save_hazard_viz_plots(
    daily = hazard_viz_daily_fixture(),
    output_dir = file.path("legacy", "path"),
    location_name = "ignored"
  )

  validation_dir <- file.path("output", "validation")
  expected_files <- file.path(
    validation_dir,
    c(
      "Saba_annual_event_count_probability.png",
      "Saba_intensity_duration.png",
      "Saba_monthly_events_boxplot.png",
      "Saba_wind_distribution.png",
      "Saba_wind_timeseries.png",
      "St_Eustatius_annual_event_count_probability.png",
      "St_Eustatius_intensity_duration.png",
      "St_Eustatius_monthly_events_boxplot.png",
      "St_Eustatius_wind_distribution.png",
      "St_Eustatius_wind_timeseries.png",
      "rate_comparison.png"
    )
  )

  expect_true(dir.exists(validation_dir))
  expect_setequal(unname(result$paths), expected_files)
  expect_true(all(file.exists(expected_files)))
  expect_true(file.exists(file.path(validation_dir, "rate_comparison.png")))
  expect_false(dir.exists(file.path(validation_dir, "comparison")))
  expect_length(list.dirs(validation_dir, recursive = FALSE, full.names = TRUE), 0L)
})

test_that("save_hazard_viz_plots uses deterministic validation filenames", {
  tmp_dir <- withr::local_tempdir()
  withr::local_dir(tmp_dir)

  first_result <- save_hazard_viz_plots(
    daily = hazard_viz_daily_fixture(),
    output_dir = "unused",
    location_name = "ignored"
  )
  second_result <- save_hazard_viz_plots(
    daily = hazard_viz_daily_fixture(),
    output_dir = file.path("another", "unused"),
    location_name = "still_ignored"
  )

  expect_identical(first_result$paths, second_result$paths)
})
