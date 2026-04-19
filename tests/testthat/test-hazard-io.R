daily_io_fixture <- function() {
  list(
    Saba = tibble::tibble(
      sim_year = c(1L, 1L),
      doy = c(1L, 2L),
      wind_kt = c(0, 42),
      surge_m = c(NA_real_, 0.8),
      event_id = c(NA_character_, "storm_a_y2000_1"),
      pressure_hpa = c(NA_real_, 990),
      pressure_deficit_hpa = c(NA_real_, 23),
      rmw_km = c(NA_real_, 35),
      damage_intensity = c(0, 0.2),
      damage_rate = c(0, 0.01),
      cum_damage = c(0, 0.01)
    ),
    St_Martin = tibble::tibble(
      sim_year = c(1L, 1L),
      doy = c(1L, 2L),
      wind_kt = c(5, 18),
      surge_m = c(0.1, 0.2),
      event_id = c(NA_character_, NA_character_),
      pressure_hpa = c(NA_real_, NA_real_),
      pressure_deficit_hpa = c(NA_real_, NA_real_),
      rmw_km = c(NA_real_, NA_real_),
      damage_intensity = c(0, 0),
      damage_rate = c(0, 0),
      cum_damage = c(0, 0)
    )
  )
}

test_that("save_daily_hazard_csvs writes one CSV per location with scenario and location in the filename", {
  out_dir <- file.path(withr::local_tempdir(), "output", "raw")

  paths <- save_daily_hazard_csvs(
    daily = daily_io_fixture(),
    scenario = "knmi_Ld",
    out_dir = out_dir
  )

  expect_setequal(names(paths), c("Saba", "St_Martin"))
  expect_true(dir.exists(out_dir))
  expect_match(basename(paths[["Saba"]]), "^daily_knmi_Ld__Saba\\.csv$")
  expect_match(basename(paths[["St_Martin"]]), "^daily_knmi_Ld__St_Martin\\.csv$")
  expect_true(all(file.exists(unname(paths))))

  saba_written <- utils::read.csv(paths[["Saba"]], stringsAsFactors = FALSE)
  expect_identical(
    names(saba_written),
    c(
      "sim_year", "doy", "wind_kt", "surge_m", "event_id",
      "pressure_hpa", "pressure_deficit_hpa", "rmw_km",
      "damage_intensity", "damage_rate", "cum_damage"
    )
  )
  expect_equal(nrow(saba_written), 2L)
})
