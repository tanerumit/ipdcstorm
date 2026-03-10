downscale_track_fixture <- function() {
  tibble::tibble(
    SID = c("storm_ts", "storm_ts", "storm_hur", "storm_hur"),
    iso_time = as.POSIXct(
      c("2000-08-01 00:00:00", "2000-08-01 06:00:00", "2000-09-10 00:00:00", "2000-09-10 06:00:00"),
      tz = "UTC"
    )
  )
}

downscale_event_fixture <- function() {
  tibble::tibble(
    location = c("Saba", "Saba"),
    storm_id = c("storm_ts", "storm_hur"),
    peak_wind_kt = c(45, 85),
    storm_intensity_kt = c(55, 110),
    min_pressure_hpa = c(995, 955),
    pressure_deficit_hpa = c(15, 55),
    rmw_mean_km = c(40, 24),
    start_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-09-10 00:00:00"), tz = "UTC"),
    end_time = as.POSIXct(c("2000-08-02 00:00:00", "2000-09-12 00:00:00"), tz = "UTC"),
    n_points = c(2L, 2L)
  )
}

downscale_out_fixture <- function() {
  list(
    cfg = list(resampling_method = "stratified"),
    trackpoints = list(Saba = downscale_track_fixture()),
    events = downscale_event_fixture(),
    sim = tibble::tibble(
      location = "Saba",
      sim_year = 1L,
      n_ts = 1L,
      n_hur = 1L
    )
  )
}

test_that("downscale library and samplers expose expected public behavior", {
  lib <- build_event_library(
    track_df = downscale_track_fixture(),
    event_df = downscale_event_fixture(),
    storm_classes = c("TD", "TS", "HUR"),
    seed = 10
  )
  lib_from_out <- build_event_library_from_out(
    out = downscale_out_fixture(),
    location = "Saba",
    seed = 10
  )

  sampled <- sample_events_for_year_extended(
    lib = lib,
    year = 2001,
    n_ts = 1L,
    n_hur = 1L,
    seed = 11
  )

  expect_true(is.list(lib))
  expect_true(is.function(lib$sample_doy))
  expect_true(is.function(lib$sample_event))
  expect_true(nrow(lib$events) >= 2)
  expect_identical(levels(lib$events$severity), c("TD", "TS", "HUR"))
  expect_equal(nrow(lib_from_out$events), nrow(lib$events))
  expect_equal(sort(sampled$event_class), c("HUR", "TS"))
  expect_true(all(sampled$dur_days >= 1L))
})

test_that("daily hazard helpers classify events and aggregate daily diagnostics", {
  sampled_events <- tibble::tibble(
    severity = c("TS", "HUR"),
    start_date = as.Date(c("2001-06-01", "2001-09-01")),
    dur_days = c(2L, 3L),
    V_peak = c(45, 90),
    event_id = c("evt_ts", "evt_hur"),
    event_class = c("TS", "HUR"),
    Pc_min_hPa = c(995, 955),
    dP_max_hPa = c(18, 55),
    RMW_mean_km = c(35, 22)
  )
  daily <- generate_daily_year_extended(2001, sampled_events, pulse_shape = "triangle")
  daily_aug <- dplyr::mutate(
    daily,
    location = "Saba",
    sim_year = 1L,
    scenario = "baseline",
    wind_gust_kt = .data$wind_kt * 1.2,
    surge_m = dplyr::if_else(is.finite(.data$pressure_hpa), 0.14 * (1013 - .data$pressure_hpa), NA_real_)
  )
  disrupted <- disruption_flags(daily_aug, thr_port = 40, thr_infra = 70, thr_surge = 5)
  peaks <- peak_wind_by_year(daily_aug)

  expect_equal(event_pulse(3, 60, shape = "triangle"), c(30, 60, 30))
  expect_true(any(is_tc_day(daily)))
  expect_true(any(is_hur_day(daily)))
  expect_true(any(exposure_hours(daily_aug, threshold_kt = 40) == 24))
  expect_true(all(c("port_disrupt", "infra_disrupt", "surge_disrupt") %in% names(disrupted)))
  expect_equal(peaks$peak_wind_kt, max(daily_aug$wind_kt))
})

test_that("daily hazard generation and damage helpers return public output schema", {
  out <- downscale_out_fixture()
  daily_list <- generate_daily_hazard_impact(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    gust_factor = 1.2,
    seed = 20
  )
  daily <- daily_list$Saba
  damage_aug <- add_damage_forcing(tibble::tibble(wind_kt = c(20, 60, 140)), dmax = 0.05)
  powerlaw <- damage_rate_from_wind(c(20, 80, 200), d_max = 0.08)

  expect_equal(names(daily_list), "Saba")
  expect_true(all(c(
    "location", "sim_year", "scenario", "date", "wind_kt", "wind_gust_kt",
    "surge_m", "event_id", "event_class", "pressure_hpa", "pressure_deficit_hpa",
    "rmw_km", "damage_intensity", "damage_rate", "cum_damage"
  ) %in% names(daily)))
  expect_identical(attr(daily, "gust_factor"), 1.2)
  expect_true(any(daily$wind_kt > 0))
  expect_equal(damage_aug$damage_intensity[c(1, 3)], c(0, 1))
  expect_true(all(powerlaw >= 0 & powerlaw <= 0.08, na.rm = TRUE))
})
