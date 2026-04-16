downscale_track_fixture <- function() {
  tibble::tibble(
    SID = c(
      "storm_ts", "storm_ts",
      "storm_ts2", "storm_ts2",
      "storm_hur", "storm_hur",
      "storm_hur2", "storm_hur2"
    ),
    iso_time = as.POSIXct(
      c(
        "2000-08-01 00:00:00", "2000-08-01 06:00:00",
        "2000-07-15 00:00:00", "2000-07-15 06:00:00",
        "2000-09-10 00:00:00", "2000-09-10 06:00:00",
        "2000-10-05 00:00:00", "2000-10-05 06:00:00"
      ),
      tz = "UTC"
    )
  )
}

downscale_event_fixture <- function() {
  tibble::tibble(
    location = c("Saba", "Saba", "Saba", "Saba"),
    storm_id = c("storm_ts", "storm_ts2", "storm_hur", "storm_hur2"),
    peak_wind_kt = c(45, 52, 85, 95),
    storm_intensity_kt = c(55, 62, 110, 120),
    min_pressure_hpa = c(995, 990, 955, 945),
    pressure_deficit_hpa = c(15, 22, 55, 65),
    rmw_mean_km = c(40, 35, 24, 28),
    start_time = as.POSIXct(
      c("2000-08-01 00:00:00", "2000-07-15 00:00:00",
        "2000-09-10 00:00:00", "2000-10-05 00:00:00"),
      tz = "UTC"
    ),
    end_time = as.POSIXct(
      c("2000-08-02 00:00:00", "2000-07-17 00:00:00",
        "2000-09-12 00:00:00", "2000-10-08 00:00:00"),
      tz = "UTC"
    ),
    n_points = c(2L, 2L, 2L, 2L)
  )
}

downscale_out_fixture <- function() {
  list(
    cfg = list(resampling_method = "stratified"),
    trackpoints = list(Saba = downscale_track_fixture()),
    events = downscale_event_fixture(),
    sim = tibble::tibble(
      location = rep("Saba", 5L),
      sim_year = 1L:5L,
      n_ts = c(1L, 1L, 1L, 1L, 1L),
      n_hur = c(1L, 1L, 1L, 1L, 1L)
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
  lib_from_out <- ipdcstorm:::build_event_library_from_out(
    out = downscale_out_fixture(),
    location = "Saba",
    seed = 10
  )

  sampled <- ipdcstorm:::sample_events_for_year_extended(
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
  daily <- ipdcstorm:::generate_daily_year_extended(2001, sampled_events, pulse_shape = "triangle")
  daily_aug <- dplyr::mutate(
    daily,
    location = "Saba",
    sim_year = 1L,
    scenario = "baseline",
    wind_gust_kt = .data$wind_kt * 1.2,
    surge_m = dplyr::if_else(is.finite(.data$pressure_hpa), 0.14 * (1013 - .data$pressure_hpa), NA_real_)
  )
  disrupted <- ipdcstorm:::disruption_flags(daily_aug, thr_port = 40, thr_infra = 70, thr_surge = 5)
  peaks <- ipdcstorm:::peak_wind_by_year(daily_aug)

  expect_equal(ipdcstorm:::event_pulse(3, 60, shape = "triangle"), c(30, 60, 30))
  expect_true(any(ipdcstorm:::is_tc_day(daily)))
  expect_true(any(ipdcstorm:::is_hur_day(daily)))
  expect_true(any(ipdcstorm:::exposure_hours(daily_aug, threshold_kt = 40) == 24))
  expect_true(all(c("port_disrupt", "infra_disrupt", "surge_disrupt") %in% names(disrupted)))
  expect_equal(peaks$peak_wind_kt, max(daily_aug$wind_kt))
})

test_that("daily hazard generation and damage helpers return public output schema", {
  out <- downscale_out_fixture()
  daily_list <- generate_daily_hazard_impact_spatial(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    gust_factor = 1.2,
    seed = 20
  )
  daily <- daily_list$Saba
  damage_aug <- ipdcstorm:::add_damage_forcing(tibble::tibble(wind_kt = c(20, 60, 140)), dmax = 0.05)
  powerlaw <- ipdcstorm:::damage_rate_from_wind(c(20, 80, 200), d_max = 0.08)

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

test_that("daily hazard generation accepts damage defaults for both methods", {
  out <- downscale_out_fixture()

  intensity <- generate_daily_hazard_impact_spatial(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    damage = list(method = "intensity"),
    seed = 20
  )$Saba
  powerlaw <- generate_daily_hazard_impact_spatial(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    damage = list(method = "powerlaw"),
    seed = 20
  )$Saba

  expect_equal(
    intensity$damage_rate,
    ipdcstorm:::add_damage_forcing(tibble::tibble(wind_kt = intensity$wind_kt))$damage_rate
  )
  expect_equal(
    powerlaw$damage_rate,
    ipdcstorm:::damage_rate_from_wind(powerlaw$wind_kt)
  )
  expect_equal(
    powerlaw$damage_intensity,
    ipdcstorm:::add_damage_forcing(
      tibble::tibble(wind_kt = powerlaw$wind_kt),
      p = 3
    )$damage_intensity
  )
})

test_that("daily hazard generation honors damage overrides numerically", {
  out <- downscale_out_fixture()

  intensity <- generate_daily_hazard_impact_spatial(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    damage = list(method = "intensity", V0 = 40, V1 = 100, p = 2, dmax = 0.05),
    seed = 20
  )$Saba
  powerlaw <- generate_daily_hazard_impact_spatial(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    damage = list(method = "powerlaw", thr = 30, V_ref = 70, d_ref = 0.02, p = 2, d_max = 0.08),
    seed = 20
  )$Saba

  expected_intensity <- ipdcstorm:::add_damage_forcing(
    tibble::tibble(wind_kt = intensity$wind_kt),
    V0 = 40,
    V1 = 100,
    p = 2,
    dmax = 0.05
  )
  expect_equal(intensity$damage_intensity, expected_intensity$damage_intensity)
  expect_equal(intensity$damage_rate, expected_intensity$damage_rate)
  expect_equal(
    powerlaw$damage_rate,
    ipdcstorm:::damage_rate_from_wind(
      powerlaw$wind_kt,
      thr = 30,
      V_ref = 70,
      d_ref = 0.02,
      p = 2,
      d_max = 0.08
    )
  )
  expect_equal(
    powerlaw$damage_intensity,
    ipdcstorm:::add_damage_forcing(
      tibble::tibble(wind_kt = powerlaw$wind_kt),
      V0 = 34,
      V1 = 120,
      p = 2
    )$damage_intensity
  )
})

test_that("daily hazard generation rejects invalid damage specifications", {
  out <- downscale_out_fixture()
  base_args <- list(
    out = out,
    location = "Saba",
    sim_years = 1L,
    year0 = 2001L,
    seed = 20
  )

  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = "intensity"))),
    "`damage` must be a named list.",
    fixed = TRUE
  )
  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = list()))),
    "`damage$method` must be a single non-empty character string.",
    fixed = TRUE
  )
  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = list(method = "bad")))),
    "`damage$method` must be one of: intensity, powerlaw.",
    fixed = TRUE
  )
  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = list(method = "intensity", thr = 34)))),
    "unsupported field(s) for method 'intensity': thr",
    fixed = TRUE
  )
  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = list(method = "powerlaw", V0 = 34)))),
    "unsupported field(s) for method 'powerlaw': V0",
    fixed = TRUE
  )
  expect_error(
    do.call(generate_daily_hazard_impact_spatial, c(base_args, list(damage = list(method = "intensity", p = c(2, 3))))),
    "`damage$p` must be a single finite numeric value.",
    fixed = TRUE
  )
})

# =============================================================================
# Reproducibility tests
# =============================================================================

test_that("generate_daily_hazard_impact_spatial is reproducible with an explicit seed", {
  out <- downscale_out_fixture()

  run1 <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L, seed = 42L
  )
  run2 <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L, seed = 42L
  )

  expect_identical(run1$Saba$wind_kt, run2$Saba$wind_kt)
  expect_identical(run1$Saba$event_id, run2$Saba$event_id)
  expect_identical(run1$Saba$damage_rate, run2$Saba$damage_rate)
})

test_that("generate_daily_hazard_impact_spatial produces different output for different seeds", {
  out <- downscale_out_fixture()

  run_a <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:5, year0 = 2001L, seed = 1L
  )
  run_b <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:5, year0 = 2001L, seed = 2L
  )

  # Wind series should differ (at least in some years) when seeds differ
  expect_false(identical(run_a$Saba$wind_kt, run_b$Saba$wind_kt))
})

test_that("generate_daily_hazard_impact_spatial inherits seed from out$run_metadata", {
  out <- downscale_out_fixture()

  # Attach a seed in run_metadata, matching what run_hazard_model() stores
  out$run_metadata <- list(seed = 77L)

  run_inherited <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L
    # seed = NULL (default) — should pick up out$run_metadata$seed
  )
  run_explicit <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L, seed = 77L
  )

  expect_identical(run_inherited$Saba$wind_kt, run_explicit$Saba$wind_kt)
  expect_identical(run_inherited$Saba$event_id, run_explicit$Saba$event_id)
})

test_that("generate_daily_hazard_impact_spatial falls back to seed 1 when run_metadata is absent", {
  out <- downscale_out_fixture()
  # No run_metadata — should fall back to seed 1L

  run_default <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L
  )
  run_seed1 <- generate_daily_hazard_impact_spatial(
    out = out, location = "Saba", sim_years = 1:3, year0 = 2001L, seed = 1L
  )

  expect_identical(run_default$Saba$wind_kt, run_seed1$Saba$wind_kt)
})

test_that("generate_daily_hazard_impact_spatial rejects invalid seed values", {
  out <- downscale_out_fixture()

  expect_error(
    generate_daily_hazard_impact_spatial(out, "Saba", seed = "abc"),
    "seed must be NULL or a single finite numeric value.",
    fixed = TRUE
  )
  expect_error(
    generate_daily_hazard_impact_spatial(out, "Saba", seed = c(1L, 2L)),
    "seed must be NULL or a single finite numeric value.",
    fixed = TRUE
  )
  expect_error(
    generate_daily_hazard_impact_spatial(out, "Saba", seed = Inf),
    "seed must be NULL or a single finite numeric value.",
    fixed = TRUE
  )
})
