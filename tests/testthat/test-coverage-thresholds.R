coverage_daily_fixture <- function() {
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

test_that("read_ibtracs_clean validates inputs and reads cleaned-schema CSVs", {
  expect_error(ipdcstorm::read_ibtracs_clean(NA_character_), "File not found")
  expect_error(ipdcstorm::read_ibtracs_clean("does-not-exist.csv"), "File not found")

  tmp_csv <- tempfile(fileext = ".csv")
  readr::write_csv(
    tibble::tibble(
      SID = c("AL01", "AL01", "EP01"),
      SEASON = c(2020L, 2020L, 2021L),
      BASIN = c("na", "na", "ep"),
      iso_time = c(
        "2020-08-01 00:00:00",
        "2020-08-01 06:00:00",
        "2021-09-02 00:00:00"
      ),
      lat = c("18.0", "18.5", "12.0"),
      lon = c("-63.0", "-63.2", "-110.0"),
      wind_kt = c("40", "55", "60"),
      r34_ne_nm = c("60", "65", "70"),
      r34_se_nm = c("50", "55", "60"),
      r34_sw_nm = c("40", "45", "50"),
      r34_nw_nm = c("30", "35", "40"),
      pres_hpa = c("995", "", "980"),
      pres_source = c("USA", "", "WMO"),
      usa_pres_hpa = c("995", "", ""),
      wmo_pres_hpa = c("", "998", "980"),
      storm_name = c("ALPHA", "ALPHA", "BETA"),
      storm_status = c("TS", "TS", "HU")
    ),
    tmp_csv
  )

  out <- suppressMessages(
    ipdcstorm::read_ibtracs_clean(tmp_csv, basin = "NA", season = 2020L, verbose = TRUE)
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 2)
  expect_true(all(out$BASIN == "NA"))
  expect_true(inherits(out$iso_time, "POSIXct"))
  expect_equal(out$storm_name[[1]], "ALPHA")
  expect_true(all(c("pres_hpa", "pres_source", "storm_status") %in% names(out)))
})

test_that("read_ibtracs_clean reads legacy USA-field CSVs with fallbacks", {
  tmp_csv <- tempfile(fileext = ".csv")
  readr::write_csv(
    tibble::tibble(
      SID = c("AL01", "AL01", "AL02"),
      SEASON = c(2020L, 2020L, 2021L),
      BASIN = c("NA", "NA", "NA"),
      ISO_TIME = c(
        "2020-08-01 00:00:00",
        "2020-08-01 06:00:00",
        "2021-09-02 00:00:00"
      ),
      USA_LAT = c("18.0", "", "19.2"),
      USA_LON = c("-63.0", "", "-62.0"),
      LAT = c("18.0", "18.4", "19.2"),
      LON = c("-63.0", "-63.4", "-62.0"),
      USA_WIND = c("45", "55", "80"),
      USA_R34_NE = c("60", "65", "90"),
      USA_R34_SE = c("50", "55", "85"),
      USA_R34_SW = c("40", "45", "80"),
      USA_R34_NW = c("30", "35", "75"),
      USA_R50_NE = c("20", "25", "40"),
      USA_R50_SE = c("20", "25", "35"),
      USA_R50_SW = c("20", "25", "30"),
      USA_R50_NW = c("20", "25", "25"),
      USA_R64_NE = c("", "", "20"),
      USA_R64_SE = c("", "", "20"),
      USA_R64_SW = c("", "", "15"),
      USA_R64_NW = c("", "", "15"),
      USA_PRES = c("995", "", "970"),
      WMO_PRES = c("", "998", ""),
      USA_POCI = c("1010", "1008", "1005"),
      USA_ROCI = c("180", "170", "150"),
      USA_RMW = c("20", "999", "25"),
      USA_STATUS = c("TS", "TS", "HU"),
      NAME = c("ALPHA", "ALPHA", "BETA"),
      STORM_SPEED = c("12", "10", "14"),
      STORM_DIR = c("270", "275", "300")
    ),
    tmp_csv
  )

  out <- ipdcstorm::read_ibtracs_clean(
    tmp_csv,
    basin = "NA",
    season = c(2020L, 2021L),
    keep_all = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(out), 3)
  expect_true(any(out$lat == 18.4, na.rm = TRUE))
  expect_true(any(out$lon == -63.4, na.rm = TRUE))
  expect_true(all(c("USA", "WMO") %in% out$pres_source))
  expect_true(any(is.na(out$rmw_km)))
  expect_true(all(c("r50_ne_nm", "r64_ne_nm", "storm_dir_deg") %in% names(out)))

  expect_error(
    ipdcstorm::read_ibtracs_clean(tmp_csv, basin = "EP", season = 2020L, verbose = FALSE),
    "After filtering, 0 rows remain"
  )
})

test_that("hazard viz helpers and plot constructors cover major branches", {
  daily <- coverage_daily_fixture()

  daily_prepped <- ipdcstorm:::prep_daily(daily)
  expect_true(all(c("doy", "month", "year") %in% names(daily_prepped)))
  expect_equal(daily_prepped$doy[[1]], 213L)
  expect_equal(daily_prepped$month[[1]], 8L)
  expect_equal(daily_prepped$year[[1]], 2001L)

  events <- ipdcstorm:::prep_events(daily)
  expect_equal(sort(unique(events$event_class)), c("HUR", "TS"))
  expect_true(all(events$dur_days >= 2L))
  expect_true(all(c("start_doy", "start_month", "max_wind_kt") %in% names(events)))

  p_wind <- ipdcstorm::plot_wind_timeseries(daily, events = events, show_thresholds = FALSE)
  p_starts <- ipdcstorm::plot_seasonality_doy(daily, metric = "starts", facet_class = FALSE, binwidth = 14)
  p_days <- ipdcstorm::plot_seasonality_doy(daily, metric = "event_days", facet_class = TRUE, binwidth = 7)
  p_month <- ipdcstorm::plot_monthly_events(daily)
  p_doy_raw <- ipdcstorm::plot_doy_wind(daily, smooth = FALSE)
  p_doy_smooth <- ipdcstorm::plot_doy_wind(daily, smooth = TRUE, span = 0.3)
  p_quant_log <- ipdcstorm::plot_monthly_quantiles(daily, log_scale = TRUE)
  p_annual_events <- ipdcstorm::plot_annual_counts(daily, metric = "events", show_poisson = TRUE)
  p_annual_days <- ipdcstorm::plot_annual_counts(daily, metric = "days", show_poisson = FALSE)
  p_scatter <- ipdcstorm::plot_intensity_duration(events = events)
  p_density <- ipdcstorm::plot_wind_distribution(daily, type = "density", log_scale = TRUE)
  p_hist <- ipdcstorm::plot_wind_distribution(daily, type = "histogram", log_scale = FALSE)
  p_rate <- ipdcstorm:::.plot_rate_comparison(daily)
  p_return_block <- ipdcstorm::plot_return_levels(daily, block_maxima = TRUE)
  p_return_pot <- ipdcstorm::plot_return_levels(daily, block_maxima = FALSE, threshold = 40)

  for (plot_obj in list(
    p_wind, p_starts, p_days, p_month, p_doy_raw, p_doy_smooth,
    p_quant_log, p_annual_events, p_annual_days, p_scatter,
    p_density, p_hist, p_rate, p_return_block, p_return_pot
  )) {
    expect_s3_class(plot_obj, "ggplot")
  }

  expect_equal(length(p_wind$layers), 2L)
  expect_true(inherits(p_quant_log$scales$scales[[2]], "ScaleContinuousPosition"))
  expect_match(p_annual_events$labels$subtitle, "Poisson probabilities")
  expect_match(p_annual_days$labels$title, "Event-Days")
  expect_equal(ipdcstorm:::.hazard_viz_location_id(" St. Eustatius "), "St_Eustatius")

  expect_error(ipdcstorm::plot_intensity_duration(), "is not TRUE")
  expect_error(ipdcstorm::plot_return_levels(daily, block_maxima = FALSE), "threshold")
})

test_that("climate data helpers and scenario helpers handle edge cases", {
  tmp_csv <- tempfile(fileext = ".csv")
  readr::write_csv(
    tibble::tibble(year = c(2000L, 2001L, 2002L), sst_mdr_aso = c(26.5, 26.8, NA_real_)),
    tmp_csv
  )

  csv_out <- ipdcstorm::read_mdr_sst_csv(tmp_csv)
  expect_equal(csv_out$year, c(2000L, 2001L))
  expect_equal(csv_out$sst_mdr_aso, c(26.5, 26.8))

  bad_csv <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(x = 1:2), bad_csv)
  expect_error(ipdcstorm::read_mdr_sst_csv(bad_csv), "Year column")

  expect_warning(
    anom <- ipdcstorm::compute_sst_anomaly(
      tibble::tibble(
        year = 2000:2008,
        sst_mdr_aso = seq(26, 27.6, length.out = 9)
      ),
      baseline_years = 2000:2008
    ),
    "Fewer than 10 years overlap"
  )
  expect_true(all(c("sst_clim", "sst_anomaly") %in% names(anom)))

  expect_equal(
    ipdcstorm::compute_p_hur_base(tibble::tibble(storm_class = c("TS", "HUR"), lambda = c(4, 1))),
    0.2
  )
  expect_equal(
    ipdcstorm::compute_p_hur_base(tibble::tibble(storm_class = c("TS", "HUR"), lambda = c(0, 0))),
    0.5
  )

  expect_true(all(ipdcstorm::sst_scenario_info("knmi23")$source == "knmi23"))
  expect_true(all(ipdcstorm::sst_scenario_info("ipcc_ar6")$source == "ipcc_ar6"))
  expect_true("knmi_Hn" %in% ipdcstorm::sst_scenario_info("all")$scenario)
  expect_equal(ipdcstorm::get_scenario_delta("knmi_Hn", target_year = 2020), 0)
  expect_equal(ipdcstorm::get_scenario_delta("knmi_Hn", target_year = 2050), 0.8)
  expect_equal(ipdcstorm::get_scenario_delta("knmi_Hn", target_year = 2150), 2.6)
  expect_error(ipdcstorm::get_scenario_delta("ssp245", target_year = 2050, baseline_years = 1:5), "length >= 10")
  expect_equal(ipdcstorm:::.validate_baseline_years(1991:2020), 1991:2020)
  expect_equal(ipdcstorm:::.resolve_climate_target(2050)$target_year, 2050)
  expect_null(ipdcstorm:::.resolve_climate_target()$target_year)
  expect_error(ipdcstorm:::.resolve_climate_target(NA_real_), "single finite numeric year")

  knmi <- ipdcstorm::knmi_scenario_info()
  expect_equal(nrow(knmi), 4)
  expect_setequal(knmi$variant, c("dry", "wet"))

  base_params <- ipdcstorm::default_cc_params()
  custom_params <- ipdcstorm::knmi_cc_params("knmi_Ln", base_params = base_params)
  expect_equal(custom_params$precip_scale, 0.10)
  expect_equal(custom_params$v_scale, base_params$v_scale)
  expect_error(ipdcstorm::knmi_cc_params("unknown"), "Unknown KNMI scenario")
})

test_that("climate fallback helpers cover insufficient-data and file validation paths", {
  annual_counts <- tibble::tibble(
    year = rep(2000:2003, each = 2),
    storm_class = rep(c("TS", "HUR"), times = 4),
    n_events = c(3L, 1L, 4L, 1L, 2L, 1L, 5L, 2L)
  )
  sst_df <- tibble::tibble(
    year = 2000:2003,
    sst_anomaly = c(-0.1, 0, 0.2, 0.3)
  )

  expect_warning(
    beta_info <- ipdcstorm::estimate_beta_sst(
      annual_counts = annual_counts,
      sst_df = sst_df,
      min_year = 2000L,
      verbose = FALSE
    ),
    "Only 4 years of overlapping data"
  )
  expect_equal(beta_info$method, "literature_fallback")
  expect_equal(beta_info$beta_sst, 0.6)
  expect_true(beta_info$guardrail$triggered)

  expect_error(
    ipdcstorm::read_mdr_sst_ersst("unused.nc"),
    "ERSST NetCDF not found"
  )
})

test_that("climate config normalizers handle direct and perturb edge cases", {
  expect_equal(
    ipdcstorm:::.normalize_climate_perturb(NULL, scenario = "stationary"),
    list(state = "disabled", params = NULL)
  )
  expect_equal(
    ipdcstorm:::.normalize_climate_perturb(list(), scenario = "stationary")$state,
    "default"
  )
  expect_error(
    ipdcstorm:::.normalize_climate_perturb("bad", scenario = "stationary"),
    "named list"
  )

  cfg_direct <- ipdcstorm:::.normalize_climate_cfg(list(
    scenario = NA_character_,
    sst_source = "builtin",
    baseline_years = 1991:2020,
    start_year = 1970L,
    delta_sst = 0.4,
    sensitivity_mode = "fixed",
    k_beta = 0,
    k_gamma = 0,
    target_year = NULL,
    perturb = NULL
  ))
  expect_identical(cfg_direct$input_mode, "direct_delta_sst")
  expect_identical(cfg_direct$perturb_state, "disabled")

  expect_error(
    ipdcstorm:::.normalize_climate_cfg(list(
      enabled = FALSE,
      scenario = "stationary",
      sst_source = "builtin",
      baseline_years = 1991:2020,
      start_year = 1970L,
      delta_sst = NULL,
      sensitivity_mode = "fixed",
      k_beta = 0,
      k_gamma = 0,
      target_year = NULL,
      perturb = NULL
    )),
    "Climate-off mode has been removed"
  )
})

test_that("downscale daily helpers and damage utilities handle common edge cases", {
  daily <- tibble::tibble(
    location = c("Saba", "Saba", "Statia"),
    sim_year = c(1L, 1L, 2L),
    scenario = c("baseline", "baseline", "future"),
    wind_kt = c(20, 55, NA_real_),
    wind_gust_kt = c(24, 66, 12),
    surge_m = c(0.2, 0.9, 0.1),
    event_class = c(NA, "TS", "HUR")
  )

  flags <- ipdcstorm:::disruption_flags(
    daily,
    thr_port = 34,
    thr_infra = 60,
    thr_surge = 0.5,
    use_gust = TRUE
  )
  expect_equal(flags$port_disrupt, c(FALSE, TRUE, FALSE))
  expect_equal(flags$infra_disrupt, c(FALSE, TRUE, FALSE))
  expect_equal(flags$surge_disrupt, c(FALSE, TRUE, FALSE))
  expect_equal(ipdcstorm:::is_tc_day(daily), c(FALSE, TRUE, TRUE))
  expect_equal(ipdcstorm:::is_hur_day(daily), c(NA, FALSE, TRUE))
  expect_equal(ipdcstorm:::exposure_hours(daily, threshold_kt = 34), c(0, 24, NA))

  expect_warning(
    peak <- ipdcstorm:::peak_wind_by_year(daily),
    "no non-missing arguments to max"
  )
  expect_equal(peak$peak_wind_kt[[1]], 55)
  expect_true(is.na(peak$peak_wind_kt[[2]]))

  pulse_tri <- ipdcstorm:::event_pulse(dur_days = 3, V_peak = 60, shape = "triangle")
  pulse_cos <- ipdcstorm:::event_pulse(dur_days = 0, V_peak = 60)
  expect_length(pulse_tri, 3L)
  expect_equal(max(pulse_tri), 60)
  expect_length(pulse_cos, 0L)

  dmg <- ipdcstorm:::add_damage_forcing(tibble::tibble(wind_kt = c(20, 40, 100)))
  expect_true(all(c("damage_intensity", "damage_rate") %in% names(dmg)))
  expect_equal(ipdcstorm:::damage_rate_from_wind(c(20, 80, 200)), c(0, 0.03, 0.10))
  expect_error(ipdcstorm:::damage_rate_from_wind(c(50), thr = 80, V_ref = 80), "V_ref must be > thr")

  spec <- ipdcstorm:::.validate_damage_spec(list(method = "powerlaw", p = 4))
  expect_equal(spec$method, "powerlaw")
  expect_equal(spec$p, 4)
  expect_error(ipdcstorm:::.validate_damage_spec(list()), "damage\\$method")
  expect_error(ipdcstorm:::.validate_damage_spec(list(method = "intensity", bad = 1)), "unsupported field")
})

test_that("downscale event-library wrappers and samplers produce deterministic schemas", {
  out <- list(
    trackpoints = list(
      Saba = tibble::tibble(
        SID = c("A", "A", "B", "B"),
        iso_time = as.POSIXct(c(
          "2000-08-01 00:00:00", "2000-08-01 06:00:00",
          "2000-09-10 00:00:00", "2000-09-10 06:00:00"
        ), tz = "UTC")
      )
    ),
    events = tibble::tibble(
      location = c("Saba", "Saba"),
      SID = c("A", "B"),
      peak_wind_kt = c(45, 80),
      storm_intensity_kt = c(45, 80),
      min_pressure_hpa = c(995, 970),
      rmw_mean_km = c(40, 25),
      start_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-09-10 00:00:00"), tz = "UTC")
    )
  )

  lib <- ipdcstorm:::build_event_library_from_out(out, location = "Saba", seed = 1L)
  sampled <- ipdcstorm:::sample_events_for_year_extended(lib, year = 2001, n_ts = 1, n_hur = 1, seed = 2L)

  expect_true(all(c("events", "sample_doy", "sample_event") %in% names(lib)))
  expect_equal(nrow(sampled), 2)
  expect_true(all(c("event_id", "event_class", "Pc_min_hPa", "RMW_mean_km") %in% names(sampled)))
  expect_error(ipdcstorm:::build_event_library_from_out(out, location = "Statia"), "no entry")
})

test_that("validation helper functions cover config, sampling, and empty-output paths", {
  expect_warning(
    cfg_warn <- ipdcstorm::make_validation_cfg(conf_level = 0.995),
    "outside the typical range"
  )
  expect_s3_class(cfg_warn, "validation_cfg")
  expect_error(ipdcstorm::make_validation_cfg(holdout_years = 0), "holdout_years must be >= 1")
  expect_error(ipdcstorm::make_validation_cfg(return_periods = 1), "must all be > 1")

  expect_equal(
    ipdcstorm:::.resolve_validation_n_sim(
      ipdcstorm::make_validation_cfg(n_sim = 120L),
      list()
    )$source,
    "validation_cfg$n_sim"
  )
  expect_equal(
    ipdcstorm:::.resolve_validation_n_sim(
      ipdcstorm::make_validation_cfg(n_sim = NULL),
      list(config = list(n_sim = 150L))
    )$n_sim,
    150L
  )
  expect_error(
    ipdcstorm:::.resolve_validation_n_sim(
      ipdcstorm::make_validation_cfg(n_sim = NULL),
      list(config = list(n_sim = 50L))
    ),
    ">= 100"
  )

  expect_equal(ipdcstorm:::.hindcast_sampler_mode("bounded"), "bounded")
  expect_error(ipdcstorm:::.hindcast_sampler_mode("bad"), "arg")

  kde_fallback <- ipdcstorm:::.fit_intensity_kde(c(40, 55), lower = 34, upper = 64)
  kde_empty <- ipdcstorm:::.fit_intensity_kde(numeric(0), lower = 34, upper = 64)
  expect_equal(kde_fallback$method, "fallback")
  expect_equal(kde_empty$pool_mean, 44)
  expect_length(ipdcstorm:::.sample_intensity_kde(kde_empty, 0), 0L)

  rl_small <- ipdcstorm:::compute_return_levels(c(1, 2), return_periods = c(2, 5))
  rl_gev_small <- ipdcstorm:::compute_return_levels_gev(c(0, 0, 10, 0), return_periods = c(2, 5))
  expect_true(all(is.na(rl_small)))
  expect_true(all(is.na(rl_gev_small$return_levels)))

  tmp_dir <- tempfile()
  tmp_csv <- tempfile(fileext = ".csv")
  expect_equal(ipdcstorm:::.validate_dir_create(tmp_dir), tmp_dir)
  ipdcstorm:::.validate_write_csv(tibble::tibble(x = 1:2), tmp_csv)
  expect_true(dir.exists(tmp_dir))
  expect_true(file.exists(tmp_csv))

  empty_val <- list(
    hindcast = list(comparison = tibble::tibble(), per_island = list()),
    rate_check = tibble::tibble(),
    wind_field = tibble::tibble(model_V_site_kt = numeric(0))
  )
  expect_null(ipdcstorm:::plot_hindcast_validation(empty_val))
  expect_null(ipdcstorm:::plot_rate_validation(empty_val))
  expect_null(ipdcstorm:::plot_wind_field_validation(empty_val))
  expect_null(ipdcstorm:::plot_bias_diagnostics(empty_val))
  expect_null(ipdcstorm:::plot_qq_validation(empty_val))
  expect_null(ipdcstorm:::plot_cdf_comparison(empty_val))
})

test_that("downscale hazard-impact wrapper executes intensity and powerlaw paths", {
  out_stub <- list(
    cfg = list(resampling_method = "stratified"),
    events = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("TS", "HUR"),
      SID = c("TS_1", "HUR_1")
    ),
    trackpoints = list(
      Saba = tibble::tibble(
        SID = c("TS_1", "HUR_1"),
        iso_time = as.POSIXct(
          c("2000-08-01 00:00:00", "2000-09-01 00:00:00"),
          tz = "UTC"
        )
      )
    ),
    sim = tibble::tibble(
      location = c("Saba", "Saba"),
      sim_year = c(1L, 2L),
      n_ts = c(1L, 0L),
      n_hur = c(0L, 1L)
    ),
    fit = structure(
      tibble::tibble(beta_sst = 0),
      delta_sst = 0.5,
      perturb = ipdcstorm::default_cc_params()
    )
  )

  sample_lib <- list(
    events = tibble::tibble(
      SID           = c("TS_1", "HUR_1"),
      storm_class   = c("TS", "HUR"),
      V_site_max_kt = c(45, 80),
      doy           = c(220L, 250L),
      dur_days      = c(3L, 4L),
      Pc_min_hPa    = c(995, 970),
      dP_max_hPa    = c(18, 40),
      RMW_mean_km   = c(40, 25)
    ),
    sample_doy = function(sev) if (identical(sev, "TS")) 220L else 250L,
    sample_event = function(sev) {
      tibble::tibble(
        SID = paste0(sev, "_1"),
        dur_days = if (identical(sev, "TS")) 3L else 4L,
        V_site_max_kt = if (identical(sev, "TS")) 45 else 80,
        Pc_min_hPa = if (identical(sev, "TS")) 995 else 970,
        dP_max_hPa = if (identical(sev, "TS")) 18 else 40,
        RMW_mean_km = if (identical(sev, "TS")) 40 else 25
      )
    }
  )

  local_mocked_bindings(
    build_event_library = function(...) sample_lib,
    .package = "ipdcstorm"
  )

  intensity <- ipdcstorm::generate_daily_hazard_impact_spatial(
    out = out_stub,
    location = "Saba",
    sim_years = 1:2,
    year0 = 2000,
    gust_factor = 1.2,
    damage = list(method = "intensity", V0 = 34, V1 = 120, p = 2, dmax = 0.03),
    pulse_shape = "triangle",
    scenario = "baseline",
    seed = 4
  )
  powerlaw <- ipdcstorm::generate_daily_hazard_impact_spatial(
    out = out_stub,
    location = "Saba",
    sim_years = 1:2,
    year0 = 2000,
    gust_factor = 1.1,
    damage = list(method = "powerlaw", thr = 34, V_ref = 80, d_ref = 0.03, p = 2, d_max = 0.08),
    pulse_shape = "cosine",
    scenario = "future",
    seed = 4
  )

  expect_true("Saba" %in% names(intensity))
  expect_true("Saba" %in% names(powerlaw))
  expect_true(all(c("wind_gust_kt", "damage_intensity", "damage_rate", "cum_damage") %in% names(intensity$Saba)))
  expect_true(all(intensity$Saba$scenario == "baseline"))
  expect_true(all(powerlaw$Saba$scenario == "future"))
  expect_equal(attr(intensity$Saba, "gust_factor"), 1.2)
  expect_true(any(!is.na(intensity$Saba$event_id)))
  expect_true(max(powerlaw$Saba$damage_rate, na.rm = TRUE) <= 0.08)
})

test_that("climate intensity and perturbation helpers cover operational branches", {
  annual_counts <- tibble::tibble(
    year = rep(2000:2011, each = 2),
    storm_class = rep(c("TS", "HUR"), times = 12),
    n_events = c(
      4L, 1L, 5L, 1L, 6L, 2L, 4L, 1L, 5L, 2L, 6L, 2L,
      7L, 3L, 6L, 2L, 7L, 3L, 8L, 3L, 7L, 2L, 8L, 4L
    )
  )
  sst_df <- tibble::tibble(
    year = 2000:2011,
    sst_anomaly = seq(-0.3, 0.8, length.out = 12)
  )

  gamma_info <- ipdcstorm::estimate_gamma_intensity(
    annual_counts = annual_counts,
    sst_df = sst_df,
    min_year = 2000L,
    gamma_prior = 0.05,
    verbose = FALSE
  )
  split <- ipdcstorm::compute_severity_split(
    lambda_table = tibble::tibble(storm_class = c("TS", "HUR"), lambda = c(0.8, 0.2)),
    sst_anomaly = c(0, 0.5, 1),
    gamma = gamma_info$gamma,
    p_hur_base = gamma_info$p_hur_base
  )
  perturbed <- ipdcstorm::perturb_event(
    events = tibble::tibble(
      V_peak = c(50, 80),
      dur_days = c(2L, 3L),
      RMW_mean_km = c(40, 25),
      storm_speed_kt = c(12, 15),
      precip_mm = c(10, 20)
    ),
    delta_sst = 0.5,
    cc_params = ipdcstorm::default_cc_params()
  )

  expect_equal(gamma_info$method, "binomial_glm")
  expect_true(all(c("p_hur", "lam_TS", "lam_HUR") %in% names(split)))
  expect_equal(split$lam_TS + split$lam_HUR, rep(1, 3))
  expect_true(all(c("precip_scaling", "delta_sst") %in% names(perturbed)))
  expect_true(all(perturbed$delta_sst == 0.5))
})

test_that("validation distribution helpers cover fitted branches beyond empty fallbacks", {
  fit <- ipdcstorm:::fit_gev_lmom(c(20, 25, 30, 35, 40, 45, 55))
  kde <- ipdcstorm:::.fit_intensity_kde(c(40, 45, 50, 55, 60, 62), lower = 34, upper = 64)
  set.seed(7)
  kde_draws <- ipdcstorm:::.sample_intensity_kde(kde, 25)
  gev_rl <- ipdcstorm:::compute_return_levels_gev(c(0, 20, 30, 40, 50, 65, 80, 90), return_periods = c(2, 5))
  set.seed(7)
  boot <- ipdcstorm:::bootstrap_return_level_ci(
    annual_max = c(0, 20, 30, 40, 50, 65, 80, 90),
    return_periods = c(2, 5),
    n_boot = 20
  )

  expect_true(isTRUE(fit$converged))
  expect_equal(kde$method, "kde")
  expect_length(kde_draws, 25L)
  expect_true(all(kde_draws >= 34 & kde_draws <= 64))
  expect_true(all(is.finite(gev_rl$return_levels)))
  expect_true(all(c("sim_median", "sim_ci_lo", "sim_ci_hi") %in% names(boot)))
})
