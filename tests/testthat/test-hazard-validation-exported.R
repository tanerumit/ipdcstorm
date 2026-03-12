validation_out_fixture <- function() {
  list(
    config = list(n_sim = 150L),
    cfg = list(n_sim = 150L),
    events = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("TS", "HUR"),
      storm_id = c("storm_ts", "2017242N16333"),
      year = c(2000L, 2017L),
      peak_wind_kt = c(50, 82),
      storm_intensity_kt = c(60, 125)
    ),
    rates = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("TS", "HUR"),
      lambda = c(1.0, 0.25),
      n_years = c(100L, 100L)
    ),
    lambda_scalers = tibble::tibble(
      location = c("Saba", "Saba"),
      storm_class = c("TS", "HUR"),
      lambda_model_raw = c(1.0, 0.25),
      lambda_ref = c(2.0, 0.6),
      expected_ratio = c(0.55, 0.30),
      lambda_target = c(1.1, 0.18),
      lambda_scale = c(1.1, 0.72),
      lambda_adj = c(1.1, 0.18),
      scale_status = c("ok", "ok"),
      scale_clamped = c(FALSE, FALSE)
    ),
    fit = tibble::tibble(beta_sst = 0, gamma_intensity = 0),
    trackpoints = list(
      Saba = tibble::tibble(
        SID = "2017242N16333",
        site_wind_kt = c(78, 80),
        dist_km = c(25, 20)
      )
    )
  )
}

validation_plot_fixture <- function() {
  list(
    hindcast = list(
      comparison = tibble::tibble(
        location = "Saba",
        return_period = 5,
        obs_full_rl = 82,
        sim_rl = 80,
        obs_ci_lo = 70,
        obs_ci_hi = 92,
        obs_in_ci = TRUE,
        bias_pct = -2
      ),
      per_island = list(
        Saba = list(
          obs_annual_max = tibble::tibble(V_max_kt = c(0, 45, 60, 75, 90, 100)),
          sim_annual_max = c(0, 40, 62, 72, 88, 96),
          gev_fit = list(
            gev_fit = list(mu = 70, sigma = 10, xi = 0.1),
            p_zero = 1 / 6
          ),
          obs_gev = list(
            gev_fit = list(xi = 0.1),
            n_nonzero = 5
          ),
          train_years = 2000:2004,
          test_years = 2005:2006
        )
      )
    ),
    rate_check = tibble::tibble(
      location = "Saba",
      storm_class = "TS",
      lambda_ref = 0.55,
      lambda_adj = 0.60,
      flag = "OK"
    ),
    wind_field = tibble::tibble(
      location = "Saba",
      storm_name = "IRMA",
      obs_1min_equiv_kt = 80,
      model_V_site_kt = 78
    )
  )
}

hindcast_retention_fixture <- function() {
  list(
    events = tibble::tibble(
      location = rep("Saba", 4),
      storm_id = c("s1", "s2", "s3", "s4"),
      year = c(2001L, 2002L, 2003L, 2004L),
      storm_class = c("TD", "TS", "HUR", "TS"),
      peak_wind_kt = c(33, 36, 70, 45)
    ),
    trackpoints = tibble::tibble(
      SID = c("s1", "s2", "s3", "s4"),
      iso_time = as.POSIXct(
        c("2001-08-01 00:00:00", "2002-08-01 00:00:00", "2003-08-01 00:00:00", "2004-08-01 00:00:00"),
        tz = "UTC"
      ),
      site_wind_kt = c(33, 36, 70, 45),
      R34_source = c("none", "partial", "climo", "observed"),
      R34_eff_km = c(NA_real_, 120, 150, 180),
      RMW_used_km = c(20, 25, 30, 35)
    ),
    metadata = list(
      model_seed = 11L,
      validation_seed = 22L,
      data_id = "ibtracs-fixture|rows=4",
      parameter_id = "params-fixture",
      lambda_scaler_id = "lambda-fixture"
    )
  )
}

test_that("validation configs and return-level helpers expose stable public results", {
  cfg <- make_validation_cfg(n_sim = 120L, save_plots = FALSE, save_tables = FALSE)
  rl_emp <- compute_return_levels(c(0, 30, 40, 50, 60, 70), return_periods = c(2, 5))
  gev_fit <- fit_gev_lmom(c(20, 25, 30, 35, 40, 45))
  rl_gev <- compute_return_levels_gev(c(0, 0, 35, 40, 60, 80, 90), return_periods = c(2, 5))

  expect_s3_class(cfg, "validation_cfg")
  expect_output(print(cfg), "Validation configuration")
  expect_true(cfg$advanced$hindcast_use_raw_rates)
  expect_equal(names(rl_emp), c("RL_2yr", "RL_5yr"))
  expect_true(isTRUE(gev_fit$converged))
  expect_true(all(c("return_levels", "p_zero", "n_total", "n_nonzero") %in% names(rl_gev)))
})

test_that("hindcast bias decomposition exposes frequency, intensity, and metadata fields", {
  decomp <- ipdcstorm:::.compute_hindcast_bias_decomposition(
    obs_annual_max = c(0, 40, 60, 0),
    sim_annual_max = c(0, 50, 70, 30),
    location = "Saba",
    metadata = list(
      model_seed = 101L,
      validation_seed = 202L,
      data_id = "ibtracs-a",
      parameter_id = "params-a",
      lambda_scaler_id = "lambda-a"
    )
  )

  expect_equal(decomp$dominant_component, "frequency")
  expect_equal(decomp$freq_contrib_kt, 12.5, tolerance = 1e-12)
  expect_equal(decomp$intensity_contrib_kt, 0, tolerance = 1e-12)
  expect_equal(decomp$interaction_contrib_kt, 0, tolerance = 1e-12)
  expect_equal(decomp$model_seed, 101L)
  expect_equal(decomp$validation_seed, 202L)
  expect_equal(decomp$data_id, "ibtracs-a")
  expect_equal(decomp$parameter_id, "params-a")
  expect_equal(decomp$lambda_scaler_id, "lambda-a")
})

test_that("hindcast retention diagnostics stratify by R34 source and summarize threshold retention", {
  fix <- hindcast_retention_fixture()

  diag <- ipdcstorm:::.summarize_hindcast_retention(
    events_island = fix$events,
    trackpoints_island = fix$trackpoints,
    train_years = 2001:2002,
    test_years = 2003:2004,
    location = "Saba",
    metadata = fix$metadata
  )

  expect_equal(sort(unique(diag$r34_source$peak_r34_source)), c("climo", "none", "observed", "partial"))

  train_summary <- dplyr::filter(diag$summary, period == "train")
  test_summary <- dplyr::filter(diag$summary, period == "test")
  expect_equal(train_summary$zero_event_years, 1L)
  expect_equal(train_summary$ts_years, 1L)
  expect_equal(train_summary$hur_years, 0L)
  expect_equal(train_summary$near_ts_threshold_years, 1L)
  expect_match(train_summary$near_ts_definition, "annual_max_kt")
  expect_equal(test_summary$zero_event_years, 0L)
  expect_equal(test_summary$ts_years, 1L)
  expect_equal(test_summary$hur_years, 1L)

  partial_row <- dplyr::filter(diag$r34_source, period == "train", peak_r34_source == "partial")
  none_row <- dplyr::filter(diag$r34_source, period == "train", peak_r34_source == "none")
  climo_row <- dplyr::filter(diag$r34_source, period == "test", peak_r34_source == "climo")
  observed_row <- dplyr::filter(diag$r34_source, period == "test", peak_r34_source == "observed")
  expect_equal(partial_row$n_tsplus_events, 1L)
  expect_equal(partial_row$n_annual_max_years, 1L)
  expect_equal(none_row$n_tsplus_events, 0L)
  expect_equal(climo_row$n_hur_events, 1L)
  expect_equal(observed_row$n_tsplus_events, 1L)

  expect_true(all(diag$summary$model_seed == 11L))
  expect_true(all(diag$r34_source$validation_seed == 22L))
  expect_true(all(diag$event_provenance$data_id == "ibtracs-fixture|rows=4"))
  expect_true(all(diag$yearly$parameter_id == "params-fixture"))
  expect_true(all(diag$yearly$lambda_scaler_id == "lambda-fixture"))
})

test_that("wind-mode retention comparison stays deterministic by source and period", {
  fix <- hindcast_retention_fixture()
  legacy <- ipdcstorm:::.summarize_hindcast_retention(
    events_island = fix$events,
    trackpoints_island = fix$trackpoints,
    train_years = 2001:2002,
    test_years = 2003:2004,
    location = "Saba",
    metadata = fix$metadata
  )$event_provenance

  diagnostic <- legacy
  diagnostic$retained_tsplus[diagnostic$storm_id == "s1"] <- TRUE
  diagnostic$retained_hur[diagnostic$storm_id == "s4"] <- TRUE

  cmp <- ipdcstorm:::.compare_hindcast_retention_by_wind(
    legacy_events = legacy,
    diagnostic_events = diagnostic
  )

  none_row <- dplyr::filter(cmp, location == "Saba", period == "train", peak_r34_source == "none")
  observed_row <- dplyr::filter(cmp, location == "Saba", period == "test", peak_r34_source == "observed")
  expect_equal(none_row$delta_tsplus_events, 1L)
  expect_equal(none_row$tsplus_flip_events, 1L)
  expect_equal(observed_row$delta_hur_events, 1L)
  expect_equal(observed_row$hur_flip_events, 1L)
  expect_true(all(cmp$model_seed[cmp$location == "Saba"] == 11L))
})

test_that("empirical hindcast intensity sampling can be bounded for diagnostics", {
  fit <- ipdcstorm:::.fit_intensity_kde(c(40, 55, 63), lower = 34, upper = 64)

  old_opt <- options(ipdcstorm.hindcast_sampler_mode = "bounded")
  on.exit(options(old_opt), add = TRUE)

  set.seed(42)
  draws_a <- ipdcstorm:::.sample_intensity_kde(fit, 50)
  set.seed(42)
  draws_b <- ipdcstorm:::.sample_intensity_kde(fit, 50)

  expect_equal(draws_a, draws_b)
  expect_true(all(draws_a >= min(fit$pool)))
  expect_lte(max(draws_a), max(fit$pool))
})

test_that("hindcast attribution grid records mode metadata from workspace reruns", {
  out_stub <- validation_out_fixture()
  hc_stub <- list(
    comparison = tibble::tibble(
      location = c("Saba", "Saba", "Statia", "Statia", "St_Martin", "St_Martin"),
      return_period = c(5, 10, 5, 10, 5, 10),
      obs_full_rl = c(80, 90, 75, 85, 95, 105),
      obs_test_rl = c(78, 88, 72, 82, 92, 102),
      sim_rl = c(82, 91, 79, 87, 98, 110),
      sim_median = c(82, 91, 79, 87, 98, 110),
      sim_ci_lo = NA_real_,
      sim_ci_hi = NA_real_,
      obs_ci_lo = NA_real_,
      obs_ci_hi = NA_real_,
      obs_ci_method = "none",
      obs_in_ci = NA,
      model_in_obs_ci = NA,
      obs_in_model_ci = NA,
      bias_pct = c(2, 1, 5, 2, 3, 5)
    ),
    per_island = list(
      Saba = list(skipped = FALSE, gev_fit = list(gev_fit = list(xi = 0.12)), obs_gev = list(gev_fit = list(xi = 0.08))),
      Statia = list(skipped = FALSE, gev_fit = list(gev_fit = list(xi = 0.14)), obs_gev = list(gev_fit = list(xi = 0.09))),
      St_Martin = list(skipped = FALSE, gev_fit = list(gev_fit = list(xi = 0.10)), obs_gev = list(gev_fit = list(xi = 0.07)))
    )
  )

  local_mocked_bindings(
    run_hazard_model = function(cfg, targets, storm_classes, climate = NULL, seed = NULL, verbose = FALSE) {
      out_stub$run_metadata <- list(
        seed = seed,
        ibtracs_data_id = "ibtracs_fixture|rows=2",
        parameter_id = "params-fixture"
      )
      out_stub$lambda_scaler_id <- "lambda-fixture"
      out_stub
    },
    .validate_hindcast_all = function(out, ...) hc_stub,
    .package = "ipdcstorm"
  )

  grid <- ipdcstorm:::.run_hindcast_attribution_grid(
    cfg = make_hazard_cfg(simulation_years = 120L),
    targets = tibble::tibble(
      location = c("Saba", "Statia", "St_Martin"),
      lat = c(17.6350, 17.4890, 18.0708),
      lon = c(-63.2300, -62.9740, -63.0501)
    ),
    validation_cfg = make_validation_cfg(
      holdout_years = 10L,
      n_sim = 120L,
      return_periods = c(5, 10),
      seed = 99L,
      save_plots = FALSE,
      save_tables = FALSE
    ),
    model_seed = 77L
  )

  expect_true(all(c("wind_rate_grid", "sampler_grid", "summary", "metadata") %in% names(grid)))
  expect_true(all(c("baseline_case_id", "baseline_diagnostics", "case_diagnostics", "wind_retention_comparison") %in% names(grid)))
  expect_true(all(c("case_id", "data_id", "parameter_id", "model_seed", "validation_seed") %in% names(grid$metadata)))
  expect_true(all(c("rl_bias_rp5", "rl_bias_rp10", "sim_xi", "obs_xi") %in% names(grid$summary)))
  expect_equal(sort(unique(grid$wind_rate_grid$annual_rate_mode)), c("adjusted", "raw"))
  expect_equal(sort(unique(grid$sampler_grid$sampler_mode)), c("bounded", "legacy"))
  expect_equal(grid$baseline_case_id, "wind=legacy|rate=raw|sampler=legacy")
})

test_that("bootstrap and reference data helpers return expected schema", {
  set.seed(1)
  ci <- bootstrap_return_level_ci(
    annual_max = c(0, 20, 30, 40, 50, 60, 75),
    return_periods = c(2, 5),
    n_boot = 20
  )
  ref_rates <- get_reference_rates()
  obs <- get_wind_observations()

  expect_equal(ci$return_period, c(2, 5))
  expect_true(all(c("sim_ci_lo", "sim_ci_hi", "sim_lo_50", "sim_hi_50") %in% names(ci)))
  expect_true(all(c("region", "storm_class", "lambda_ref") %in% names(ref_rates)))
  expect_true(all(c("storm_sid", "target_island", "obs_quality") %in% names(obs)))
})

test_that("rate and wind-field validation work with lightweight fixtures", {
  out <- validation_out_fixture()
  obs_tbl <- tibble::tibble(
    storm_sid = "2017242N16333",
    storm_name = "IRMA",
    year = 2017L,
    target_island = "Saba",
    station = "test_station",
    obs_wind_kt = 80,
    obs_type = "1min_sust",
    obs_quality = "A",
    obs_source = "fixture",
    notes = "fixture"
  )

  expect_message(
    rates <- validate_rates(out),
    "\\[Rate Check\\] Summary:"
  )
  expect_message(
    wf <- validate_wind_field(out, obs_table = obs_tbl),
    "\\[Wind Field Check\\]"
  )

  expect_true(all(c("raw_ratio", "adj_ratio", "flag") %in% names(rates)))
  expect_equal(wf$model_V_site_kt, 82)
  expect_true(is.logical(wf$bias_ok))
})

test_that("plot exporters save validation artifacts", {
  val <- validation_plot_fixture()
  cfg <- make_validation_cfg(
    n_sim = 120L,
    out_dir = file.path(tempdir(), "validation-plots"),
    save_plots = FALSE,
    save_tables = FALSE
  )

  paths_hindcast <- plot_hindcast_validation(val, cfg = cfg)
  path_rate <- plot_rate_validation(val, cfg = cfg)
  paths_wind <- plot_wind_field_validation(val, cfg = cfg)
  paths_bias <- plot_bias_diagnostics(val, cfg = cfg)
  path_qq <- plot_qq_validation(val, cfg = cfg)
  path_cdf <- plot_cdf_comparison(val, cfg = cfg)

  expect_true(file.exists(unname(paths_hindcast[[1]])))
  expect_true(file.exists(path_rate))
  expect_true(file.exists(unname(paths_wind[[1]])))
  expect_true(file.exists(unname(paths_bias[[1]])))
  expect_true(file.exists(path_qq))
  expect_true(file.exists(path_cdf))
})

test_that("validation wrappers return structured outputs and normalize forwarded storm classes", {
  out_stub <- validation_out_fixture()
  cfg <- make_validation_cfg(n_sim = NULL, save_plots = FALSE, save_tables = FALSE)

  expect_message(
    val <- run_validation_suite(out_stub, cfg = cfg),
    "HAZARD MODEL VALIDATION SUITE"
  )
  expect_true(all(c("hindcast", "rate_check", "wind_field", "summary", "artifacts") %in% names(val)))

  local_mocked_bindings(
    run_hazard_model = function(cfg, targets, storm_classes, climate = NULL) {
      list(cfg = cfg, targets = targets, storm_classes = storm_classes, climate = climate)
    },
    run_validation_suite = function(out, cfg) {
      list(summary = tibble::tibble(location = "Saba"), artifacts = list(plots = list(), tables = list()), out = out, cfg = cfg)
    },
    .package = "ipdcstorm"
  )

  result <- validate_hazard_model(
    cfg = make_hazard_cfg(simulation_years = 120L),
    targets = tibble::tibble(location = "Saba", lat = 17.63, lon = -63.23, search_radius_km = 200),
    validation_cfg = make_validation_cfg(n_sim = 120L, save_plots = FALSE, save_tables = FALSE),
    storm_classes = "TS"
  )

  expect_equal(result$out$storm_classes, "TS")
  expect_true("summary" %in% names(result$val))
})
