# =============================================================================
# Tests for spatial-coherence downscale extensions:
#   - .ev_dur_days / .ev_V_peak / .ev_safe_num / .ev_doy  (event-row helpers)
#   - make_background_wind_cfg() validation + structure
#   - .generate_background_wind_year() output properties
#   - background_wind integration in generate_daily_hazard_impact_spatial()
#   - pinned_sids validation + functional guarantee
#   - perturb = TRUE / FALSE shortcuts in make_climate_cfg()
#   - KNMI warning when perturb is disabled
# =============================================================================

# ---------------------------------------------------------------------------
# Shared minimal out fixture compatible with generate_daily_hazard_impact_spatial
# ---------------------------------------------------------------------------

.new_feat_out <- function() {
  track <- tibble::tibble(
    SID = c("storm_ts", "storm_ts", "storm_hur", "storm_hur"),
    iso_time = as.POSIXct(
      c("2000-08-01 00:00:00", "2000-08-01 06:00:00",
        "2000-09-10 00:00:00", "2000-09-10 06:00:00"),
      tz = "UTC"
    )
  )
  events <- tibble::tibble(
    location         = c("Saba", "Saba"),
    storm_id         = c("storm_ts", "storm_hur"),
    peak_wind_kt     = c(45, 85),
    storm_intensity_kt = c(55, 110),
    min_pressure_hpa = c(995, 955),
    pressure_deficit_hpa = c(15, 55),
    rmw_mean_km      = c(40, 24),
    start_time = as.POSIXct(c("2000-08-01 00:00:00", "2000-09-10 00:00:00"), tz = "UTC"),
    end_time   = as.POSIXct(c("2000-08-02 00:00:00", "2000-09-12 00:00:00"), tz = "UTC"),
    n_points   = c(2L, 2L)
  )
  list(
    cfg         = list(resampling_method = "stratified"),
    trackpoints = list(Saba = track),
    events      = events,
    sim         = tibble::tibble(
      location = "Saba",
      sim_year = 1L,
      n_ts     = 1L,
      n_hur    = 1L
    )
  )
}

# Helper: run spatial gen, suppressing common spurious warnings from tiny fixtures
.run_spatial <- function(...) suppressWarnings(suppressMessages(generate_daily_hazard_impact_spatial(...)))

# =============================================================================
# 1) Module-level helpers: .ev_dur_days
# =============================================================================

test_that(".ev_dur_days returns dur_days when present and valid", {
  row <- tibble::tibble(dur_days = 5L)
  expect_identical(ipdcstorm:::.ev_dur_days(row), 5L)
})

test_that(".ev_dur_days falls back to start_time/end_time difference", {
  row <- tibble::tibble(
    start_time = as.POSIXct("2000-09-01", tz = "UTC"),
    end_time   = as.POSIXct("2000-09-04", tz = "UTC")  # 3 days apart → floor(3)+1 = 4
  )
  expect_identical(ipdcstorm:::.ev_dur_days(row), 4L)
})

test_that(".ev_dur_days falls back to n_points / 4 ceiling", {
  row <- tibble::tibble(n_points = 10L)
  expect_identical(ipdcstorm:::.ev_dur_days(row), 3L)  # ceiling(10/4) = 3
})

test_that(".ev_dur_days returns 1L when all columns missing or invalid", {
  row <- tibble::tibble(x = 1)
  expect_identical(ipdcstorm:::.ev_dur_days(row), 1L)

  row_bad <- tibble::tibble(dur_days = -3L)
  expect_identical(ipdcstorm:::.ev_dur_days(row_bad), 1L)
})

# =============================================================================
# 2) Module-level helpers: .ev_V_peak
# =============================================================================

test_that(".ev_V_peak reads V_site_max_kt first", {
  row <- tibble::tibble(V_site_max_kt = 75, wind_max_kt = 60, peak_wind_kt = 50)
  expect_equal(ipdcstorm:::.ev_V_peak(row, "HUR"), 75)
})

test_that(".ev_V_peak falls through to wind_max_kt then peak_wind_kt", {
  row1 <- tibble::tibble(wind_max_kt = 65, peak_wind_kt = 50)
  expect_equal(ipdcstorm:::.ev_V_peak(row1, "TS"), 65)

  row2 <- tibble::tibble(peak_wind_kt = 55)
  expect_equal(ipdcstorm:::.ev_V_peak(row2, "HUR"), 55)
})

test_that(".ev_V_peak uses severity-aware fallback when all columns are missing/invalid", {
  row <- tibble::tibble(x = 1)
  expect_equal(ipdcstorm:::.ev_V_peak(row, "HUR"), 80)
  expect_equal(ipdcstorm:::.ev_V_peak(row, "TS"),  40)
  expect_equal(ipdcstorm:::.ev_V_peak(row, "TD"),  25)
})

test_that(".ev_V_peak uses fallback when stored value is zero or non-finite", {
  row <- tibble::tibble(V_site_max_kt = 0)
  expect_equal(ipdcstorm:::.ev_V_peak(row, "HUR"), 80)

  row_na <- tibble::tibble(V_site_max_kt = NA_real_)
  expect_equal(ipdcstorm:::.ev_V_peak(row_na, "TS"), 40)
})

# =============================================================================
# 3) Module-level helpers: .ev_safe_num and .ev_doy
# =============================================================================

test_that(".ev_safe_num returns finite numeric or NA", {
  row <- tibble::tibble(pressure = 1013.2, bad_col = "text")
  expect_equal(ipdcstorm:::.ev_safe_num(row, "pressure"), 1013.2)
  expect_true(is.na(ipdcstorm:::.ev_safe_num(row, "missing_col")))
  expect_true(is.na(ipdcstorm:::.ev_safe_num(row, "bad_col")))
})

test_that(".ev_doy reads doy column when valid", {
  row <- tibble::tibble(doy = 200L)
  expect_identical(ipdcstorm:::.ev_doy(row), 200L)
})

test_that(".ev_doy falls back to start_time then to 220L", {
  row_time <- tibble::tibble(
    start_time = as.POSIXct("2000-08-01", tz = "UTC")  # DOY 214 in non-leap year
  )
  expect_identical(ipdcstorm:::.ev_doy(row_time), 214L)

  row_none <- tibble::tibble(x = 1)
  expect_identical(ipdcstorm:::.ev_doy(row_none), 220L)

  row_oob <- tibble::tibble(doy = 400L)
  expect_identical(ipdcstorm:::.ev_doy(row_oob), 220L)
})

# =============================================================================
# 4) make_background_wind_cfg() — validation
# =============================================================================

.weibull_single <- function(loc = "Saba") {
  stats::setNames(
    list(data.frame(month = 1:12, shape = rep(1.5, 12), scale = rep(8, 12))),
    loc
  )
}

.weibull_two <- function() {
  list(
    Saba  = data.frame(month = 1:12, shape = rep(1.5, 12), scale = rep(8, 12)),
    Statia = data.frame(month = 1:12, shape = rep(1.4, 12), scale = rep(7, 12))
  )
}

test_that("make_background_wind_cfg returns a background_wind_cfg object", {
  cfg <- make_background_wind_cfg(.weibull_single())
  expect_s3_class(cfg, "background_wind_cfg")
  expect_true(all(c("weibull_params", "cor_matrix", "chol_U", "ar1", "locations") %in% names(cfg)))
})

test_that("make_background_wind_cfg defaults to identity cor_matrix", {
  cfg <- make_background_wind_cfg(.weibull_single())
  expect_equal(cfg$cor_matrix, diag(1), ignore_attr = TRUE)
  expect_equal(cfg$ar1, 0)
})

test_that("make_background_wind_cfg accepts a valid 2x2 cor_matrix", {
  cor_mat <- matrix(c(1, 0.6, 0.6, 1), nrow = 2,
                    dimnames = list(c("Saba", "Statia"), c("Saba", "Statia")))
  cfg <- make_background_wind_cfg(.weibull_two(), cor_matrix = cor_mat)
  expect_equal(cfg$cor_matrix, cor_mat)
  expect_equal(dim(cfg$chol_U), c(2L, 2L))
})

test_that("make_background_wind_cfg rejects unnamed weibull_params", {
  bad <- list(data.frame(month = 1:12, shape = 1.5, scale = 8))
  expect_error(make_background_wind_cfg(bad), "named list")
})

test_that("make_background_wind_cfg rejects missing required columns", {
  bad <- list(Saba = data.frame(month = 1:12, shape = 1.5))  # no scale
  expect_error(make_background_wind_cfg(bad), "month, shape, scale")
})

test_that("make_background_wind_cfg rejects invalid month values", {
  bad <- list(Saba = data.frame(month = c(0, 1:11), shape = 1.5, scale = 8))
  expect_error(make_background_wind_cfg(bad), "month")
})

test_that("make_background_wind_cfg rejects non-positive shape or scale", {
  bad_shape <- list(Saba = data.frame(month = 1:12, shape = c(-1, rep(1.5, 11)), scale = 8))
  expect_error(make_background_wind_cfg(bad_shape), "shape")

  bad_scale <- list(Saba = data.frame(month = 1:12, shape = 1.5, scale = c(0, rep(8, 11))))
  expect_error(make_background_wind_cfg(bad_scale), "scale")
})

test_that("make_background_wind_cfg rejects non-symmetric cor_matrix", {
  asym <- matrix(c(1, 0.7, 0.3, 1), nrow = 2,
                 dimnames = list(c("Saba", "Statia"), c("Saba", "Statia")))
  expect_error(make_background_wind_cfg(.weibull_two(), cor_matrix = asym), "symmetric")
})

test_that("make_background_wind_cfg rejects non-unit diagonal", {
  bad_diag <- matrix(c(2, 0.6, 0.6, 2), nrow = 2,
                     dimnames = list(c("Saba", "Statia"), c("Saba", "Statia")))
  expect_error(make_background_wind_cfg(.weibull_two(), cor_matrix = bad_diag), "unit diagonal")
})

test_that("make_background_wind_cfg rejects ar1 out of [0, 1)", {
  expect_error(make_background_wind_cfg(.weibull_single(), ar1 = -0.1), "ar1")
  expect_error(make_background_wind_cfg(.weibull_single(), ar1 = 1.0),  "ar1")
  expect_error(make_background_wind_cfg(.weibull_single(), ar1 = 1.5),  "ar1")
})

# =============================================================================
# 5) .generate_background_wind_year() — output properties
# =============================================================================

test_that(".generate_background_wind_year returns 365 values for non-leap year", {
  cfg <- make_background_wind_cfg(.weibull_single())
  set.seed(1)
  bg <- ipdcstorm:::.generate_background_wind_year(2001, "Saba", cfg)
  expect_length(bg$Saba, 365L)
})

test_that(".generate_background_wind_year returns 366 values for leap year", {
  cfg <- make_background_wind_cfg(.weibull_single())
  set.seed(1)
  bg <- ipdcstorm:::.generate_background_wind_year(2000, "Saba", cfg)
  expect_length(bg$Saba, 366L)
})

test_that(".generate_background_wind_year output is non-negative", {
  cfg <- make_background_wind_cfg(.weibull_single())
  set.seed(1)
  bg <- ipdcstorm:::.generate_background_wind_year(2001, "Saba", cfg)
  expect_true(all(bg$Saba >= 0))
})

test_that(".generate_background_wind_year with correlated matrix preserves marginals", {
  cor_mat <- matrix(c(1, 0.8, 0.8, 1), nrow = 2,
                    dimnames = list(c("Saba", "Statia"), c("Saba", "Statia")))
  cfg <- make_background_wind_cfg(.weibull_two(), cor_matrix = cor_mat)
  set.seed(7)
  bg <- ipdcstorm:::.generate_background_wind_year(2001, c("Saba", "Statia"), cfg)
  # Both locations should have 365 non-negative values
  expect_length(bg$Saba,   365L)
  expect_length(bg$Statia, 365L)
  expect_true(all(bg$Saba   >= 0))
  expect_true(all(bg$Statia >= 0))
  # With high correlation, the two series should be positively correlated
  expect_gt(cor(bg$Saba, bg$Statia), 0.5)
})

test_that(".generate_background_wind_year AR(1) persistence increases autocorrelation", {
  cfg_ar0 <- make_background_wind_cfg(.weibull_single(), ar1 = 0)
  cfg_ar8 <- make_background_wind_cfg(.weibull_single(), ar1 = 0.8)
  set.seed(42)
  bg0 <- ipdcstorm:::.generate_background_wind_year(2001, "Saba", cfg_ar0)$Saba
  set.seed(42)
  bg8 <- ipdcstorm:::.generate_background_wind_year(2001, "Saba", cfg_ar8)$Saba
  # AR(1) series should have higher lag-1 autocorrelation
  acf0 <- acf(bg0, lag.max = 1, plot = FALSE)$acf[2]
  acf8 <- acf(bg8, lag.max = 1, plot = FALSE)$acf[2]
  expect_gt(acf8, acf0)
})

# =============================================================================
# 6) background_wind integration in generate_daily_hazard_impact_spatial()
# =============================================================================

test_that("background_wind parameter is validated: must be background_wind_cfg", {
  out <- .new_feat_out()
  expect_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 background_wind = list(not = "a cfg")),
    "background_wind_cfg"
  )
})

test_that("background_wind with missing location raises informative error", {
  out <- .new_feat_out()
  cfg <- make_background_wind_cfg(.weibull_single("Other"))  # wrong location name
  expect_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 background_wind = cfg),
    "Saba"  # error should name the missing location
  )
})

test_that("background_wind integration: wind_kt is non-negative with background enabled", {
  out <- .new_feat_out()
  bg_cfg <- make_background_wind_cfg(.weibull_single("Saba"), ar1 = 0.5)
  res <- .run_spatial(out, location = "Saba", sim_years = 1L,
                      background_wind = bg_cfg, seed = 10L)
  expect_true(all(res$Saba$wind_kt >= 0, na.rm = TRUE))
  expect_equal(nrow(res$Saba), 366L)  # year0=2000 is a leap year
})

test_that("background_wind raises wind floor: min wind_kt > 0 on non-storm days", {
  out <- .new_feat_out()
  bg_cfg <- make_background_wind_cfg(.weibull_single("Saba"))

  res_no_bg <- .run_spatial(out, location = "Saba", sim_years = 1L, seed = 5L)
  res_bg    <- .run_spatial(out, location = "Saba", sim_years = 1L,
                            background_wind = bg_cfg, seed = 5L)

  # With background wind, minimum daily wind must be strictly positive
  expect_gt(min(res_bg$Saba$wind_kt, na.rm = TRUE), 0)
  # Without background, many days will be exactly zero (no storm)
  expect_true(any(res_no_bg$Saba$wind_kt == 0, na.rm = TRUE))
})

# =============================================================================
# 7) pinned_sids — validation
# =============================================================================

test_that("pinned_sids validation rejects non-list input", {
  out <- .new_feat_out()
  expect_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 pinned_sids = "storm_hur"),
    "named list"
  )
})

test_that("pinned_sids validation rejects multi-element character entries", {
  out <- .new_feat_out()
  expect_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 pinned_sids = list("1" = c("storm_ts", "storm_hur"))),
    "single character SID"
  )
})

test_that("pinned_sids = NULL is accepted (no pinning)", {
  out <- .new_feat_out()
  expect_no_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 pinned_sids = NULL, seed = 1L)
  )
})

test_that("pinned_sids with NA entry is accepted (year skips pinning)", {
  out <- .new_feat_out()
  expect_no_error(
    .run_spatial(out, location = "Saba", sim_years = 1L,
                 pinned_sids = list("1" = NA_character_), seed = 1L)
  )
})

# =============================================================================
# 8) pinned_sids — functional: pinned SID must appear in output
# =============================================================================

test_that("pinned_sids guarantees the focal SID appears in the output event_id", {
  # The fixture library contains "storm_hur" as the HUR-class event.
  # With n_hur = 1, pinning "storm_hur" should always produce an event_id
  # that starts with "storm_hur" (suffix: _y<year>_<counter>).
  out <- .new_feat_out()

  res <- .run_spatial(
    out, location = "Saba", sim_years = 1L,
    pinned_sids = list("1" = "storm_hur"),
    seed = 99L
  )

  event_ids  <- stats::na.omit(res$Saba$event_id)
  hur_events <- event_ids[startsWith(event_ids, "storm_hur")]
  expect_gt(length(hur_events), 0L)
})

test_that("pinned_sids result is reproducible across equivalent calls", {
  out <- .new_feat_out()
  pins <- list("1" = "storm_hur")

  res1 <- .run_spatial(out, location = "Saba", sim_years = 1L,
                       pinned_sids = pins, seed = 7L)
  res2 <- .run_spatial(out, location = "Saba", sim_years = 1L,
                       pinned_sids = pins, seed = 7L)

  expect_identical(res1$Saba$event_id, res2$Saba$event_id)
})

# =============================================================================
# 9) perturb = TRUE / FALSE shortcuts in make_climate_cfg()
# =============================================================================

test_that("perturb = TRUE produces the same config as perturb = list()", {
  cfg_true <- make_climate_cfg(scenario = "stationary", perturb = TRUE)
  cfg_list <- make_climate_cfg(scenario = "stationary", perturb = list())
  # Both should have the same perturb_state in the resolved config
  expect_identical(cfg_true$perturb_state, cfg_list$perturb_state)
})

test_that("perturb = FALSE produces the same config as perturb = NULL", {
  cfg_false <- make_climate_cfg(scenario = "stationary", perturb = FALSE)
  cfg_null  <- make_climate_cfg(scenario = "stationary", perturb = NULL)
  expect_identical(cfg_false$perturb_state, cfg_null$perturb_state)
  expect_identical(cfg_false$perturb_state, "disabled")
})

test_that("perturb = TRUE enables perturbation (state != disabled)", {
  cfg <- make_climate_cfg(scenario = "stationary", perturb = TRUE)
  expect_false(cfg$perturb_state == "disabled")
})

# =============================================================================
# 10) KNMI warning when perturb is disabled
# =============================================================================

test_that("make_climate_cfg emits a warning for KNMI scenario with perturb = NULL", {
  expect_warning(
    make_climate_cfg(scenario = "knmi_Hd", target_year = 2100, perturb = NULL),
    regexp = "perturbation is disabled",
    ignore.case = TRUE
  )
})

test_that("make_climate_cfg emits a warning for KNMI scenario with perturb = FALSE", {
  expect_warning(
    make_climate_cfg(scenario = "knmi_Ld", target_year = 2050, perturb = FALSE),
    regexp = "perturbation is disabled",
    ignore.case = TRUE
  )
})

test_that("make_climate_cfg does NOT warn for KNMI scenario with perturb = TRUE", {
  expect_no_warning(
    make_climate_cfg(scenario = "knmi_Hd", target_year = 2100, perturb = TRUE)
  )
})

test_that("make_climate_cfg does NOT warn for stationary scenario with perturb = NULL", {
  expect_no_warning(
    make_climate_cfg(scenario = "stationary", perturb = NULL)
  )
})
