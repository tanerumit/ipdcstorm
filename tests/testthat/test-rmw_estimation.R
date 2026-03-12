rmw_fixture_track <- function(rmw_km = NA_real_,
                              wind_kt = 90,
                              lat = 18,
                              r34_mean_km = 240,
                              r50_mean_km = 150,
                              r64_mean_km = 100) {
  to_nm <- function(x) {
    ifelse(is.finite(x), x / 1.852, NA_real_)
  }

  tibble::tibble(
    SID = "AL012020",
    iso_time = as.POSIXct("2020-08-01 00:00:00", tz = "UTC"),
    lat = lat,
    lon = -63,
    dist_km = 40,
    wind_kt = wind_kt,
    rmw_km = rmw_km,
    r34_ne_nm = to_nm(r34_mean_km),
    r34_se_nm = to_nm(r34_mean_km),
    r34_sw_nm = to_nm(r34_mean_km),
    r34_nw_nm = to_nm(r34_mean_km),
    r50_ne_nm = to_nm(r50_mean_km),
    r50_se_nm = to_nm(r50_mean_km),
    r50_sw_nm = to_nm(r50_mean_km),
    r50_nw_nm = to_nm(r50_mean_km),
    r64_ne_nm = to_nm(r64_mean_km),
    r64_se_nm = to_nm(r64_mean_km),
    r64_sw_nm = to_nm(r64_mean_km),
    r64_nw_nm = to_nm(r64_mean_km)
  )
}

test_that("observed rmw overrides inferred values when valid", {
  track <- rmw_fixture_track(rmw_km = 42)
  out <- compute_site_winds_full(track, target_lat = 18.5, target_lon = -62.5)

  expect_equal(out$RMW_km, 42)
})

test_that("invalid observed rmw is ignored and inference is used", {
  inferred <- ipdcstorm:::.cap_inferred_rmw_km(0.6517550 * 100, 90)
  resolved <- ipdcstorm:::.resolve_trackpoint_rmw_km(
    rmw_obs_km = c(NA_real_, 0, 200),
    R64_mean_km = c(100, 100, 100),
    R50_mean_km = c(150, 150, 150),
    R34_mean_km = c(240, 240, 240),
    Vmax_kt = c(90, 90, 90),
    lat = c(18, 18, 18)
  )

  expect_equal(resolved, rep(inferred, 3))
})

test_that("radii to rmw mapping uses fixed calibrated coefficients", {
  mapped <- ipdcstorm:::.estimate_rmw_from_mean_radii(
    R64_mean_km = c(NA_real_, 120, NA_real_),
    R50_mean_km = c(150, 200, NA_real_),
    R34_mean_km = c(300, 300, 300)
  )

  expect_equal(
    mapped,
    c(0.6676334 * 150, 0.6517550 * 120, 0.4106665 * 300),
    tolerance = 1e-8
  )
})

test_that("knaff fallback respects shared min and max caps", {
  rmw_km <- estimate_RMW_knaff(
    Vmax_kt = c(30, 160),
    lat = c(80, 35)
  )

  expect_equal(rmw_km, c(140, 50), tolerance = 1e-8)
})

test_that("rmw resolution is deterministic", {
  args <- list(
    rmw_obs_km = c(NA_real_, 38),
    R64_mean_km = c(100, 100),
    R50_mean_km = c(150, 150),
    R34_mean_km = c(240, 240),
    Vmax_kt = c(95, 95),
    lat = c(18, 18)
  )

  first <- do.call(ipdcstorm:::.resolve_trackpoint_rmw_km, args)
  second <- do.call(ipdcstorm:::.resolve_trackpoint_rmw_km, args)

  expect_equal(first, second)
})

test_that("internal option can disable R34 calibration for diagnostics", {
  args <- list(
    Vmax_kt = 80,
    r_km = 280,
    R34_km = 220,
    RMW_km = 35,
    lat = 18
  )

  default_val <- do.call(ipdcstorm:::.estimate_site_wind_holland, args)

  old_opt <- getOption("ipdcstorm.disable_r34_calibration")
  on.exit(options(ipdcstorm.disable_r34_calibration = old_opt), add = TRUE)
  options(ipdcstorm.disable_r34_calibration = TRUE)

  no_cal_val <- do.call(ipdcstorm:::.estimate_site_wind_holland, args)

  expect_lt(no_cal_val, default_val)
})

test_that("missing R34 falls back to climatology and remains deterministic", {
  args <- list(
    Vmax_kt = 80,
    r_km = 280,
    R34_km = 220,
    RMW_km = 35,
    lat = 18
  )

  observed_val <- do.call(ipdcstorm:::.estimate_site_wind_holland, args)
  fallback_args <- utils::modifyList(args, list(R34_km = NA_real_))
  fallback_first <- do.call(ipdcstorm:::.estimate_site_wind_holland, fallback_args)
  fallback_second <- do.call(ipdcstorm:::.estimate_site_wind_holland, fallback_args)

  expect_true(is.finite(fallback_first))
  expect_identical(fallback_first, fallback_second)
  expect_lte(fallback_first, observed_val)
})

test_that("climatological R34 fallback leaves near-core winds unchanged in the tested regime", {
  args <- list(
    Vmax_kt = 130,
    r_km = 38.5,
    R34_km = 220,
    RMW_km = 35,
    lat = 18
  )

  observed_val <- do.call(ipdcstorm:::.estimate_site_wind_holland, args)
  fallback_val <- do.call(
    ipdcstorm:::.estimate_site_wind_holland,
    utils::modifyList(args, list(R34_km = NA_real_))
  )

  expect_equal(fallback_val, observed_val, tolerance = 1e-10)
})

test_that("compute_site_winds_full records current R34 fallback diagnostics", {
  to_nm <- function(x) ifelse(is.finite(x), x / 1.852, NA_real_)
  track <- tibble::tibble(
    SID = rep("AL012020", 3),
    iso_time = as.POSIXct(
      c("2020-08-01 00:00:00", "2020-08-01 06:00:00", "2020-08-01 12:00:00"),
      tz = "UTC"
    ),
    lat = c(18, 18, 18),
    lon = c(-63, -63, -63),
    dist_km = c(40, 40, 40),
    wind_kt = c(90, 90, 90),
    rmw_km = c(35, 35, 35),
    storm_speed_kt = c(10, 10, 10),
    r34_ne_nm = c(to_nm(180), NA_real_, NA_real_),
    r34_se_nm = c(to_nm(180), to_nm(220), NA_real_),
    r34_sw_nm = c(to_nm(180), to_nm(200), NA_real_),
    r34_nw_nm = c(to_nm(180), NA_real_, NA_real_),
    r50_ne_nm = c(to_nm(120), to_nm(120), to_nm(120)),
    r50_se_nm = c(to_nm(120), to_nm(120), to_nm(120)),
    r50_sw_nm = c(to_nm(120), to_nm(120), to_nm(120)),
    r50_nw_nm = c(to_nm(120), to_nm(120), to_nm(120)),
    r64_ne_nm = c(to_nm(80), to_nm(80), to_nm(80)),
    r64_se_nm = c(to_nm(80), to_nm(80), to_nm(80)),
    r64_sw_nm = c(to_nm(80), to_nm(80), to_nm(80)),
    r64_nw_nm = c(to_nm(80), to_nm(80), to_nm(80))
  )

  out <- compute_site_winds_full(track, target_lat = 18.5, target_lon = -62.5)

  expect_equal(out$R34_missing, c(FALSE, FALSE, TRUE))
  expect_equal(out$R34_is_climo, c(FALSE, FALSE, TRUE))
  expect_true(all(is.finite(out$R34_eff_km[1:2])))
  expect_true(is.finite(out$R34_eff_km[3]))
  expect_equal(out$R34_eff_km[1], 180, tolerance = 1e-8)
  expect_equal(out$R34_eff_km[2], 210, tolerance = 1e-8)
  expect_gt(out$R34_eff_km[3], 0)
})
