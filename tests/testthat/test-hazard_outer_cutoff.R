test_that("observed R34 uses a 1.5x Holland outer cutoff", {
  old_opt <- options(ipdcstorm.wind_field_mode = NULL)
  on.exit(options(old_opt), add = TRUE)

  R34_km <- c(80, 120)
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_source = "observed"
  )

  expect_equal(R_outer_km / R34_km, c(1.5, 1.5))
})

test_that("default ordinary runs stay on the validated legacy wind path", {
  old_opt <- options(
    ipdcstorm.wind_field_mode = NULL,
    ipdcstorm.hindcast_sampler_mode = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  expect_equal(getOption("ipdcstorm.wind_field_mode", "legacy"), "legacy")
  expect_equal(ipdcstorm:::.hindcast_sampler_mode(), "legacy")
  expect_equal(
    ipdcstorm:::.holland_outer_cutoff_multipliers(),
    list(observed = 1.5, partial = 1.5, climo = 1.5)
  )
})

test_that("partial and climatology R34 default to the legacy 1.5x Holland outer cutoff", {
  R34_km <- c(80, 120)
  R_outer_partial <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_source = "partial"
  )
  R_outer_climo <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_source = "climo"
  )

  expect_equal(R_outer_partial / R34_km, c(1.5, 1.5))
  expect_equal(R_outer_climo / R34_km, c(1.5, 1.5))
})

test_that("observed R34 decays immediately beyond the 1.5x cutoff", {
  Vmax_kt <- 50
  R34_km <- 100
  RMW_km <- 140
  r_km <- 160
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_source = "observed"
  )

  wind_no_decay <- ipdcstorm:::.estimate_site_wind_holland(
    Vmax_kt = Vmax_kt,
    r_km = r_km,
    R34_km = 200,
    RMW_km = RMW_km,
    lat = 18
  )
  wind_with_decay <- ipdcstorm:::.estimate_site_wind_holland(
    Vmax_kt = Vmax_kt,
    r_km = r_km,
    R34_km = R34_km,
    RMW_km = RMW_km,
    lat = 18
  )
  expected_decay <- exp(-2 * (r_km - R_outer_km) / R_outer_km)

  expect_gt(r_km, R_outer_km)
  expect_equal(wind_with_decay, wind_no_decay * expected_decay, tolerance = 1e-10)
})

test_that("climatology R34 decays immediately beyond the legacy 1.5x cutoff", {
  Vmax_kt <- 34
  RMW_km <- 200
  lat <- 25
  R34_fallback_km <- ipdcstorm:::estimate_R34_climo(Vmax_kt, lat = lat)
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_fallback_km,
    R34_source = "climo"
  )
  delta_km <- 3
  r0_km <- R_outer_km
  r1_km <- ceiling(R_outer_km + delta_km)

  # Compute wind at the cutoff (no decay applied at r == R_outer),
  # then verify the beyond-cutoff exponential taper is applied.
  wind_at_cutoff <- ipdcstorm:::.estimate_site_wind_holland(
    Vmax_kt = Vmax_kt,
    r_km = r0_km,
    R34_km = NA_real_,
    RMW_km = RMW_km,
    lat = lat
  )
  wind_with_decay <- ipdcstorm:::.estimate_site_wind_holland(
    Vmax_kt = Vmax_kt,
    r_km = r1_km,
    R34_km = NA_real_,
    RMW_km = RMW_km,
    lat = lat
  )
  expected_decay <- exp(-2 * (r1_km - r0_km) / R_outer_km)

  expect_gt(r1_km, R_outer_km)
  expect_equal(wind_with_decay, wind_at_cutoff * expected_decay, tolerance = 0.5)
})
