test_that("observed R34 uses a 1.5x Holland outer cutoff", {
  R34_km <- c(80, 120)
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_is_fallback = FALSE
  )

  expect_equal(R_outer_km / R34_km, c(1.5, 1.5))
})

test_that("fallback R34 uses a 1.8x Holland outer cutoff", {
  R34_km <- c(80, 120)
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_is_fallback = TRUE
  )

  expect_equal(R_outer_km / R34_km, c(1.8, 1.8))
})

test_that("observed R34 decays immediately beyond the 1.5x cutoff", {
  Vmax_kt <- 50
  R34_km <- 100
  RMW_km <- 140
  r_km <- 160
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_km,
    R34_is_fallback = FALSE
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

test_that("fallback R34 decays immediately beyond the 1.8x cutoff", {
  Vmax_kt <- 34
  RMW_km <- 200
  lat <- 25
  R34_fallback_km <- ipdcstorm:::estimate_R34_climo(Vmax_kt, lat = lat)
  R_outer_km <- ipdcstorm:::.resolve_holland_outer_cutoff_km(
    R34_km = R34_fallback_km,
    R34_is_fallback = TRUE
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
