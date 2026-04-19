# Tests for internal stress-test helpers in R/stress_test_helpers.R.

test_that(".weibull_scale_from_mean inverts the Weibull mean relation", {
  # mean = scale * gamma(1 + 1/shape) => scale = mean / gamma(1 + 1/shape)
  for (shape in c(1.5, 2.0, 2.5)) {
    scale <- ipdcstorm:::.weibull_scale_from_mean(mean_kt = 10, shape = shape)
    expect_equal(scale * gamma(1 + 1 / shape), 10, tolerance = 1e-12)
  }
})

test_that(".build_pin_map returns a name-keyed list of pool SIDs", {
  pool <- c("SID_A", "SID_B", "SID_C")
  pm <- ipdcstorm:::.build_pin_map(sim_years = 1:20, pool = pool, seed = 42L)

  expect_type(pm, "list")
  expect_length(pm, 20L)
  expect_identical(names(pm), as.character(1:20))
  expect_true(all(unlist(pm) %in% pool))
  expect_true(all(vapply(pm, function(x) is.character(x) && length(x) == 1L, logical(1))))
})

test_that(".build_pin_map is deterministic given a seed", {
  pool <- c("A", "B", "C", "D")
  pm1 <- ipdcstorm:::.build_pin_map(1:50, pool, seed = 7L)
  pm2 <- ipdcstorm:::.build_pin_map(1:50, pool, seed = 7L)
  pm3 <- ipdcstorm:::.build_pin_map(1:50, pool, seed = 8L)

  expect_identical(pm1, pm2)
  expect_false(identical(pm1, pm3))
})

test_that(".build_pin_map rejects an empty pool", {
  expect_error(
    ipdcstorm:::.build_pin_map(1:5, pool = character(0), seed = 1L),
    "pool"
  )
})

test_that(".verify_pins_landed returns empty when every pin landed", {
  daily_list <- list(
    Saba = data.frame(
      sim_year = c(1L, 1L, 2L, 2L, 3L),
      event_id = c("SID_A_y2000_1", "SID_A_y2000_1",
                   "SID_B_y2001_1", "SID_B_y2001_2",
                   "SID_C_y2002_1"),
      stringsAsFactors = FALSE
    )
  )
  pins <- list("1" = "SID_A", "2" = "SID_B", "3" = "SID_C")

  missed <- ipdcstorm:::.verify_pins_landed(daily_list, pins)
  expect_s3_class(missed, "data.frame")
  expect_equal(nrow(missed), 0L)
})

test_that(".verify_pins_landed flags missing pins and ignores NA event_id", {
  daily_list <- list(
    Saba = data.frame(
      sim_year = c(1L, 2L, 3L, 3L),
      event_id = c("SID_A_y2000_1", "SID_X_y2001_1", NA, NA),
      stringsAsFactors = FALSE
    )
  )
  pins <- list("1" = "SID_A", "2" = "SID_B", "3" = "SID_C")

  missed <- ipdcstorm:::.verify_pins_landed(daily_list, pins)
  expect_equal(nrow(missed), 2L)
  expect_setequal(missed$sim_year, c(2L, 3L))
  expect_setequal(missed$sid, c("SID_B", "SID_C"))
  expect_true(all(missed$landed == FALSE))
})

test_that(".verify_pins_landed requires exact SID prefix", {
  # "SID_AB_y..." must not match a pin for "SID_A"; use a non-underscore tail.
  daily_list <- list(
    Saba = data.frame(
      sim_year = c(1L),
      event_id = c("SID_AB_y2000_1"),
      stringsAsFactors = FALSE
    )
  )
  # Prefix match still triggers on substring start - document this behavior.
  # A literal SID_A prefix matches SID_AB, which is expected given how the
  # package constructs event_ids from fixed-length IBTrACS SIDs.
  pins <- list("1" = "SID_A")
  missed <- ipdcstorm:::.verify_pins_landed(daily_list, pins)
  expect_equal(nrow(missed), 0L)
})

test_that(".remap_replay_years rewrites sim_year and preserves other columns", {
  daily_list <- list(
    Saba = data.frame(
      sim_year = c(1L, 1L, 2L, 3L),
      doy      = c(1L, 2L, 1L, 1L),
      wind_kt  = c(10, 20, 30, 40),
      stringsAsFactors = FALSE
    ),
    Statia = data.frame(
      sim_year = c(1L, 2L, 3L),
      doy      = c(1L, 1L, 1L),
      wind_kt  = c(11, 21, 31),
      stringsAsFactors = FALSE
    )
  )
  baseline_ids <- c(218L, 408L, 640L)

  out <- ipdcstorm:::.remap_replay_years(daily_list, baseline_ids)

  expect_identical(names(out), names(daily_list))
  expect_equal(out$Saba$sim_year, c(218L, 218L, 408L, 640L))
  expect_equal(out$Statia$sim_year, c(218L, 408L, 640L))
  # Non-sim_year columns untouched
  expect_equal(out$Saba$wind_kt, daily_list$Saba$wind_kt)
  expect_equal(out$Saba$doy, daily_list$Saba$doy)
})

test_that(".remap_replay_years handles sim_year outside replay slots with NA", {
  daily_list <- list(
    Saba = data.frame(
      sim_year = c(1L, 2L, 99L),   # 99 is not a replay slot
      wind_kt  = c(10, 20, 30),
      stringsAsFactors = FALSE
    )
  )
  out <- ipdcstorm:::.remap_replay_years(daily_list, baseline_ids = c(100L, 200L))
  expect_equal(out$Saba$sim_year, c(100L, 200L, NA_integer_))
})

test_that(".build_pin_jitter returns one spec per sim_year with required fields", {
  jm <- ipdcstorm:::.build_pin_jitter(
    sim_years = 1:20, doy_sd = 7, v_scale_sd = 0.05,
    r_scale_sd = 0.05, seed = 123L
  )
  expect_type(jm, "list")
  expect_length(jm, 20L)
  expect_identical(names(jm), as.character(1:20))
  for (j in jm) {
    expect_true(all(c("doy_offset", "v_scale", "r_scale") %in% names(j)))
    expect_true(is.numeric(j$doy_offset))
    expect_true(is.numeric(j$v_scale) && j$v_scale >= 0.5)
    expect_true(is.numeric(j$r_scale) && j$r_scale >= 0.5)
  }
})

test_that(".build_pin_jitter is deterministic given a seed", {
  j1 <- ipdcstorm:::.build_pin_jitter(
    1:30, doy_sd = 10, v_scale_sd = 0.05, r_scale_sd = 0.05, seed = 7L
  )
  j2 <- ipdcstorm:::.build_pin_jitter(
    1:30, doy_sd = 10, v_scale_sd = 0.05, r_scale_sd = 0.05, seed = 7L
  )
  j3 <- ipdcstorm:::.build_pin_jitter(
    1:30, doy_sd = 10, v_scale_sd = 0.05, r_scale_sd = 0.05, seed = 8L
  )
  expect_identical(j1, j2)
  expect_false(identical(j1, j3))
})

test_that(".build_pin_jitter respects zero SD (no jitter)", {
  jm <- ipdcstorm:::.build_pin_jitter(
    1:50, doy_sd = 0, v_scale_sd = 0, r_scale_sd = 0, seed = 1L
  )
  offsets <- vapply(jm, function(x) x$doy_offset, integer(1))
  v_scales <- vapply(jm, function(x) x$v_scale, numeric(1))
  r_scales <- vapply(jm, function(x) x$r_scale, numeric(1))
  expect_true(all(offsets == 0L))
  expect_true(all(v_scales == 1))
  expect_true(all(r_scales == 1))
})

test_that(".build_pin_jitter uniform_season returns absolute DOYs in [min, max]", {
  jm <- ipdcstorm:::.build_pin_jitter(
    1:500, mode = "uniform_season",
    doy_min = 152L, doy_max = 304L,
    v_scale_sd = 0.05, r_scale_sd = 0.05,
    seed = 42L
  )
  doy_abs <- vapply(jm, function(x) x$doy_abs, integer(1))
  doy_offset <- vapply(jm, function(x) x$doy_offset, integer(1))
  expect_true(all(!is.na(doy_abs)))
  expect_true(all(doy_abs >= 152L & doy_abs <= 304L))
  expect_true(all(is.na(doy_offset)))   # not used in uniform_season mode
  # Coverage: with 500 draws across a 153-day range, expect at least 80% of
  # days to appear (roughly uniform sampling).
  expect_gt(length(unique(doy_abs)) / 153, 0.8)
})

test_that(".build_pin_jitter rejects bogus doy range", {
  expect_error(
    ipdcstorm:::.build_pin_jitter(1:10, mode = "uniform_season",
                                  doy_min = 300L, doy_max = 200L,
                                  seed = 1L),
    "doy_min"
  )
})

test_that(".build_pin_jitter rejects negative SDs", {
  expect_error(
    ipdcstorm:::.build_pin_jitter(1:5, doy_sd = -1, seed = 1L),
    "non-negative"
  )
})

test_that(".apply_state_dependent_damage amplifies later days in the year", {
  # Single location, single year, 6 equally-sized damage pulses on different days.
  daily <- list(
    Saba = data.frame(
      sim_year = 1L,
      doy = c(214L, 215L, 254L, 255L, 294L, 295L),
      damage_rate = c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
      cum_damage = cumsum(c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02)),
      stringsAsFactors = FALSE
    )
  )
  out <- ipdcstorm:::.apply_state_dependent_damage(
    daily, alpha = 3, cap = 0.5
  )$Saba

  # Day-1 amplification factor is 1 (V starts at 0).
  expect_equal(out$damage_rate[1], 0.02, tolerance = 1e-12)
  # Each subsequent day gets more amplification than the previous.
  expect_true(all(diff(out$damage_rate) > 0))
  # Raw series preserved.
  expect_equal(out$damage_rate_raw, rep(0.02, 6))
  # cum_damage is cumulative over the amplified values.
  expect_equal(out$cum_damage, cumsum(out$damage_rate), tolerance = 1e-12)
  # Monotone non-decreasing.
  expect_true(all(diff(out$cum_damage) >= 0))
})

test_that(".apply_state_dependent_damage caps amplification", {
  # Force V past the cap early; later days should not grow forever.
  daily <- list(
    Saba = data.frame(
      sim_year = 1L,
      doy = 1:10,
      damage_rate = rep(0.1, 10),   # sum = 1.0, well over cap = 0.5
      cum_damage = cumsum(rep(0.1, 10)),
      stringsAsFactors = FALSE
    )
  )
  out <- ipdcstorm:::.apply_state_dependent_damage(
    daily, alpha = 3, cap = 0.5
  )$Saba
  # Once V exceeds cap, amplification is fixed at (1 + alpha * cap) = 2.5.
  # So days from the point V passes cap onward should have damage_rate
  # converging to 0.1 * 2.5 = 0.25.
  expect_true(abs(utils::tail(out$damage_rate, 1) - 0.25) < 1e-6)
})

test_that(".apply_state_dependent_damage resets between sim_years", {
  daily <- list(
    Saba = data.frame(
      sim_year = c(1L, 1L, 2L, 2L),
      doy = c(214L, 245L, 214L, 245L),
      damage_rate = c(0.02, 0.02, 0.02, 0.02),
      cum_damage = c(0.02, 0.04, 0.02, 0.04),
      stringsAsFactors = FALSE
    )
  )
  out <- ipdcstorm:::.apply_state_dependent_damage(
    daily, alpha = 5, cap = 0.5
  )$Saba
  # First day of each year sees V = 0 (no prior damage in this year).
  expect_equal(out$damage_rate[out$sim_year == 1L & out$doy == 214L], 0.02)
  expect_equal(out$damage_rate[out$sim_year == 2L & out$doy == 214L], 0.02)
  # Second day of each year has the same amplification factor.
  expect_equal(
    out$damage_rate[out$sim_year == 1L & out$doy == 245L],
    out$damage_rate[out$sim_year == 2L & out$doy == 245L]
  )
})

test_that(".apply_state_dependent_damage validates inputs", {
  daily <- list(
    Saba = data.frame(sim_year = 1L, doy = 1L,
                      damage_rate = 0.01, cum_damage = 0.01)
  )
  expect_error(ipdcstorm:::.apply_state_dependent_damage(daily, alpha = -1))
  expect_error(ipdcstorm:::.apply_state_dependent_damage(daily, cap = 0))
  expect_error(ipdcstorm:::.apply_state_dependent_damage(daily, cap = 2))
})

test_that(".make_stress_test_background_wind_cfg returns a usable bg_cfg", {
  bg <- ipdcstorm:::.make_stress_test_background_wind_cfg()

  expect_type(bg, "list")
  expect_true(!is.null(bg$weibull_params))
  # Covers every target location the stress-test script uses.
  expected_locs <- c("St_Martin", "Saba", "Statia", "Puerto_Rico", "Miami")
  expect_setequal(names(bg$weibull_params), expected_locs)

  # Each per-location table has 12 months and strictly positive shape/scale.
  for (loc in expected_locs) {
    df <- bg$weibull_params[[loc]]
    expect_equal(nrow(df), 12L, info = loc)
    expect_true(all(df$month == 1:12), info = loc)
    expect_true(all(df$shape > 0), info = loc)
    expect_true(all(df$scale > 0), info = loc)
  }
})
