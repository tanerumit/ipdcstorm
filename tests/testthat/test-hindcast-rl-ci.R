test_that("observed RL CI returns finite ordered bounds for moderate sample", {
  set.seed(101)
  x <- ifelse(stats::runif(160) < 0.2, 0, ipdcstorm:::.rgev(160, mu = 52, sigma = 11, xi = 0.08))
  ci <- ipdcstorm:::.compute_obs_return_level_ci(
    annual_max = x,
    return_periods = c(5, 10, 25, 50),
    xi_bounds = c(-0.3, 0.4),
    n_boot_fallback = 220L,
    seed = 123
  )

  expect_true(all(is.finite(ci$sim_ci_lo)))
  expect_true(all(is.finite(ci$sim_ci_hi)))
  expect_true(all(ci$sim_ci_lo <= ci$sim_ci_hi))
  expect_true(all(ci$ci_method %in% c("delta", "bootstrap")))
})

test_that("small-to-moderate sample still yields CI via delta or bootstrap", {
  set.seed(202)
  x <- ifelse(stats::runif(20) < 0.1, 0, ipdcstorm:::.rgev(20, mu = 48, sigma = 10, xi = 0.05))
  ci <- ipdcstorm:::.compute_obs_return_level_ci(
    annual_max = x,
    return_periods = c(5, 10, 25, 50),
    xi_bounds = c(-0.3, 0.4),
    n_boot_fallback = 220L,
    seed = 456
  )

  expect_true(all(is.finite(ci$sim_ci_lo)))
  expect_true(all(is.finite(ci$sim_ci_hi)))
  expect_true(all(ci$sim_ci_lo <= ci$sim_ci_hi))
})

test_that("very small sample returns unavailable CI and NA classification", {
  x <- c(0, 0, 40, 42, 43)
  ci <- ipdcstorm:::.compute_obs_return_level_ci(
    annual_max = x,
    return_periods = c(5, 10, 25, 50),
    xi_bounds = c(-0.3, 0.4),
    n_boot_fallback = 220L,
    seed = 789
  )

  expect_true(all(is.na(ci$sim_ci_lo)))
  expect_true(all(is.na(ci$sim_ci_hi)))
  expect_true(all(ci$ci_method == "unavailable"))

  sim_rl <- rep(50, length(ci$return_period))
  outside_flag <- ifelse(
    is.finite(ci$sim_ci_lo) & is.finite(ci$sim_ci_hi),
    sim_rl < ci$sim_ci_lo | sim_rl > ci$sim_ci_hi,
    NA
  )
  expect_true(all(is.na(outside_flag)))
})

test_that("bootstrap fallback is used deterministically when delta path is blocked", {
  set.seed(303)
  x <- ifelse(stats::runif(80) < 0.15, 0, ipdcstorm:::.rgev(80, mu = 50, sigma = 9, xi = 0.06))

  ci1 <- ipdcstorm:::.compute_obs_return_level_ci(
    annual_max = x,
    return_periods = c(5, 10, 25, 50),
    xi_bounds = c(-0.3, 0.4),
    n_boot_fallback = 220L,
    seed = 999,
    force_bootstrap = TRUE
  )
  ci2 <- ipdcstorm:::.compute_obs_return_level_ci(
    annual_max = x,
    return_periods = c(5, 10, 25, 50),
    xi_bounds = c(-0.3, 0.4),
    n_boot_fallback = 220L,
    seed = 999,
    force_bootstrap = TRUE
  )

  expect_true(all(ci1$ci_method == "bootstrap"))
  expect_true(all(is.finite(ci1$sim_ci_lo)))
  expect_true(all(is.finite(ci1$sim_ci_hi)))
  expect_equal(ci1$sim_ci_lo, ci2$sim_ci_lo)
  expect_equal(ci1$sim_ci_hi, ci2$sim_ci_hi)
})
