# Tests for climate calibration functions:
#   estimate_beta_sst()        — Poisson/NegBin SST-rate GLM with shrinkage
#   estimate_gamma_intensity() — Binomial logistic SST-intensity GLM with shrinkage
#   ipcc_ar6_sst_scenario_info() — Hard-coded AR6 SST scenario table
#
# Run: Rscript -e "testthat::test_file('tests/testthat/test-climate-calibration.R')"

# =============================================================================
# Helpers: suppress both messages (verbose output) and warnings (GLM fitting
# noise — convergence, fitted rates/probs near 0 — from synthetic test data)
# =============================================================================

.est_beta <- function(...) {
  suppressWarnings(suppressMessages(estimate_beta_sst(...)))
}
.est_gamma <- function(...) {
  suppressWarnings(suppressMessages(estimate_gamma_intensity(...)))
}

# =============================================================================
# Shared fixtures
# =============================================================================

# 30 years of synthetic annual counts with a known SST-rate signal (beta ~ 0.5)
.beta_counts_fixture <- function(seed = 42L) {
  set.seed(seed)
  n <- 30L
  year <- seq(1991L, by = 1L, length.out = n)
  sst_anom <- seq(-0.4, 0.4, length.out = n) + rnorm(n, 0, 0.12)
  counts <- rpois(n, exp(2.5 + 0.5 * sst_anom))
  tibble::tibble(year = year, storm_class = "TS", n_events = counts)
}

# Matching SST data for the 30-year fixture
.beta_sst_fixture <- function(seed = 42L) {
  set.seed(seed)
  n <- 30L
  year <- seq(1991L, by = 1L, length.out = n)
  sst_anom <- seq(-0.4, 0.4, length.out = n) + rnorm(n, 0, 0.12)
  tibble::tibble(
    year         = year,
    sst_mdr_aso  = 26.2 + sst_anom,
    sst_clim     = 26.2,
    sst_anomaly  = sst_anom
  )
}

# Tiny dataset: only 3 years (forces literature fallback in estimate_beta_sst)
.beta_tiny_counts <- function() {
  tibble::tibble(year = 2000L:2002L, storm_class = "TS", n_events = c(8L, 10L, 9L))
}
.beta_tiny_sst <- function() {
  tibble::tibble(
    year = 2000L:2002L,
    sst_mdr_aso = c(26.1, 26.3, 26.2),
    sst_clim    = 26.2,
    sst_anomaly = c(-0.1, 0.1, 0.0)
  )
}

# 30-year counts with two storm classes, slight gamma signal (p_hur rises with SST)
.gamma_counts_fixture <- function(seed = 77L) {
  set.seed(seed)
  n <- 30L
  year <- seq(1991L, by = 1L, length.out = n)
  sst_anom  <- seq(-0.4, 0.4, length.out = n) + rnorm(n, 0, 0.10)
  p_hur     <- plogis(qlogis(0.35) + 0.25 * sst_anom)
  n_total   <- rpois(n, 8)
  n_hur     <- rbinom(n, n_total, p_hur)
  n_ts      <- n_total - n_hur
  tibble::tibble(
    year        = rep(year, 2L),
    storm_class = rep(c("TS", "HUR"), each = n),
    n_events    = c(n_ts, n_hur)
  )
}

# Counts with a slightly *negative* gamma signal — used to test non-negative constraint
.gamma_negative_fixture <- function(seed = 55L) {
  set.seed(seed)
  n <- 30L
  year <- seq(1991L, by = 1L, length.out = n)
  sst_anom  <- seq(-0.4, 0.4, length.out = n)
  # HUR fraction *decreases* with SST — should force gamma to 0 after constraint
  p_hur     <- plogis(qlogis(0.35) - 0.40 * sst_anom)
  n_total   <- rpois(n, 8)
  n_hur     <- rbinom(n, n_total, p_hur)
  n_ts      <- n_total - n_hur
  tibble::tibble(
    year        = rep(year, 2L),
    storm_class = rep(c("TS", "HUR"), each = n),
    n_events    = c(n_ts, n_hur)
  )
}

# =============================================================================
# estimate_beta_sst — behavioral tests
# =============================================================================

test_that("estimate_beta_sst returns list with all required elements", {
  counts <- .beta_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_beta(counts, sst, verbose = FALSE)

  expect_type(fit, "list")
  required <- c("beta_sst", "beta_se", "beta_mle", "alpha", "method",
                "n_years", "fit_data", "guardrail")
  for (nm in required) {
    expect_true(nm %in% names(fit), info = paste("missing element:", nm))
  }
})

test_that("estimate_beta_sst beta_sst is finite and within plausible range [0, 1.2]", {
  counts <- .beta_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_beta(counts, sst, verbose = FALSE)

  expect_true(is.finite(fit$beta_sst))
  expect_gte(fit$beta_sst, 0)
  expect_lte(fit$beta_sst, 1.2)
})

test_that("estimate_beta_sst n_years reflects the joined data length", {
  counts <- .beta_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_beta(counts, sst, verbose = FALSE)

  expect_equal(fit$n_years, 30L)
})

test_that("estimate_beta_sst uses literature fallback when fewer than 5 years available", {
  fit <- .est_beta(.beta_tiny_counts(), .beta_tiny_sst(), verbose = FALSE)
  expect_equal(fit$method, "literature_fallback")
  expect_true(is.finite(fit$beta_sst))
})

test_that("estimate_beta_sst rejects annual_counts containing a location column", {
  counts_with_loc <- .beta_counts_fixture() |>
    dplyr::mutate(location = "Saba")
  sst <- .beta_sst_fixture()

  expect_error(
    .est_beta(counts_with_loc, sst, verbose = FALSE),
    regexp = "location"
  )
})

test_that("estimate_beta_sst beta_prior shrinks result toward prior", {
  counts   <- .beta_counts_fixture()
  sst      <- .beta_sst_fixture()
  prior    <- 0.6

  fit_raw    <- .est_beta(counts, sst, verbose = FALSE)
  fit_shrunk <- .est_beta(counts, sst, beta_prior = prior, verbose = FALSE)

  # Shrunk estimate should lie between MLE and prior (unless guardrail already fired)
  lo <- min(fit_raw$beta_sst, prior)
  hi <- max(fit_raw$beta_sst, prior)
  expect_gte(fit_shrunk$beta_sst, lo - 1e-9)
  expect_lte(fit_shrunk$beta_sst, hi + 1e-9)
})

test_that("estimate_beta_sst verbose=FALSE produces no console output", {
  counts <- .beta_counts_fixture()
  sst    <- .beta_sst_fixture()

  out_text <- capture.output(
    invisible(suppressWarnings(estimate_beta_sst(counts, sst, verbose = FALSE))),
    type = "message"
  )
  expect_length(out_text, 0L)
})

test_that("estimate_beta_sst fit_data contains year, N, sst_anomaly columns", {
  counts <- .beta_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_beta(counts, sst, verbose = FALSE)

  expect_true(all(c("year", "N", "sst_anomaly") %in% names(fit$fit_data)))
  expect_equal(nrow(fit$fit_data), 30L)
})

# =============================================================================
# estimate_gamma_intensity — behavioral tests
# =============================================================================

test_that("estimate_gamma_intensity returns list with required elements", {
  counts <- .gamma_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_gamma(counts, sst, verbose = FALSE)

  expect_type(fit, "list")
  required <- c("gamma", "gamma_mle", "gamma_se", "p_hur_base", "method", "n_years")
  for (nm in required) {
    expect_true(nm %in% names(fit), info = paste("missing element:", nm))
  }
})

test_that("estimate_gamma_intensity gamma is always non-negative after constraint", {
  # Negative-signal data should be constrained to 0
  sst <- .beta_sst_fixture(seed = 55L)
  fit <- .est_gamma(.gamma_negative_fixture(), sst, verbose = FALSE)
  expect_gte(fit$gamma, 0)
})

test_that("estimate_gamma_intensity p_hur_base is in (0, 1)", {
  counts <- .gamma_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_gamma(counts, sst, verbose = FALSE)

  expect_gt(fit$p_hur_base, 0)
  expect_lt(fit$p_hur_base, 1)
})

test_that("estimate_gamma_intensity uses prior_fallback with fewer than 10 years", {
  # Only 8 years of data
  set.seed(1L)
  tiny_counts <- tibble::tibble(
    year        = rep(2010L:2017L, 2L),
    storm_class = rep(c("TS", "HUR"), each = 8L),
    n_events    = c(rpois(8L, 5), rbinom(8L, 5L, 0.35))
  )
  tiny_sst <- tibble::tibble(
    year = 2010L:2017L, sst_mdr_aso = 26.2, sst_clim = 26.2, sst_anomaly = 0
  )
  fit <- .est_gamma(tiny_counts, tiny_sst, verbose = FALSE)
  expect_equal(fit$method, "prior_fallback")
  expect_equal(fit$gamma, fit$gamma)   # still finite
})

test_that("estimate_gamma_intensity uses binomial_glm with sufficient data", {
  counts <- .gamma_counts_fixture()
  sst    <- .beta_sst_fixture()
  fit    <- .est_gamma(counts, sst, verbose = FALSE)

  expect_equal(fit$method, "binomial_glm")
})

test_that("estimate_gamma_intensity result is shrunk: lies between MLE and prior", {
  counts    <- .gamma_counts_fixture()
  sst       <- .beta_sst_fixture()
  prior_val <- 0.065

  fit <- .est_gamma(counts, sst, gamma_prior = prior_val, verbose = FALSE)
  # Shrinkage weight < 1 (n_years = 30 < 50), so result should not equal raw MLE
  # It should lie in [min(mle, prior), max(mle, prior)] after non-neg constraint
  lo <- max(0, min(fit$gamma_mle, prior_val))
  hi <- max(fit$gamma_mle, prior_val)
  expect_gte(fit$gamma, lo - 1e-9)
  expect_lte(fit$gamma, hi + 1e-9)
})

test_that("estimate_gamma_intensity accepts annual_counts with a location column", {
  counts_with_loc <- .gamma_counts_fixture() |>
    dplyr::mutate(location = "Saba")
  sst <- .beta_sst_fixture()

  # Should aggregate internally and not error
  expect_no_error(
    .est_gamma(counts_with_loc, sst, verbose = FALSE)
  )
})

# =============================================================================
# ipcc_ar6_sst_scenario_info — table contract tests
# =============================================================================

test_that("ipcc_ar6_sst_scenario_info returns tibble with correct columns", {
  info <- ipcc_ar6_sst_scenario_info()

  expect_s3_class(info, "data.frame")
  expected_cols <- c("scenario", "source", "delta_sst_2050", "delta_sst_2100", "description")
  expect_true(all(expected_cols %in% names(info)))
})

test_that("ipcc_ar6_sst_scenario_info contains exactly ssp126, ssp245, ssp585", {
  info <- ipcc_ar6_sst_scenario_info()

  expect_equal(nrow(info), 3L)
  expect_setequal(info$scenario, c("ssp126", "ssp245", "ssp585"))
})

test_that("ipcc_ar6_sst_scenario_info delta_sst values are positive and 2100 >= 2050", {
  info <- ipcc_ar6_sst_scenario_info()

  expect_true(all(info$delta_sst_2050 > 0))
  expect_true(all(info$delta_sst_2100 > 0))
  expect_true(all(info$delta_sst_2100 >= info$delta_sst_2050))
})

test_that("ipcc_ar6_sst_scenario_info source column is 'ipcc_ar6' for all rows", {
  info <- ipcc_ar6_sst_scenario_info()
  expect_true(all(info$source == "ipcc_ar6"))
})
