# Tests for hazard query functions:
#   lookup_storm_id()        — filter historical events from run_hazard_model output
#   query_storm_track_years() — find sim_years where a given SID was sampled
#   query_impact_years()     — find sim_years exceeding a reference impact level
#   query_aftermath_impact() — compute post-event impact window metrics
#
# Run: Rscript -e "testthat::test_file('tests/testthat/test-hazard-query.R')"

# =============================================================================
# Helper: query_aftermath_impact iterates over every (location, sim_year) pair,
# including calm years where the anchor storm is absent.  That can trigger
# base-R "no non-missing arguments to max; returning -Inf" warnings internally.
# Suppress these at the call site so they don't pollute the test output.
# =============================================================================

.aftermath <- function(...) suppressWarnings(query_aftermath_impact(...))

# =============================================================================
# Shared fixtures
# =============================================================================

# Minimal out$events tibble (output of run_hazard_model)
.events_fixture <- function() {
  tibble::tibble(
    storm_id          = c("2017242N16333", "2005280N15305", "1985243N13315"),
    location          = c("Saba",          "Saba",          "Saba"),
    start_time        = as.POSIXct(
                          c("2017-09-06 12:00:00",
                            "2005-10-07 18:00:00",
                            "1985-08-31 00:00:00"), tz = "UTC"),
    peak_wind_kt      = c(82.0, 45.0, 71.0),
    storm_intensity_kt = c(125.0, 65.0, 90.0)
  )
}

# Minimal hazard model output list (only events needed for lookup_storm_id)
.out_fixture <- function() list(events = .events_fixture())

# Daily hazard tibble with event_id, cum_damage, damage_rate columns.
# Three years, one location ("Saba"):
#   Year 1: IRMA-like HUR (82 kt) + follow-up TS within aftermath window
#   Year 2: moderate HUR (67 kt), no aftermath activity
#   Year 3: calm year
.daily_fixture <- function() {
  make_year <- function(yr, events) {
    wind <- numeric(365L)
    eid  <- rep(NA_character_, 365L)

    for (ev in events) {
      wind[ev$days] <- ev$wind
      eid[ev$days]  <- paste0(ev$sid, "_y", yr, "_e", ev$en)
    }

    dmg_r <- pmax(0, (wind / 70)^3 * 0.01)
    tibble::tibble(
      location    = "Saba",
      sim_year    = as.integer(yr),
      date        = seq(as.Date("2000-01-01") + (yr - 1L) * 365L,
                        by = "day", length.out = 365L),
      wind_kt     = wind,
      event_id    = eid,
      damage_rate = dmg_r,
      cum_damage  = cumsum(dmg_r)
    )
  }

  dplyr::bind_rows(
    make_year(1L, list(
      list(sid = "2017242N16333", en = 1L, days = 250:252, wind = c(72, 82, 68)),
      list(sid = "2005280N15305", en = 2L, days = 270:272, wind = c(38, 44, 36))
    )),
    make_year(2L, list(
      list(sid = "2017242N16333", en = 1L, days = 280:282, wind = c(65, 75, 60))
    )),
    make_year(3L, list())
  )
}

# =============================================================================
# lookup_storm_id
# =============================================================================

test_that("lookup_storm_id returns tibble with required columns", {
  res <- lookup_storm_id(.out_fixture())

  expect_s3_class(res, "data.frame")
  expected <- c("storm_id", "location", "start_time", "peak_wind_kt", "storm_intensity_kt")
  expect_true(all(expected %in% names(res)))
})

test_that("lookup_storm_id with no filters returns all events", {
  res <- lookup_storm_id(.out_fixture())
  expect_equal(nrow(res), 3L)
})

test_that("lookup_storm_id year filter returns only matching years", {
  res <- lookup_storm_id(.out_fixture(), year = 2017L)

  expect_equal(nrow(res), 1L)
  expect_equal(res$storm_id, "2017242N16333")
})

test_that("lookup_storm_id location filter returns only that location", {
  out <- .out_fixture()
  out$events <- dplyr::bind_rows(
    out$events,
    dplyr::mutate(out$events[1L, ], location = "Statia")
  )
  res <- lookup_storm_id(out, location = "Saba")

  expect_true(all(res$location == "Saba"))
})

test_that("lookup_storm_id min_wind_kt filter returns only storms above threshold", {
  res <- lookup_storm_id(.out_fixture(), min_wind_kt = 65)

  expect_true(all(res$peak_wind_kt >= 65))
  expect_equal(nrow(res), 2L)
})

test_that("lookup_storm_id returns empty tibble with correct schema when no matches", {
  res <- suppressMessages(lookup_storm_id(.out_fixture(), year = 1900L))

  expect_equal(nrow(res), 0L)
  expect_true("storm_id" %in% names(res))
  expect_true("peak_wind_kt" %in% names(res))
})

test_that("lookup_storm_id stops when out$events is absent", {
  expect_error(lookup_storm_id(list()), regexp = "events")
})

# =============================================================================
# query_storm_track_years
# =============================================================================

test_that("query_storm_track_years returns location/sim_year tibble for matched SID", {
  res <- query_storm_track_years(.daily_fixture(), storm_id = "2017242N16333")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("location", "sim_year") %in% names(res)))
  # IRMA appears in years 1 and 2
  expect_setequal(res$sim_year, c(1L, 2L))
})

test_that("query_storm_track_years returns empty tibble for SID not in data", {
  res <- suppressMessages(
    query_storm_track_years(.daily_fixture(), storm_id = "1900999N00000")
  )

  expect_equal(nrow(res), 0L)
  expect_true(all(c("location", "sim_year") %in% names(res)))
})

test_that("query_storm_track_years stops on non-character storm_id", {
  expect_error(
    query_storm_track_years(.daily_fixture(), storm_id = 2017242),
    regexp = "storm_id"
  )
})

test_that("query_storm_track_years location filter restricts output", {
  # Two-location daily data
  daily_two <- dplyr::bind_rows(
    .daily_fixture(),
    dplyr::mutate(.daily_fixture(), location = "Statia")
  )
  res <- query_storm_track_years(daily_two,
                                 storm_id = "2017242N16333",
                                 location = "Saba")

  expect_true(all(res$location == "Saba"))
})

# =============================================================================
# query_impact_years
# =============================================================================

test_that("query_impact_years returns tibble with location, sim_year and metric column", {
  res <- query_impact_years(
    .daily_fixture(),
    storm_id  = "2017242N16333",
    metric    = "peak_wind_kt",
    threshold = 65
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("location", "sim_year") %in% names(res)))
  expect_true("peak_wind_kt" %in% names(res))
})

test_that("query_impact_years explicit threshold filters correctly", {
  res <- query_impact_years(
    .daily_fixture(),
    storm_id  = "2017242N16333",
    metric    = "peak_wind_kt",
    threshold = 65
  )

  # Only years with peak_wind >= 65 kt (years 1 and 2)
  expect_true(all(res$peak_wind_kt >= 65))
  expect_setequal(res$sim_year, c(1L, 2L))
})

test_that("query_impact_years returns empty tibble when threshold is very high", {
  res <- query_impact_years(
    .daily_fixture(),
    storm_id  = "2017242N16333",
    metric    = "peak_wind_kt",
    threshold = 200
  )

  expect_equal(nrow(res), 0L)
})

test_that("query_impact_years percentile gate restricts to upper quantile", {
  res_all  <- query_impact_years(
    .daily_fixture(), storm_id = "2017242N16333",
    metric    = "peak_wind_kt",
    threshold = 0           # very low explicit threshold
  )
  res_pct  <- query_impact_years(
    .daily_fixture(), storm_id = "2017242N16333",
    metric     = "peak_wind_kt",
    threshold  = 0,
    percentile = 0.6        # additional percentile gate
  )

  # Percentile gate should reduce or maintain (not increase) the selection
  expect_lte(nrow(res_pct), nrow(res_all))
})

test_that("query_impact_years min_threshold acts as absolute floor", {
  res <- query_impact_years(
    .daily_fixture(),
    storm_id      = "2017242N16333",
    metric        = "peak_wind_kt",
    threshold     = 0,
    min_threshold = 70       # floor at 70 kt
  )

  expect_true(all(res$peak_wind_kt >= 70))
})

test_that("query_impact_years stops on invalid storm_id type", {
  expect_error(
    query_impact_years(.daily_fixture(), storm_id = 123),
    regexp = "storm_id"
  )
})

# =============================================================================
# query_aftermath_impact
# =============================================================================

test_that("query_aftermath_impact returns tibble with required columns", {
  res <- .aftermath(.daily_fixture(), storm_id = "2017242N16333")

  expect_s3_class(res, "data.frame")
  expected <- c("location", "sim_year", "primary_event_id",
                "aftermath_peak_wind_kt", "aftermath_n_events",
                "aftermath_n_hur_days", "aftermath_cum_damage",
                "aftermath_max_damage_rate", "aftermath_rank")
  expect_true(all(expected %in% names(res)),
    info = paste("missing:", paste(setdiff(expected, names(res)), collapse = ", "))
  )
})

test_that("query_aftermath_impact aftermath metrics are non-negative", {
  res <- .aftermath(.daily_fixture(), storm_id = "2017242N16333")

  expect_true(all(res$aftermath_peak_wind_kt    >= 0))
  expect_true(all(res$aftermath_n_events         >= 0))
  expect_true(all(res$aftermath_n_hur_days        >= 0))
  expect_true(all(res$aftermath_cum_damage        >= 0))
  expect_true(all(res$aftermath_max_damage_rate   >= 0))
})

test_that("query_aftermath_impact aftermath_rank=1 has highest cum_damage within location", {
  res  <- .aftermath(.daily_fixture(), storm_id = "2017242N16333")
  top  <- res[res$aftermath_rank == 1L, ]

  expect_equal(nrow(top), 1L)
  expect_equal(top$aftermath_cum_damage, max(res$aftermath_cum_damage))
})

test_that("query_aftermath_impact year 1 captures follow-up storm within window", {
  # Year 1: IRMA ends day 252; follow-up TS at days 270-272 is within 60-day window
  res <- .aftermath(.daily_fixture(), storm_id = "2017242N16333", window_days = 60L)
  yr1 <- res[res$sim_year == 1L, ]

  expect_gte(yr1$aftermath_n_events,  1L)
  expect_gt( yr1$aftermath_peak_wind_kt, 0)
})

test_that("query_aftermath_impact stops on non-positive window_days", {
  expect_error(
    query_aftermath_impact(.daily_fixture(), window_days = 0L),
    regexp = "window_days|positive"
  )
})
