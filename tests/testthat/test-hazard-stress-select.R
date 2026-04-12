# Tests for stress-year selection pipeline:
#   compute_stress_year_metrics() — annual metrics from daily hazard tibble
#   aggregate_stress_metrics()    — composite score via rank normalisation + weights
#   select_stress_years()         — stratified / diverse / top selection of k years
#
# Run: Rscript -e "testthat::test_file('tests/testthat/test-hazard-stress-select.R')"

# =============================================================================
# Shared fixture: minimal daily hazard tibble (generate_daily_hazard_impact schema)
# =============================================================================

# Five synthetic years, one location ("Saba").
# Year 1: IRMA-like HUR event (82 kt, 3 days)
# Year 2: moderate HUR (68 kt, 2 days)
# Year 3: TS only (42 kt, 2 days)
# Year 4: calm — no TC activity
# Year 5: severe HUR + follow-up TS (90 kt + 38 kt, two distinct events)

.stress_daily_fixture <- function() {
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
      list(sid = "2017242N16333", en = 1L, days = 250:252, wind = c(72, 82, 68))
    )),
    make_year(2L, list(
      list(sid = "2005280N15305", en = 1L, days = 280:281, wind = c(68, 65))
    )),
    make_year(3L, list(
      list(sid = "2005280N15305", en = 1L, days = 200:201, wind = c(42, 38))
    )),
    make_year(4L, list()),    # calm year
    make_year(5L, list(
      list(sid = "2017242N16333", en = 1L, days = 245:247, wind = c(80, 90, 78)),
      list(sid = "2005280N15305", en = 2L, days = 265:266, wind = c(38, 35))
    ))
  )
}

# Convenience: pre-computed metrics and scored metrics
.stress_metrics <- function() compute_stress_year_metrics(.stress_daily_fixture())
.stress_scored  <- function() aggregate_stress_metrics(.stress_metrics())

# =============================================================================
# compute_stress_year_metrics
# =============================================================================

test_that("compute_stress_year_metrics returns tibble with all required columns", {
  m <- .stress_metrics()

  expect_s3_class(m, "data.frame")
  expected <- c("location", "sim_year", "peak_wind_kt", "n_ts_days",
                "n_hur_days", "n_events", "max_event_dur_days",
                "cum_damage", "max_damage_rate")
  expect_true(all(expected %in% names(m)), info = paste(
    "missing:", paste(setdiff(expected, names(m)), collapse = ", ")
  ))
})

test_that("compute_stress_year_metrics calm year produces all-zero wind/event metrics", {
  m    <- .stress_metrics()
  calm <- m[m$sim_year == 4L, ]

  expect_equal(nrow(calm), 1L)
  expect_equal(calm$peak_wind_kt,    0)
  expect_equal(calm$n_ts_days,       0L)
  expect_equal(calm$n_hur_days,      0L)
  expect_equal(calm$n_events,        0L)
  expect_equal(calm$max_event_dur_days, 0L)
})

test_that("compute_stress_year_metrics peak_wind_kt matches maximum in year", {
  m  <- .stress_metrics()
  yr1 <- m[m$sim_year == 1L, ]
  yr5 <- m[m$sim_year == 5L, ]

  expect_equal(yr1$peak_wind_kt, 82)
  expect_equal(yr5$peak_wind_kt, 90)
})

test_that("compute_stress_year_metrics n_events counts distinct storm events", {
  m   <- .stress_metrics()
  yr5 <- m[m$sim_year == 5L, ]

  expect_equal(yr5$n_events, 2L)   # HUR + follow-up TS
})

test_that("compute_stress_year_metrics integer sim_years filter restricts output", {
  m_all      <- compute_stress_year_metrics(.stress_daily_fixture())
  m_filtered <- compute_stress_year_metrics(.stress_daily_fixture(), sim_years = c(1L, 3L))

  expect_setequal(m_filtered$sim_year, c(1L, 3L))
  expect_equal(nrow(m_filtered), 2L)
})

test_that("compute_stress_year_metrics stops on missing required column", {
  bad <- dplyr::select(.stress_daily_fixture(), -wind_kt)
  expect_error(compute_stress_year_metrics(bad), regexp = "wind_kt")
})

test_that("compute_stress_year_metrics returns empty tibble with correct schema when filter removes all rows", {
  result <- suppressWarnings(
    compute_stress_year_metrics(.stress_daily_fixture(), sim_years = 999L)
  )

  expect_equal(nrow(result), 0L)
  expect_true("sim_year" %in% names(result))
  expect_true("peak_wind_kt" %in% names(result))
})

# =============================================================================
# aggregate_stress_metrics
# =============================================================================

test_that("aggregate_stress_metrics adds composite_score and composite_rank columns", {
  s <- .stress_scored()

  expect_true("composite_score" %in% names(s))
  expect_true("composite_rank"  %in% names(s))
})

test_that("aggregate_stress_metrics composite_score is in [0, 1]", {
  s <- .stress_scored()

  expect_gte(min(s$composite_score), 0 - 1e-9)
  expect_lte(max(s$composite_score), 1 + 1e-9)
})

test_that("aggregate_stress_metrics composite_rank 1 corresponds to highest composite_score", {
  s      <- .stress_scored()
  top_yr <- s[s$composite_rank == 1L, ]

  expect_equal(nrow(top_yr), 1L)
  expect_equal(top_yr$composite_score, max(s$composite_score))
})

test_that("aggregate_stress_metrics constant metric is normalised to 0.5", {
  m          <- .stress_metrics()
  m$n_ts_days <- 1L   # constant across all years — range = 0

  s <- aggregate_stress_metrics(m, metrics_used = "n_ts_days")
  expect_true(all(abs(s$composite_score - 0.5) < 1e-9))
})

test_that("aggregate_stress_metrics named weights are applied in the correct order", {
  m <- .stress_metrics()

  # Down-weight everything except peak_wind_kt
  w  <- c(peak_wind_kt = 1, n_hur_days = 0, n_events = 0,
          max_event_dur_days = 0, cum_damage = 0, max_damage_rate = 0)
  s  <- aggregate_stress_metrics(m, weights = w)

  # Year 5 has highest peak_wind_kt (90 kt) — should rank first
  expect_equal(s$composite_rank[s$sim_year == 5L], 1L)
})

test_that("aggregate_stress_metrics portfolio mode removes location column", {
  m <- .stress_metrics()
  m_two <- dplyr::bind_rows(
    m,
    dplyr::mutate(m, location = "Statia")
  )

  s <- aggregate_stress_metrics(m_two, location_agg = "mean")
  expect_false("location" %in% names(s))
  expect_equal(nrow(s), 5L)          # 5 years, aggregated across 2 locations
})

test_that("aggregate_stress_metrics stops on negative weights", {
  m <- .stress_metrics()
  expect_error(
    aggregate_stress_metrics(m, weights = c(peak_wind_kt = -1, n_hur_days = 1)),
    regexp = "negative|non-negative"
  )
})

# =============================================================================
# select_stress_years
# =============================================================================

test_that("select_stress_years returns exactly k rows", {
  s <- .stress_scored()

  for (k in c(1L, 3L, 5L)) {
    sel <- select_stress_years(s, k = k, method = "top")
    expect_equal(nrow(sel), k, info = paste("k =", k))
  }
})

test_that("select_stress_years method=top returns k highest composite_score rows", {
  s   <- .stress_scored()
  sel <- select_stress_years(s, k = 2L, method = "top")

  top2_expected <- s[order(-s$composite_score), ][1:2, ]$sim_year
  expect_setequal(sel$sim_year, top2_expected)
})

test_that("select_stress_years method=stratified spans the full composite_score range", {
  s   <- .stress_scored()
  sel <- select_stress_years(s, k = 3L, method = "stratified")

  # Lowest score tier should be represented (selection_rank exists)
  expect_equal(nrow(sel), 3L)
  expect_true(min(sel$composite_score) < max(sel$composite_score))
})

test_that("select_stress_years method=diverse selects spread-out years", {
  s   <- .stress_scored()
  sel <- select_stress_years(s, k = 3L, method = "diverse")

  # All selected years should be distinct
  expect_equal(nrow(sel), 3L)
  expect_equal(length(unique(sel$sim_year)), 3L)
})

test_that("select_stress_years selection_rank is an integer sequence 1:k", {
  s   <- .stress_scored()
  sel <- select_stress_years(s, k = 3L)

  expect_setequal(sel$selection_rank, 1L:3L)
})

test_that("select_stress_years returns all rows when nrow <= k", {
  s   <- .stress_scored()
  sel <- select_stress_years(s, k = 10L)    # more than 5 rows available

  expect_equal(nrow(sel), nrow(s))
})
