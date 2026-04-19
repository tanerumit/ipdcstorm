# Tests for spatially-coherent daily hazard generation:
#   .spatial_gen() — shared SID pool across locations
#
# This function requires a full run_hazard_model() output. The fixture is
# cached inside a local environment to avoid re-running the model per test.
#
# Run: Rscript -e "testthat::test_file('tests/testthat/test-hazard-downscale-spatial.R')"

# =============================================================================
# Cached fixture: two-location hazard model run (demo IBTrACS, 30 sim years)
# =============================================================================

.spatial_out <- local({
  .cache <- NULL
  function() {
    if (is.null(.cache)) {
      data_path <- system.file("extdata", "ibtracs_demo.csv", package = "ipdcstorm")
      targets   <- tibble::tribble(
        ~name,    ~lat,     ~lon,
        "Saba",   17.635,  -63.230,
        "Statia", 17.489,  -62.974
      )
      cfg    <- make_hazard_cfg(
        data_path        = data_path,
        simulation_years = 30L,
        climate          = make_climate_cfg(scenario = "stationary")
      )
      .cache <<- suppressWarnings(
        run_hazard_model(cfg, targets = targets, seed = 42L, verbose = FALSE)
      )
    }
    .cache
  }
})

# Helper: suppress warnings from internal event-sampling/damage steps
# (e.g. empty Poisson draws, pressure NA coercions in very short sim runs).
.spatial_gen <- function(...) suppressWarnings(suppressMessages(generate_daily_hazard_impact_spatial(...)))

# =============================================================================
# generate_daily_hazard_impact_spatial
# =============================================================================

test_that("generate_daily_hazard_impact_spatial returns named list with one element per location", {
  out <- .spatial_out()
  res <- .spatial_gen(
    out, location = c("Saba", "Statia"),
    sim_years = 1L:5L, seed = 42L
  )

  expect_type(res, "list")
  expect_setequal(names(res), c("Saba", "Statia"))
})

test_that("generate_daily_hazard_impact_spatial each element has required columns", {
  out <- .spatial_out()
  res <- .spatial_gen(
    out, location = c("Saba", "Statia"),
    sim_years = 1L:3L, seed = 42L
  )

  required <- c(
    "sim_year", "doy", "wind_kt", "surge_m", "event_id",
    "pressure_hpa", "pressure_deficit_hpa", "rmw_km",
    "damage_intensity", "damage_rate", "cum_damage"
  )
  for (loc in names(res)) {
    for (col in required) {
      expect_true(col %in% names(res[[loc]]),
        info = paste("location:", loc, "| missing column:", col))
    }
    expect_false(any(c("location", "scenario", "wind_gust_kt", "event_class") %in% names(res[[loc]])))
  }
})

test_that("generate_daily_hazard_impact_spatial is reproducible with the same seed", {
  out  <- .spatial_out()
  res1 <- .spatial_gen(
    out, location = c("Saba", "Statia"),
    sim_years = 1L:5L, seed = 42L
  )
  res2 <- .spatial_gen(
    out, location = c("Saba", "Statia"),
    sim_years = 1L:5L, seed = 42L
  )

  expect_equal(res1$Saba$wind_kt,   res2$Saba$wind_kt)
  expect_equal(res1$Statia$wind_kt, res2$Statia$wind_kt)
})

test_that("generate_daily_hazard_impact_spatial covers all requested sim_years", {
  out <- .spatial_out()
  res <- .spatial_gen(
    out, location = "Saba",
    sim_years = c(1L, 5L, 10L), seed = 42L
  )

  expect_setequal(unique(res$Saba$sim_year), c(1L, 5L, 10L))
})

test_that("generate_daily_hazard_impact_spatial stores gust_factor as metadata", {
  out    <- .spatial_out()
  gf     <- 1.25
  res    <- .spatial_gen(
    out, location = "Saba",
    sim_years = 1L:5L, gust_factor = gf, seed = 42L
  )

  d <- res$Saba
  expect_identical(attr(d, "gust_factor"), gf)
  expect_false("wind_gust_kt" %in% names(d))
})

test_that("generate_daily_hazard_impact_spatial cum_damage is non-decreasing within each year", {
  out <- .spatial_out()
  res <- .spatial_gen(
    out, location = "Saba",
    sim_years = 1L:10L, seed = 42L
  )

  d <- dplyr::arrange(res$Saba, sim_year, doy)
  by_yr <- split(d, d$sim_year)
  all_mono <- vapply(by_yr, function(df) all(diff(df$cum_damage) >= -1e-12), logical(1L))
  expect_true(all(all_mono), info = "cum_damage is not monotone in some years")
})

test_that("generate_daily_hazard_impact_spatial shares storm SIDs across locations", {
  out <- .spatial_out()
  res <- .spatial_gen(
    out, location = c("Saba", "Statia"),
    sim_years = 1L:30L, seed = 42L
  )

  # Extract SID portion (strip the "_y{yr}_e{n}" suffix)
  extract_sids <- function(event_ids) {
    unique(sub("_y[0-9]+_e[0-9]+$", "", stats::na.omit(event_ids)))
  }
  sids_saba   <- extract_sids(res$Saba$event_id)
  sids_statia <- extract_sids(res$Statia$event_id)

  # Spatial coherence: at least some SIDs should appear in both locations
  shared <- intersect(sids_saba, sids_statia)
  expect_gt(length(shared), 0L)
})

test_that("generate_daily_hazard_impact_spatial resolves duplicate library SIDs consistently", {
  out <- .spatial_out()
  dup_sid <- out$events$storm_id[1L]
  dup_row <- out$events[1L, , drop = FALSE]
  dup_row$peak_wind_kt <- dup_row$peak_wind_kt + 25
  dup_row$storm_intensity_kt <- dup_row$storm_intensity_kt + 25
  dup_row$min_pressure_hpa <- dup_row$min_pressure_hpa - 10
  out$events <- dplyr::bind_rows(out$events, dup_row)

  res <- .spatial_gen(
    out, location = "Saba",
    sim_years = 1L:30L, seed = 42L
  )$Saba

  dup_days <- dplyr::filter(res, !is.na(.data$event_id), grepl(paste0("^", dup_sid, "_y"), .data$event_id))
  if (nrow(dup_days) > 0L) {
    expect_true(any(dup_days$pressure_hpa == dup_row$min_pressure_hpa, na.rm = TRUE))
  } else {
    succeed()
  }
})
