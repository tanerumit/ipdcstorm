core_track_fixture <- function() {
  tibble::tibble(
    SID = c("AL012020", "AL012020", "AL022020"),
    iso_time = as.POSIXct(
      c("2020-08-01 00:00:00", "2020-08-01 06:00:00", "2020-09-05 00:00:00"),
      tz = "UTC"
    ),
    lat = c(18, 18.4, 19),
    lon = c(-63.4, -63.0, -62.2),
    dist_km = c(30, 45, NA_real_),
    wind_kt = c(90, 95, 45),
    pres_hpa = c(960, 955, 1002),
    poci_hpa = c(1010, 1010, 1012),
    rmw_km = c(35, 35, NA_real_),
    r34_ne_nm = c(90, 90, 55),
    r34_se_nm = c(90, 90, 55),
    r34_sw_nm = c(90, 90, 55),
    r34_nw_nm = c(90, 90, 55),
    r50_ne_nm = c(60, 60, 30),
    r50_se_nm = c(60, 60, 30),
    r50_sw_nm = c(60, 60, 30),
    r50_nw_nm = c(60, 60, 30),
    r64_ne_nm = c(40, 40, 20),
    r64_se_nm = c(40, 40, 20),
    r64_sw_nm = c(40, 40, 20),
    r64_nw_nm = c(40, 40, 20),
    storm_speed_kt = c(12, 12, 8)
  )
}

test_that("geometry helpers return expected lengths and preserve missing values", {
  d <- dist_to_target(c(18, NA_real_), c(-63, -63.2), 18.2, -63.1)
  b <- calculate_bearing(c(18, NA_real_), c(-63, -63.2), 18.2, -63.1)

  expect_length(d, 2)
  expect_length(b, 2)
  expect_true(is.finite(d[1]))
  expect_true(is.finite(b[1]))
  expect_true(is.na(d[2]))
  expect_true(is.na(b[2]))
})

test_that("core rate and event helpers produce deterministic summaries", {
  track <- core_track_fixture()
  heading <- compute_storm_heading(track)
  winds <- compute_site_winds_full(track, target_lat = 18.25, target_lon = -63.05)
  events <- make_storm_events(
    dplyr::mutate(
      winds,
      storm_id = SID,
      storm_class = classify_severity(c(70, 50, 20))[match(SID, c("AL012020", "AL012020", "AL022020"))]
    )
  )

  counts <- compute_annual_counts(
    tibble::tibble(
      year = c(2020L, 2020L, 2021L),
      storm_class = c("TS", "HUR", "TS"),
      storm_id = c("a", "b", "c")
    )
  )
  lambda_tbl <- compute_lambda_table(counts)
  annual_counts <- get_annual_counts(list(
    events = tibble::tibble(
      location = c("Saba", "Saba", "Saba"),
      year = c(2020L, 2020L, 2021L),
      storm_class = c("TS", "HUR", "TS"),
      storm_id = c("a", "b", "c")
    )
  ))
  k_hat <- estimate_k_hat(counts)

  expect_true("heading_deg" %in% names(heading))
  expect_true(all(is.finite(heading$heading_deg[1:2])))
  expect_true(all(c("V_site_kt", "RMW_used_km", "R34_eff_km") %in% names(winds)))
  expect_true(all(winds$V_site_kt <= winds$Vmax_kt, na.rm = TRUE))

  expect_equal(classify_severity(c(20, 40, 80, NA_real_)), c("TD", "TS", "HUR", "unknown"))
  expect_equal(sort(unique(events$storm_id)), c("AL012020", "AL022020"))
  expect_true(all(c("year", "peak_wind_kt", "storm_intensity_kt") %in% names(events)))

  expect_equal(nrow(counts), 4)
  expect_true(all(c("lambda", "prob_annual", "prob_none") %in% names(lambda_tbl)))
  expect_equal(sort(unique(annual_counts$storm_class)), c("HUR", "TS"))
  expect_true(is.list(k_hat))
  expect_true(is.finite(k_hat$k_hat))
})

test_that("climatology size estimators return bounded positive values", {
  r34 <- estimate_R34_climo(c(20, 50, 100), lat = c(18, 18, 18))
  rmw <- estimate_RMW_knaff(c(30, 90, 160), lat = c(18, 18, 30))

  expect_true(is.na(r34[1]))
  expect_true(all(r34[2:3] > 0))
  expect_true(all(is.finite(rmw)))
  expect_true(all(rmw >= 8))
})
