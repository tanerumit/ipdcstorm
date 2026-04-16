test_that("generate_daily_year_extended classifies daily event_class from realized pulse peak", {
  sampled_events <- tibble::tibble(
    start_date = as.Date(c("2001-01-10", "2001-02-10", "2001-03-10")),
    dur_days = c(1L, 1L, 1L),
    V_peak = c(80, 50, 20),
    event_id = c("evt_hur", "evt_ts", "evt_td"),
    event_class = c("TS", "HUR", "TS"),
    Pc_min_hPa = c(980, 995, 1005),
    dP_max_hPa = c(30, 18, 8),
    RMW_mean_km = c(20, 35, 50)
  )

  daily <- ipdcstorm:::generate_daily_year_extended(
    year = 2001,
    sampled_events = sampled_events,
    pulse_shape = "cosine"
  )

  events <- daily |>
    dplyr::filter(!is.na(.data$event_id)) |>
    dplyr::group_by(.data$event_id) |>
    dplyr::summarise(
      event_class = dplyr::first(stats::na.omit(.data$event_class)),
      realized_peak_kt = max(.data$wind_kt, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$event_id)

  expect_equal(events$event_id, c("evt_hur", "evt_td", "evt_ts"))
  expect_equal(events$event_class, c("HUR", "TD", "TS"))
  expect_equal(events$realized_peak_kt, c(80, 20, 50))
  expect_false(any(events$event_class == "TS" & events$realized_peak_kt >= 64))
  expect_false(any(events$event_class == "HUR" & events$realized_peak_kt < 64))
})
