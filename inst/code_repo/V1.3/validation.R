# Validation check
daily_saba |>
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::group_by(sim_year) |>
  dplyr::summarise(
    any_TS  = any(wind_kt >= 34),
    any_HUR = any(wind_kt >= 64),
    .groups = "drop"
  ) |>
  dplyr::summarise(
    p_any_TS  = mean(any_TS),
    p_any_HUR = mean(any_HUR)
  )

#p_any_TS  In x% of years, the island experiences ≥ 1 TS
#p_any_HUR In y% of years, the island experiences ≥ 1 hurricane