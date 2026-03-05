# Assumes you already have:
# - out <- run_hazard_model(...)
# - targets tibble with name/lat/lon
# - compute_site_winds_full() available (pkgload::load_all("."))

library(dplyr)

site <- "Miami"
lat <- targets$lat[targets$name == site][1]
lon <- targets$lon[targets$name == site][1]

dat_loc <- out$trackpoints[[site]]
tmp <- compute_site_winds_full(dat_loc, target_lat = lat, target_lon = lon)

# ---- 1) Extract the annual-max driver row for a specific year (Miami 2018) ----
driver_2018 <- tmp %>%
  mutate(year = as.integer(format(iso_time, "%Y"))) %>%
  filter(year == 2018, is.finite(V_site_kt)) %>%
  slice_max(order_by = V_site_kt, n = 1, with_ties = FALSE)

print(driver_2018, width = Inf)

# If you want a cleaner view with the most relevant columns:
driver_2018 %>%
  select(
    SID, iso_time, lat, lon, year,
    V_site_kt, V_site_symmetric_kt, Vmax_kt, storm_speed_kt, heading_deg,
    dist_km, bearing_to_target, quadrant,
    rmw_km, RMW_km, RMW_used_km,
    r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm,
    R34_nm, R34_km, R34_missing, R34_is_climo, R34_eff_km,
    r50_ne_nm, r50_se_nm, r50_sw_nm, r50_nw_nm, R50_nm, R50_km,
    r64_ne_nm, r64_se_nm, r64_sw_nm, r64_nw_nm, R64_nm, R64_km
  ) %>%
  print(width = Inf)

# ---- 2) Any other "suspicious" rows (your choice of rules) ----
# Example rule: R34_eff_km unexpectedly small for a TS/HUR storm
suspicious <- tmp %>%
  mutate(year = as.integer(format(iso_time, "%Y"))) %>%
  filter(
    is.finite(Vmax_kt), Vmax_kt >= 34,
    is.finite(R34_eff_km), R34_eff_km > 0,
    R34_eff_km < 120  # tune threshold
  ) %>%
  arrange(R34_eff_km) %>%
  select(
    SID, iso_time, year, Vmax_kt, V_site_kt, dist_km, quadrant,
    R34_eff_km, R34_km, R34_is_climo,
    r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm
  )

print(suspicious, n = 50, width = Inf)

# Optional: focus on the annual-max driver rows per year, then flag suspicious years
drivers_by_year <- tmp %>%
  mutate(year = as.integer(format(iso_time, "%Y"))) %>%
  filter(is.finite(year), is.finite(V_site_kt)) %>%
  group_by(year) %>%
  slice_max(order_by = V_site_kt, n = 1, with_ties = FALSE) %>%
  ungroup()

drivers_suspicious <- drivers_by_year %>%
  filter(is.finite(Vmax_kt), Vmax_kt >= 34,
         is.finite(R34_eff_km), R34_eff_km > 0,
         R34_eff_km < 120) %>%
  select(
    year, SID, iso_time, Vmax_kt, V_site_kt, dist_km, quadrant,
    R34_eff_km, R34_km, R34_is_climo,
    r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm
  ) %>%
  arrange(year)

print(drivers_suspicious, n = nrow(drivers_suspicious), width = Inf)
