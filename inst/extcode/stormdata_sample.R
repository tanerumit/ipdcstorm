
# Storm data for Saba
haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  lat1r <- lat1 * pi / 180
  lat2r <- lat2 * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon/2)^2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}

ibtracs_box$dist_km <- haversine_km(
  ibtracs_box$lon_num, ibtracs_box$lat_num,
  saba_lon, saba_lat
)

ibtracs_box <- ibtracs_box |> filter(dist_km <= 800)


ibtracs_all <- read_ibtracs_clean(
  "data/ibtracs/ibtracs.NA.list.v04r01.csv",
  basin = "NA", keep_all = FALSE, verbose = TRUE
)

# Filter to 1970+, within 800 km of Saba
ibtracs_box <- ibtracs_all |>
  filter(as.integer(SEASON) >= 1970) |>
  mutate(dist_km = haversine_km(lon, lat, saba_lon, saba_lat)) |>
  filter(dist_km <= 800)

# Keep full tracks for heading computation
sids_near_saba <- unique(ibtracs_box$SID)
ibtracs_demo <- ibtracs_all |> filter(SID %in% sids_near_saba)

# Save — this is ready for make_hazard_cfg(data_path = ...)
write_csv(ibtracs_demo, "data/ibtracs_saba_demo.csv")
