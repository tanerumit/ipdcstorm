

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(geosphere)   # for spatial distance
library(lubridate)


# Coordinates (decimal degrees)
targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)

# Read-in dataset
ibtracs <- read_csv("./data/ibtracs.NA.list.v04r01.csv")

# Select relevant columns
ibtracs_sub <- ibtracs %>% 
  slice(-1) %>%
  select(SID, SEASON, BASIN, SUBBASIN, NAME, ISO_TIME,
    LAT, LON, WMO_WIND, USA_WIND, USA_SSHS,
    starts_with("USA_R34"), starts_with("USA_R50"), starts_with("USA_R64"),
    STORM_SPEED, STORM_DIR, DIST2LAND, LANDFALL) %>%
  filter(!is.na(LAT)) %>%
  filter(!is.na(LON)) %>%
  mutate(LON = as.numeric(LON), LAT = as.numeric(LAT)) %>%
  mutate(LON = ifelse(LON > 180, LON - 360, LON)) %>%  # convert 0–360 → −180–180
  mutate(ISO_TIME = parse_date_time(ISO_TIME, 
          orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S"), tz = "UTC"))
         

# Function to calculate distance to target location
dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}


# Extract events within the buffer distance
buffer_km <- 200
results <- list()

for (i in 1:nrow(targets)) {
  loc <- targets[i, ]
  
  dat_loc <- ibtracs_sub %>%
    mutate(dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon)) %>%
    filter(dist_km <= buffer_km)
  
  results[[loc$name]] <- dat_loc
}

# Hazard summary per location
hazard_summary <- lapply(results, function(df) {
  
  df %>%
    group_by(SID, SEASON, NAME) %>%
    summarise(
      min_dist_km         = suppressWarnings(min(dist_km, na.rm = TRUE)),
      max_wind_nearby     = suppressWarnings(max(WMO_WIND, na.rm = TRUE)),
      cat_max             = suppressWarnings(max(USA_SSHS, na.rm = TRUE)),
      
      # duration inside buffer
      hours_inside_buffer = sum(!is.na(dist_km) & dist_km <= buffer_km) * 6,
      # assuming 6-hour intervals
      
      # footprint: inside wind radii
      R34_max = suppressWarnings(max(c(USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW), na.rm = TRUE)),
      inside_R34 = any(dist_km <= R34_max, na.rm = TRUE),
      storm_speed_mean = suppressWarnings(mean(STORM_SPEED, na.rm = TRUE)), .groups = "drop") 
})


# “What is the probability that at least one event of type X occurs in a given year?”


# Sketch for one location (e.g. St Martin):
  

buffer_km <- 200

haz_loc <- hazard_summary[["St_Martin"]] %>%
  mutate(
    year = as.integer(SEASON),
    severity = case_when(
      max_wind_nearby < 34 | min_dist_km > buffer_km ~ "none",
      max_wind_nearby >= 34 & max_wind_nearby < 64 ~ "TS",
      cat_max %in% 1:2 & min_dist_km <= 100        ~ "Cat1_2_close",
      cat_max >= 3 & min_dist_km <= 150            ~ "Major_TC",
      TRUE                                         ~ "other"
    )
  ) %>%
  filter(!is.na(year), year >= 1970)   # choose period

# annual counts per severity
haz_loc_counts <- haz_loc %>%
  filter(severity != "none") %>%
  count(year, severity, name = "n_events") %>%
  tidyr::complete(
    year = full_seq(range(year), 1),
    severity,
    fill = list(n_events = 0)
  )

# estimate lambdas
lambda_table <- haz_loc_counts %>%
  group_by(severity) %>%
  summarise(
    lambda = mean(n_events),
    n_years = n(),
    p_at_least_one = 1 - exp(-lambda),
    p_zero         = exp(-lambda),
    .groups = "drop"
  )

lambda_table












# How close did it come?	min_dist_km
# How strong at closest approach?	max_wind_nearby
# What category was it?	cat_max
# How long did it linger?	hours_inside_buffer
# Did damaging winds reach the site?	inside_R34
# How big was it?	R34_max
# Was it fast or slow moving?	storm_speed_mean


#WMO_WIND (knots)
# USA_SSHS (derived from USA_WIND)

# Distance to cyclone track
#USA_R64_NE
#USA_R64_SE
#USA_R64_SW
#USA_R64_NW

# Essential columns
# Basin, Subbasin, USA_SSHS, USA_WIND (wind speed used for hurricane category)
# 

# BASIN - na, nORTH aTLANTIC, MM missing, 
# SUBBASIN - CS: CARRIBBEAN SEA, GULF OF MEXICO (GM)
# USA_SSHS
  # KMA_WIND
# REUNION_WIND
#
##


ibtracs_cs <- ibtracs %>% filter(SUBBASIN == "CS")

# Continuous wind speed at closest approach
# Distance to track
# Whether the asset is inside R34/R50/R64 wind radii
# Storm speed and duration