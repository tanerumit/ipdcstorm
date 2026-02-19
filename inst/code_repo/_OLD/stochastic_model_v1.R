
##############################################################
# 0. Libraries
##############################################################

library(readr)
library(dplyr)
library(lubridate)
library(geosphere)
library(tidyr)
library(purrr)

set.seed(123)


##############################################################
# 1. Read IBTrACS
##############################################################

ibtracs <- read_csv("./data/ibtracs.NA.list.v04r01.csv")


##############################################################
# 2. Target locations
##############################################################

targets <- tibble::tribble(
  ~name,        ~lat,      ~lon,
  "St_Martin",   18.0708,  -63.0501,
  "Saba",        17.6350,  -63.2300,
  "Statia",      17.4890,  -62.9740,
  "Puerto_Rico", 18.2208,  -66.5901,
  "Miami",       25.7617,  -80.1918
)


##############################################################
# 3. Select and clean IBTrACS variables
##############################################################

ib_sub <- ibtracs %>%
  select(
    SID, SEASON, BASIN, SUBBASIN, NAME, ISO_TIME,
    LAT, LON,
    WMO_WIND, USA_WIND, USA_SSHS,
    starts_with("USA_R34"), starts_with("USA_R50"), starts_with("USA_R64"),
    STORM_SPEED, STORM_DIR, DIST2LAND, LANDFALL
  ) %>%
  mutate(
    LAT = suppressWarnings(as.numeric(LAT)),
    LON = suppressWarnings(as.numeric(LON)),
    LON = ifelse(LON > 180, LON - 360, LON),
    ISO_TIME = parse_date_time(
      ISO_TIME,
      orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S"),
      tz = "UTC"
    )
  ) %>%
  filter(!is.na(LAT), !is.na(LON))



##############################################################
# 4. Distance function
##############################################################

dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}


##############################################################
# 5. Extract track points within 200 km of each location
##############################################################

buffer_km <- 200
results <- list()

for (i in 1:nrow(targets)) {
  loc <- targets[i, ]
  
  dat_loc <- ib_sub %>%
    mutate(dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon)) %>%
    filter(dist_km <= buffer_km)
  
  results[[loc$name]] <- dat_loc
}


##############################################################
# 6. Summarise to storm-level hazard metrics
##############################################################

hazard_summary <- lapply(results, function(df) {
  df %>%
    group_by(SID, SEASON, NAME) %>%
    summarise(
      min_dist_km         = suppressWarnings(min(dist_km, na.rm = TRUE)),
      max_wind_nearby     = suppressWarnings(max(WMO_WIND, na.rm = TRUE)),
      cat_max             = suppressWarnings(max(USA_SSHS, na.rm = TRUE)),
      
      hours_inside_buffer = sum(!is.na(dist_km) & dist_km <= buffer_km) * 6,
      
      R34_max = suppressWarnings(max(
        c(USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW), na.rm = TRUE
      )),
      
      inside_R34 = any(dist_km <= R34_max, na.rm = TRUE),
      
      storm_speed_mean = suppressWarnings(mean(STORM_SPEED, na.rm = TRUE))
    ) %>%
    ungroup()
})



##############################################################
# 7. Attach year and categorize severity (example categories)
##############################################################

hazard_summary <- lapply(hazard_summary, function(df) {
  df %>%
    mutate(
      year = as.integer(SEASON),
      severity = case_when(
        max_wind_nearby < 34 | min_dist_km > buffer_km ~ "none",
        max_wind_nearby >= 34 & max_wind_nearby < 64   ~ "TS",
        cat_max %in% 1:2 & min_dist_km <= 100          ~ "Cat1_2_close",
        cat_max >= 3 & min_dist_km <= 150              ~ "Major_TC",
        TRUE                                           ~ "other"
      )
    ) %>%
    filter(!is.na(year))
})



##############################################################
# 8. Select a location (example: Miami)
##############################################################

haz_loc <- hazard_summary[["Miami"]] %>%
  filter(year >= 1970)                               # restrict to reliable period



##############################################################
# 9. Annual counts by severity
##############################################################

haz_loc_counts <- haz_loc %>%
  filter(severity != "none") %>%
  count(year, severity, name = "n_events") %>%
  tidyr::complete(
    year = full_seq(range(year), 1),
    severity,
    fill = list(n_events = 0)
  )



##############################################################
# 10. Estimate Poisson λ for each severity class
##############################################################

lambda_table <- haz_loc_counts %>%
  group_by(severity) %>%
  summarise(
    lambda = mean(n_events),
    n_years = n(),
    p_at_least_one = 1 - exp(-lambda),
    p_zero         = exp(-lambda),
    .groups = "drop"
  )

print(lambda_table)



##############################################################
# 11. Monte Carlo simulation for 1000 years
##############################################################

n_years_sim <- 1000

sim_counts <- tidyr::crossing(
  sim_year = 1:n_years_sim,
  severity = lambda_table$severity) %>%
  left_join(lambda_table %>% select(severity, lambda), by = "severity") %>%
  mutate(n_events = rpois(n(), lambda))

##############################################################
# 12. Probability that at least one event of ANY severity occurs
##############################################################

p_any_event <- sim_counts %>%
  group_by(sim_year) %>%
  summarise(total = sum(n_events)) %>%
  summarise(p_at_least_one = mean(total >= 1))

print(p_any_event)

##############################################################
# DONE
##############################################################
