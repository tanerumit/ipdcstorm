
##############################################################
# Tropical Cyclone Hazard Assessment Model
# Using IBTrACS Data with Directional Wind Radii
##############################################################
#
# WORKFLOW OVERVIEW:
# 
# This script implements a stochastic hazard model for tropical cyclones
# at user-defined locations using the IBTrACS North Atlantic best-track 
# dataset. The key improvement over basic approaches is the use of 
# directional wind radii rather than assuming symmetric wind fields.
#
# For each target location (e.g., Miami, St. Martin), the workflow:
#
#   1. Extracts all storm track points (6-hourly observations) that pass 
#      within a specified buffer distance (default 200 km) of the location.
#
#   2. For each observation, calculates the bearing from the storm center 
#      to the target location and determines which meteorological quadrant 
#      (NE/SE/SW/NW) the location occupies relative to the storm.
#
#   3. Applies the appropriate directional wind radius (R34, R50, R64) for 
#      that specific quadrant, rather than assuming the maximum radius 
#      applies in all directions. This accounts for the asymmetric nature 
#      of tropical cyclone wind fields.
#
#   4. Aggregates individual track points to storm-level "events" with 
#      hazard metrics including minimum distance, maximum winds, exposure 
#      duration, and wind threshold exceedances.
#
#   5. Classifies each event into categorical severity classes based on 
#      actual wind exposure: none, Tropical Storm (34-63 kt), Category 1-2 
#      Hurricane (64+ kt), or Major Hurricane (Category 3+).
#
#   6. Counts how many events of each severity occur per year over a chosen 
#      reliable period (default: 1970-present, the satellite era).
#
#   7. Fits a Poisson model for each severity class to estimate the mean 
#      annual occurrence rate (lambda, λ) and the probability of at least 
#      one event occurring in any given year.
#
#   8. Uses Monte Carlo simulation to project the frequency and probability 
#      of events over a synthetic 1000-year period.
#
# The result is a statistically robust model of tropical cyclone occurrence 
# by severity at each location, suitable for risk assessment, planning, and 
# long-term hazard characterization. The directional radii approach provides 
# more accurate exposure estimates than methods that assume symmetric wind 
# fields, reducing false positives and improving rate estimates.
#
##############################################################

##############################################################
# 0. Libraries
##############################################################

# Load required packages for data manipulation, spatial analysis, 
# and time series handling.

library(readr)
library(dplyr)
library(lubridate)
library(geosphere)
library(tidyr)
library(purrr)
library(stringr)

set.seed(123)

################################################################################
#### HELPERS
################################################################################

to_num_quiet <- function(x) {
  x <- str_trim(as.character(x))
  
  # normalize common missing tokens
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
  
  # drop common unit/axis words (extend if your file has others)
  x <- str_replace_all(x, regex("\\b(degrees?_north|degrees?_south|degrees?_east|degrees?_west)\\b", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("\\b(deg|degrees?)\\b", ignore_case = TRUE), "")
  x <- str_trim(x)
  
  # parse numbers, but never warn; unparseable -> NA
  out <- suppressWarnings(parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
  
  # optional sentinels
  out[out %in% c(-999, -99, 999)] <- NA_real_
  
  out
}

max_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

min_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

mean_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}


##############################################################
# 4. Distance function
##############################################################

# Define a helper function to calculate the great-circle distance between 
# a storm center position (lat, lon) and a target location (t_lat, t_lon).
# Uses the Haversine formula implemented by the geosphere package, which 
# accounts for the Earth's spherical geometry. Returns distance in 
# kilometers. This function is applied to every storm track point to 
# determine proximity to each target location.

dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}


##############################################################
# 5. NEW: Bearing and quadrant functions
##############################################################

# Implement the directional wind radii methodology by defining functions 
# to determine the geometric relationship between storm center and target.
# 
# - calculate_bearing(): Computes the bearing (azimuth) from the storm 
#   center to the target location, returned in degrees (0-360°) where 
#   0° is North, 90° is East, 180° is South, and 270° is West.
#
# - get_quadrant(): Maps the bearing to one of four meteorological 
#   quadrants (NE, SE, SW, NW) following standard convention. Each 
#   quadrant spans 90° of azimuth.
#
# - get_directional_radius(): Selects the appropriate wind radius value 
#   for the determined quadrant. Tropical cyclone wind fields are 
#   asymmetric, with different radii in different directions. This 
#   function ensures we use the correct radius for the target's actual 
#   position relative to the storm, rather than assuming the maximum 
#   radius applies everywhere (which would overestimate exposure).

# Calculate bearing from storm center (lat, lon) to target (t_lat, t_lon)
# Returns bearing in degrees (0-360), where 0=North, 90=East, 180=South, 270=West

calculate_bearing <- function(lat, lon, t_lat, t_lon) {
  geosphere::bearing(cbind(lon, lat), cbind(t_lon, t_lat))
}

# Determine which quadrant the target is in relative to storm center
# Based on meteorological convention:
#   NE: 0-90 degrees (North to East)
#   SE: 90-180 degrees (East to South)
#   SW: 180-270 degrees (South to West)
#   NW: 270-360 degrees (West to North)
get_quadrant <- function(bearing) {
  # Normalize bearing to 0-360
  bearing <- (bearing + 360) %% 360
  
  case_when(
    bearing >= 0 & bearing < 90    ~ "NE",
    bearing >= 90 & bearing < 180  ~ "SE",
    bearing >= 180 & bearing < 270 ~ "SW",
    bearing >= 270 & bearing < 360 ~ "NW",
    TRUE ~ NA_character_
  )
}

# Select the appropriate wind radius based on quadrant
# Returns the wind radius (in nautical miles) for the given quadrant and threshold
get_directional_radius <- function(quadrant, r34_ne, r34_se, r34_sw, r34_nw) {
  case_when(
    quadrant == "NE" ~ r34_ne,
    quadrant == "SE" ~ r34_se,
    quadrant == "SW" ~ r34_sw,
    quadrant == "NW" ~ r34_nw,
    TRUE ~ NA_real_
  )
}




##############################################################
# 1. Read IBTrACS
##############################################################

# Import the IBTrACS (International Best Track Archive for Climate 
# Stewardship) dataset for the North Atlantic basin. This dataset 
# contains historical tropical cyclone track data with 6-hourly 
# observations including position, intensity, and wind field structure.
# Version 04r01 provides the most complete and quality-controlled record 
# available for the Atlantic basin.

ibtracs <- read_csv("./data/ibtracs.NA.list.v04r01.csv")


##############################################################
# 2. Target locations
##############################################################

# Define the geographic locations for hazard assessment. Each location 
# requires a name identifier, latitude (decimal degrees), and longitude 
# (decimal degrees, negative for West). These coordinates represent the 
# point of interest for which tropical cyclone exposure will be calculated.
# Locations can be cities, islands, ports, or any other sites where 
# storm hazard characterization is needed.

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

# Extract and prepare the necessary variables from the IBTrACS dataset.
# This includes storm identifiers (SID, NAME), temporal information 
# (SEASON, ISO_TIME), position data (LAT, LON), intensity metrics 
# (WMO_WIND, USA_SSHS), and critically, the quadrant-specific wind radii 
# (USA_R34_*, USA_R50_*, USA_R64_*) for the 34, 50, and 64-knot wind 
# thresholds. Coordinates are standardized to decimal degrees with 
# longitude converted to -180 to +180 range. All wind radii columns are 
# explicitly converted to numeric to prevent type conflicts. Records with 
# missing position data are removed as they cannot be spatially analyzed.

ib_sub <- ibtracs %>%
  select(
    SID, SEASON, BASIN, SUBBASIN, NAME, ISO_TIME,
    LAT, LON,
    WMO_WIND, USA_WIND, USA_SSHS,
    starts_with("USA_R34"), starts_with("USA_R50"), starts_with("USA_R64"),
    STORM_SPEED, STORM_DIR, DIST2LAND, LANDFALL
  ) %>%
  mutate(
    LAT = to_num_quiet(LAT),
    LON = to_num_quiet(LON),
    LON = if_else(!is.na(LON) & LON > 180, LON - 360, LON),
    
    ISO_TIME = parse_date_time(
      ISO_TIME,
      orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S"),
      tz = "UTC"
    ),
    
    across(starts_with("USA_R"), to_num_quiet)
  ) %>%
  filter(!is.na(LAT), !is.na(LON))

##############################################################
# Extract track points within 200 km of each location
#    NOW with directional wind radii
##############################################################

# For each target location, filter storm track observations to those 
# within the specified buffer distance (200 km by default). For each 
# retained observation, calculate the bearing and quadrant relative to 
# the target, then extract the directional wind radii for that specific 
# quadrant at the 34, 50, and 64-knot thresholds. Convert radii from 
# nautical miles (as stored in IBTrACS) to kilometers for consistency.
# Determine whether the target location falls within each wind threshold 
# by comparing its distance to the appropriate directional radius. This 
# approach provides accurate, location-specific wind exposure estimates 
# that respect the asymmetric structure of tropical cyclone wind fields.
# Results are stored in a list with one data frame per location.

buffer_km <- 200
results <- list()

for (i in 1:nrow(targets)) {
  loc <- targets[i, ]
  
  dat_loc <- ib_sub %>%
    mutate(
      # Calculate distance and bearing to target
      dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon),
      bearing_to_target = calculate_bearing(LAT, LON, loc$lat, loc$lon),
      quadrant = get_quadrant(bearing_to_target),
      
      # Get directional wind radii for each threshold (in nautical miles)
      R34_directional = get_directional_radius(quadrant, USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
      R50_directional = get_directional_radius(quadrant, USA_R50_NE, USA_R50_SE, USA_R50_SW, USA_R50_NW),
      R64_directional = get_directional_radius(quadrant, USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW),
      
      # Convert to km (1 nautical mile = 1.852 km)
      R34_directional_km = R34_directional * 1.852,
      R50_directional_km = R50_directional * 1.852,
      R64_directional_km = R64_directional * 1.852,
      
      # Determine if location is inside each wind threshold
      inside_R34 = !is.na(R34_directional_km) & dist_km <= R34_directional_km,
      inside_R50 = !is.na(R50_directional_km) & dist_km <= R50_directional_km,
      inside_R64 = !is.na(R64_directional_km) & dist_km <= R64_directional_km
    ) %>%
    filter(dist_km <= buffer_km)
  
  results[[loc$name]] <- dat_loc
}


##############################################################
# 7. Summarise to storm-level hazard metrics
##############################################################

# Aggregate the 6-hourly track point observations to storm-level (per SID) 
# summary metrics. For each unique storm that passed near a location, 
# compute key hazard indicators including the minimum distance achieved, 
# maximum wind intensity observed, highest Saffir-Simpson category reached, 
# and total duration within the buffer zone. Calculate wind exposure metrics 
# using the directional radii: maximum directional radius encountered, 
# hours of exposure to each wind threshold (34, 50, 64 kt), and boolean 
# flags indicating whether the location experienced each wind intensity.
# Also compute the mean storm translation speed, which can be important for 
# rainfall and surge hazards. This storm-level aggregation transforms 
# multiple track points into a single "event" record per storm, providing 
# the unit of analysis for the statistical hazard model.

hazard_summary <- lapply(results, function(df) {
  df %>%
    group_by(SID, SEASON, NAME) %>%
    summarise(
      min_dist_km         = min_na(dist_km),
      max_wind_nearby     = max_na(WMO_WIND),
      cat_max             = max_na(USA_SSHS),
      
      hours_inside_buffer = sum(!is.na(dist_km) & dist_km <= buffer_km) * 6,
      
      R34_dir_max_km = max_na(R34_directional_km),
      R50_dir_max_km = max_na(R50_directional_km),
      R64_dir_max_km = max_na(R64_directional_km),
      
      hours_inside_R34 = sum(inside_R34 %in% TRUE) * 6,
      hours_inside_R50 = sum(inside_R50 %in% TRUE) * 6,
      hours_inside_R64 = sum(inside_R64 %in% TRUE) * 6,
      
      exposed_to_34kt = any(inside_R34 %in% TRUE),
      exposed_to_50kt = any(inside_R50 %in% TRUE),
      exposed_to_64kt = any(inside_R64 %in% TRUE),
      
      storm_speed_mean = mean_na(STORM_SPEED),
      
      .groups = "drop"
    )
})



##############################################################
# 8. Attach year and categorize severity (IMPROVED categories)
##############################################################

# Extract the year from the storm season and classify each storm into 
# discrete severity categories based on actual wind exposure. The 
# classification system uses the directional wind radius analysis to 
# determine whether the location experienced tropical storm-force winds 
# (34+ kt), hurricane-force winds (64+ kt), and combines this with the 
# storm's maximum intensity category. Categories are: "none" (no significant 
# wind exposure), "TS" (tropical storm winds only, 34-63 kt), "Cat1_2" 
# (hurricane-force winds from Category 1-2 storms), and "Major_TC" 
# (hurricane-force winds from Category 3+ major hurricanes). This 
# classification is more physically meaningful than approaches that mix 
# distance and intensity criteria, as it directly reflects the wind hazard 
# experienced at the location. Records without valid year information are 
# filtered out.

hazard_summary <- lapply(hazard_summary, function(df) {
  df %>%
    mutate(
      year = as.integer(SEASON),
      
      # IMPROVED: Severity based on directional wind exposure
      severity = case_when(
        # No significant exposure
        !exposed_to_34kt ~ "none",
        
        # Tropical Storm winds (34-63 kt) based on directional radii
        exposed_to_34kt & !exposed_to_64kt ~ "TS",
        
        # Hurricane winds (64+ kt) based on directional radii
        exposed_to_64kt & cat_max %in% 1:2 ~ "Cat1_2",
        exposed_to_64kt & cat_max >= 3      ~ "Major_TC",
        
        # Catch remaining cases
        TRUE ~ "other"
      )
    ) %>%
    filter(!is.na(year))
})



##############################################################
# 9. Select a location (example: Miami)
##############################################################


# Choose which target location to analyze and restrict the temporal domain 
# to the reliable satellite era (1970 onwards by default). Data quality 
# and completeness in the IBTrACS dataset improves substantially after 1970 
# when satellite monitoring became routine. Earlier records may have 
# incomplete wind radii data or missed storms, which would bias the 
# statistical rate estimates. The 1970 cutoff balances data quality with 
# sample size for robust Poisson parameter estimation.

haz_loc <- hazard_summary[["Miami"]] %>%
  filter(year >= 1970)                               # restrict to reliable period



##############################################################
# 10. Annual counts by severity
##############################################################

# Transform the storm-level event data into an annual time series of 
# event counts by severity category. Exclude "none" severity events as 
# these represent storms that passed nearby but did not produce significant 
# winds at the location. Count the number of events in each severity class 
# for each year, then use complete() to ensure all year-severity combinations 
# are represented, filling missing combinations with zero counts. This creates 
# a complete annual matrix suitable for statistical analysis, ensuring that 
# years with no events are properly represented as zeros rather than missing 
# data, which is critical for accurate rate estimation.

haz_loc_counts <- haz_loc %>%
  filter(severity != "none") %>%
  count(year, severity, name = "n_events") %>%
  tidyr::complete(
    year = full_seq(range(year), 1),
    severity,
    fill = list(n_events = 0)
  )


##############################################################
# 11. Estimate Poisson λ for each severity class
##############################################################

# Fit a Poisson probability model to the annual event counts for each 
# severity category. The Poisson distribution is appropriate for modeling 
# rare, independent events occurring at a constant average rate over time.
# For each severity class, calculate lambda (λ), the mean annual occurrence 
# rate, which is simply the average number of events per year. From λ, 
# derive key risk metrics: the probability of experiencing at least one 
# event in any given year (1 - e^-λ) and the probability of zero events 
# (e^-λ). These parameters form the core of the statistical hazard model 
# and can be used for probabilistic risk assessment, return period 
# calculations, and Monte Carlo simulation of long-term event frequencies.

lambda_table <- haz_loc_counts %>%
  group_by(severity) %>%
  summarise(
    lambda = mean(n_events),
    n_years = n(),
    p_at_least_one = 1 - exp(-lambda),
    p_zero         = exp(-lambda),
    .groups = "drop"
  )

print("Lambda estimates by severity class (with directional wind radii):")
print(lambda_table)



##############################################################
# 12. Monte Carlo simulation for 1000 years
##############################################################

# Generate a synthetic catalog of tropical cyclone events over a 1000-year 
# simulation period using the estimated Poisson rates. For each year and 
# severity combination, draw a random number of events from a Poisson 
# distribution with the appropriate lambda parameter. This stochastic 
# simulation extends the ~50-year historical record to explore the full 
# range of possible outcomes over longer timescales. The 1000-year time frame 
# allows estimation of rare events (e.g., 100-year or 500-year return 
# periods) and characterization of inter-annual and inter-decade variability.

# The simulation assumes stationary (constant λ over time), which may not 
# hold under climate change scenarios but provides a baseline for current 
# hazard characterization.

n_years_sim <- 1000

sim_counts <- tidyr::crossing(
  sim_year = 1:n_years_sim,
  severity = lambda_table$severity) %>%
  left_join(lambda_table %>% select(severity, lambda), by = "severity") %>%
  mutate(n_events = rpois(n(), lambda))

##############################################################
# 13. Probability that at least one event of ANY severity occurs
##############################################################

# Calculate the overall annual probability of tropical cyclone impact at 
# the location, combining all severity categories. For each simulated year, 
# sum the event counts across all severity classes to determine whether at 
# least one event occurred (regardless of intensity). The proportion of 
# years with one or more events provides an empirical estimate of the annual 
# exceedance probability. This metric answers the practical question: "What 
# is the chance this location will be affected by a tropical cyclone in any 
# given year?" It aggregates the risks across all intensity levels into a 
# single summary statistic useful for planning and communication.

p_any_event <- sim_counts %>%
  group_by(sim_year) %>%
  summarise(total = sum(n_events)) %>%
  summarise(p_at_least_one = mean(total >= 1))

print("Probability of at least one TC event in any given year:")
print(p_any_event)


