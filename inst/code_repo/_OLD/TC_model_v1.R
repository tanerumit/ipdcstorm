
############################################################
# Stochastic Tropical Cyclone Hazard Model – IBTrACS (NA)
#
# OVERALL GOAL
# ------------
# This script implements a simple stochastic hazard model for
# tropical cyclones at user-defined locations using the
# IBTrACS North Atlantic best-track dataset.
#
# For each location (e.g., Miami, St. Martin), the workflow:
#   1. Extracts all storm track points that pass within a
#      specified distance (buffer) of the location.
#   2. Aggregates these points to storm-level “events” with
#      hazard metrics (min distance, max wind, etc.).
#   3. Classifies each event into categorical severity classes
#      (e.g., Tropical Storm, Cat1–2 close, Major TC) based on
#      wind intensity and distance.
#   4. Counts how many events of each severity occur per year
#      over a chosen period (e.g., 1970–present).
#   5. Fits a Poisson model for each severity class to estimate
#      the mean annual rate (lambda, λ) and the probability of
#      at least one such event occurring in any given year.
#
# The result is a compact statistical model of tropical cyclone
# occurrence by severity at each location, suitable for
# risk, planning, and Monte Carlo simulation.
#
#
# SCRIPT STRUCTURE
# ----------------
#
# 1) prepare_ibtracs()
#    -----------------
#    Input:
#      - Path to IBTrACS North Atlantic CSV file
#        (e.g., "ibtracs.NA.list.v04r01.csv").
#
#    Tasks:
#      - Reads the CSV into R.
#      - Drops the first metadata row.
#      - Selects only the columns needed for hazard modeling:
#        * Storm ID and metadata (SID, SEASON, BASIN, SUBBASIN, NAME)
#        * Time (ISO_TIME)
#        * Position (LAT, LON)
#        * Intensity (WMO_WIND, USA_WIND, USA_SSHS)
#        * Wind radii (USA_R34_*, USA_R50_*, USA_R64_*)
#        * Translation metrics (STORM_SPEED, STORM_DIR)
#        * Land interaction metrics (DIST2LAND, LANDFALL)
#      - Converts LAT/LON to numeric and normalizes longitude to
#        the [-180, 180] range (to work with geosphere).
#      - Parses ISO_TIME into a proper POSIXct datetime (UTC) using
#        multiple possible string formats.
#
#    Output:
#      - A cleaned data frame with one row per track point and
#        standardized coordinates and timestamps.
#
#
# 2) dist_to_target()
#    ----------------
#    A small helper function to compute great-circle distance:
#
#      - Inputs: storm latitude/longitude vectors and a single
#        target latitude/longitude (asset or location).
#      - Uses geosphere::distHaversine() to compute distances
#        in meters, converted to kilometers.
#
#    Used later to compute distance between each track point
#    and each target location.
#
#
# 3) extract_exposure()
#    -------------------
#    Input:
#      - ib: cleaned IBTrACS data (from prepare_ibtracs()).
#      - targets: data frame of target locations with columns:
#          * name – location identifier (e.g., "Miami")
#          * lat, lon – coordinates in decimal degrees.
#      - buffer_km: radius around each target (e.g., 200 km).
#
#    Tasks:
#      - For each target location:
#        * Computes distance from every storm track point to the
#          target (dist_km).
#        * Filters to points where dist_km <= buffer_km.
#      - Stores the resulting subset for each target in a list,
#        indexed by target name.
#
#    Output:
#      - A named list where each element is a data frame with all
#        track points passing within buffer_km km of that location.
#
#
# 4) hazard_summary()
#    -----------------
#    Input:
#      - exposure_list: output from extract_exposure().
#      - buffer_km: same buffer distance used for extraction.
#
#    Tasks (for each location):
#      - Groups track points by storm (SID, SEASON, NAME).
#      - Computes storm-level hazard metrics at the location:
#        * min_dist_km:
#            minimum distance of storm center to the location.
#        * max_wind_nearby:
#            maximum WMO_WIND observed while within the buffer.
#        * cat_max:
#            maximum USA_SSHS category observed near the location.
#        * hours_inside_buffer:
#            total time steps within the buffer, multiplied by 6 h
#            (assuming 6-hour track intervals).
#        * R34_max:
#            maximum 34-kt wind radius across all quadrants
#            (USA_R34_NE, SE, SW, NW).
#        * inside_R34:
#            TRUE/FALSE if the location ever lies within the 34-kt
#            wind footprint, using R34_max as a simple proxy.
#        * storm_speed_mean:
#            mean storm translation speed during the close approach.
#
#    Output:
#      - A list of data frames (one per location). Each row is now
#        a storm-level “event” at that location with summarized
#        hazard metrics.
#
#
# 5) categorize_events()
#    --------------------
#    Input:
#      - haz_df: storm-level hazard summary for one location
#                (one element of the hazard_summary list).
#      - buffer_km: buffer radius (for filtering out non-events).
#      - period_start: first year to include (e.g., 1970) to
#                      restrict to reliable data.
#
#    Tasks:
#      - Adds a numeric year column from SEASON.
#      - Defines a categorical “severity” label for each storm
#        based on simple rules combining intensity and distance.
#        Example default:
#          * "none":
#              weak or distant storms (max_wind_nearby < 34 kt or
#              min_dist_km > buffer_km).
#          * "TS":
#              34–64 kt within the buffer.
#          * "Cat1_2_close":
#              Category 1–2 (cat_max in 1:2) within 100 km.
#          * "Major_TC":
#              Category ≥3 within 150 km.
#          * "other":
#              everything else.
#      - Filters out storms with missing year and restricts to
#        years >= period_start.
#
#    Output:
#      - A data frame with one row per storm per location and
#        columns: year, hazard metrics, and severity class.
#
#    Note:
#      The classification logic is fully configurable. You can
#      change thresholds or category definitions to match
#      your application (e.g., operational thresholds, footprint).
#
#
# 6) fit_hazard_model()
#    -------------------
#    Input:
#      - haz_loc: categorized storm events for one location
#                 (output of categorize_events()).
#
#    Tasks:
#      - Counts how many events of each severity occur per year:
#          n_events(year, severity).
#      - Fills missing years with zero events (so every year has
#        a complete count).
#      - For each severity class:
#        * Estimates the mean annual rate:
#              lambda = mean(n_events).
#        * Computes:
#              p_at_least_one = 1 - exp(-lambda)
#              p_zero         = exp(-lambda)
#        which are the probabilities of at least one event or zero
#        events per year under a Poisson model.
#
#    Output:
#      - A lambda_table with:
#          * severity
#          * lambda      (mean annual event rate)
#          * n_years     (length of record)
#          * p_at_least_one (P[N >= 1 per year])
#          * p_zero         (P[N = 0 per year])
#
#    Interpretation:
#      This is a Poisson occurrence model per severity category,
#      answering questions like:
#        "What is the probability that at least one Major TC occurs
#         within 150 km of this location in any given year?"
#
#
# 7) Example Usage
#    --------------
#    Typical workflow:
#
#      ib       <- prepare_ibtracs(path_to_csv)
#      exposure <- extract_exposure(ib, targets, buffer_km)
#      haz_list <- hazard_summary(exposure, buffer_km)
#
#      # Choose a location, e.g. "Miami"
#      haz_miami <- categorize_events(haz_list[["Miami"]], buffer_km, period_start = 1970)
#      lambda_miami <- fit_hazard_model(haz_miami)
#
#    The resulting lambda_miami table gives the annual rates and
#    probabilities for each event category at Miami.
#
############################################################


############################################################
# PREPARE IBTRACS
############################################################

# prepare_ibtracs()
#
# PURPOSE
# -------
# Reads and preprocesses the IBTrACS North Atlantic best-track dataset
# into a clean, analysis-ready table of storm track points.
#
# This function standardizes coordinates, parses timestamps, and
# selects only those variables required for location-based hazard
# modeling.
#
# INPUT
# -----
# path : character
#   File path to the IBTrACS North Atlantic CSV file
#   (e.g., "ibtracs.NA.list.v04r01.csv").
#
# PROCESSING STEPS
# ----------------
# - Reads the CSV file using readr::read_csv().
# - Removes the first row, which contains metadata rather than data.
# - Selects a subset of columns relevant for hazard modeling:
#     * Storm identifiers and metadata: SID, SEASON, BASIN, SUBBASIN, NAME
#     * Timing: ISO_TIME
#     * Position: LAT, LON
#     * Intensity: WMO_WIND, USA_WIND, USA_SSHS
#     * Wind radii: USA_R34_*, USA_R50_*, USA_R64_*
#     * Storm motion: STORM_SPEED, STORM_DIR
#     * Land interaction: DIST2LAND, LANDFALL
# - Converts latitude and longitude to numeric.
# - Normalizes longitude values from [0, 360] to [-180, 180],
#   which is required for correct geodesic distance calculations.
# - Parses ISO_TIME into a POSIXct datetime object (UTC), allowing for
#   multiple date-time string formats used across IBTrACS versions.
#
# OUTPUT
# ------
# A data frame where each row represents a single storm track point,
# with cleaned coordinates and timestamps.
#
# NOTES
# -----
# - No spatial or temporal filtering is applied here.
# - This function is intentionally generic and reusable across
#   different locations and buffer distances.


prepare_ibtracs <- function(path) {
  
  ib <- readr::read_csv(path, show_col_types = FALSE)
  
  ib %>%
    slice(-1) %>%   # drop metadata row
    select(
      SID, SEASON, BASIN, SUBBASIN, NAME, ISO_TIME,
      LAT, LON, WMO_WIND, USA_WIND, USA_SSHS,
      starts_with("USA_R34"), starts_with("USA_R50"), starts_with("USA_R64"),
      STORM_SPEED, STORM_DIR, DIST2LAND, LANDFALL
    ) %>%
    filter(!is.na(LAT), !is.na(LON)) %>%
    mutate(
      LAT = as.numeric(LAT),
      LON = as.numeric(LON),
      LON = ifelse(LON > 180, LON - 360, LON),  # convert 0–360 → -180–180
      ISO_TIME = lubridate::parse_date_time(
        ISO_TIME, 
        orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S"),
        tz = "UTC"
      )
    )
}


############################################################
# DISTANCE HELPER
############################################################

# dist_to_target()
#
# PURPOSE
# -------
# Computes the great-circle (haversine) distance between storm
# track points and a single target location.
#
# INPUT
# -----
# lat, lon : numeric vectors
#   Latitude and longitude of storm track points (decimal degrees).
#
# t_lat, t_lon : numeric scalars
#   Latitude and longitude of the target location (decimal degrees).
#
# METHOD
# ------
# Uses geosphere::distHaversine(), which accounts for Earth curvature
# and returns distances in meters. Distances are converted to kilometers.
#
# OUTPUT
# ------
# A numeric vector of distances (km), one per storm track point.
#
# NOTES
# -----
# - Assumes coordinates are in WGS84.
# - Vectorized for efficiency.


dist_to_target <- function(lat, lon, t_lat, t_lon) {
  geosphere::distHaversine(cbind(lon, lat), cbind(t_lon, t_lat)) / 1000
}

############################################################
# EXTRACT EXPOSURE FOR MULTIPLE TARGET LOCATIONS
############################################################

# extract_exposure()
#
# PURPOSE
# -------
# Identifies all storm track points that pass within a specified
# distance of one or more target locations.
#
# INPUT
# -----
# ib : data frame
#   Cleaned IBTrACS track-point data returned by prepare_ibtracs().
#
# targets : data frame
#   Target locations with columns:
#     - name : character identifier for the location
#     - lat  : latitude (decimal degrees)
#     - lon  : longitude (decimal degrees)
#
# buffer_km : numeric
#   Radial distance threshold (km) defining potential exposure.
#
# PROCESSING STEPS
# ----------------
# - Iterates over each target location.
# - Computes the distance between every storm track point and the
#   target using dist_to_target().
# - Retains only track points with distance <= buffer_km.
# - Stores the filtered track points for each location in a list.
#
# OUTPUT
# ------
# A named list of data frames.
# Each list element corresponds to one target location and contains
# all storm track points that passed within buffer_km of that location.
#
# NOTES
# -----
# - This step does not yet define "events"; it only filters track points.
# - A storm may appear multiple times (multiple time steps) within a
#   given location’s buffer.

extract_exposure <- function(ib, targets, buffer_km) {
  
  results <- list()
  
  for (i in 1:nrow(targets)) {
    loc <- targets[i, ]
    
    dat_loc <- ib %>%
      mutate(dist_km = dist_to_target(LAT, LON, loc$lat, loc$lon)) %>%
      filter(dist_km <= buffer_km)
    
    results[[loc$name]] <- dat_loc
  }
  
  results
}

############################################################
# HAZARD SUMMARY
############################################################

# hazard_summary()
#
# PURPOSE
# -------
# Aggregates storm track points near each location into storm-level
# hazard events with summary metrics.
#
# INPUT
# -----
# exposure_list : list
#   Output from extract_exposure(), containing track points per location.
#
# buffer_km : numeric
#   Buffer distance (km) used to define proximity.
#
# PROCESSING STEPS
# ----------------
# For each location:
# - Groups track points by storm (SID, SEASON, NAME).
# - Computes storm-level hazard metrics:
#     * min_dist_km:
#         Closest approach of storm center to the location.
#     * max_wind_nearby:
#         Maximum WMO-reported sustained wind observed within the buffer.
#     * cat_max:
#         Maximum Saffir–Simpson category (USA_SSHS) observed nearby.
#     * hours_inside_buffer:
#         Total time the storm spent within the buffer,
#         assuming 6-hour track intervals.
#     * R34_max:
#         Maximum 34-kt wind radius across all quadrants.
#     * inside_R34:
#         Logical indicator of whether the location ever fell
#         inside the 34-kt wind footprint.
#     * storm_speed_mean:
#         Mean translation speed of the storm during close approach.
#
# OUTPUT
# ------
# A list of data frames (one per location).
# Each row represents a single storm event at that location,
# summarized over all nearby track points.
#
# NOTES
# -----
# - Missing wind radii or intensity data are handled safely using
#   suppressWarnings() and NA propagation.
# - This step converts track-level data into event-level data.

hazard_summary <- function(exposure_list, buffer_km) {
  
  lapply(exposure_list, function(df) {
    df %>%
      group_by(SID, SEASON, NAME) %>%
      summarise(
        min_dist_km = suppressWarnings(min(dist_km, na.rm = TRUE)),
        max_wind_nearby = suppressWarnings(max(WMO_WIND, na.rm = TRUE)),
        cat_max = suppressWarnings(max(USA_SSHS, na.rm = TRUE)),
        
        hours_inside_buffer = sum(!is.na(dist_km) & dist_km <= buffer_km) * 6,
        
        R34_max = suppressWarnings(max(
          c(USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW),
          na.rm = TRUE
        )),
        
        inside_R34 = any(dist_km <= R34_max, na.rm = TRUE),
        
        storm_speed_mean = suppressWarnings(mean(STORM_SPEED, na.rm = TRUE)),
        
        .groups = "drop"
      )
  })
}


############################################################
# CATEGORIZE EVENTS (CONFIGURABLE)
############################################################

# categorize_events()
#
# PURPOSE
# -------
# Assigns each storm event to a discrete hazard severity category
# based on intensity and proximity criteria.
#
# INPUT
# -----
# haz_df : data frame
#   Storm-level hazard summary for a single location
#   (one element of hazard_summary()).
#
# buffer_km : numeric
#   Maximum distance used to define relevance of an event.
#
# period_start : integer, default = 1970
#   First year to include in the analysis, typically chosen to
#   restrict the dataset to the satellite era or a reliable period.
#
# CLASSIFICATION LOGIC
# --------------------
# Default categories:
#   - "none":
#       Weak storms (<34 kt) or storms whose closest approach
#       exceeds buffer_km.
#   - "TS":
#       Tropical-storm-force winds (34–64 kt) within the buffer.
#   - "Cat1_2_close":
#       Category 1–2 storms with closest approach <= 100 km.
#   - "Major_TC":
#       Category ≥3 storms with closest approach <= 150 km.
#   - "other":
#       Events not fitting the above categories.
#
# PROCESSING STEPS
# ----------------
# - Converts SEASON to a numeric calendar year.
# - Applies the categorical rules using dplyr::case_when().
# - Removes events with missing year.
# - Restricts the dataset to years >= period_start.
#
# OUTPUT
# ------
# A data frame of storm events with an added "severity" column
# and a cleaned year field.
#
# NOTES
# -----
# - Severity definitions are intentionally simple and configurable.
# - Thresholds can be adapted for specific assets or operational needs.

categorize_events <- function(haz_df, buffer_km, period_start = 1970) {
  
  haz_df %>%
    mutate(
      year = as.integer(SEASON),
      
      severity = dplyr::case_when(
        max_wind_nearby < 34 | min_dist_km > buffer_km ~ "none",
        max_wind_nearby >= 34 & max_wind_nearby < 64   ~ "TS",
        cat_max %in% 1:2 & min_dist_km <= 100          ~ "Cat1_2_close",
        cat_max >= 3 & min_dist_km <= 150              ~ "Major_TC",
        TRUE                                           ~ "other"
      )
    ) %>%
    filter(!is.na(year), year >= period_start)
}


############################################################
# FIT POISSON MODEL & RETURN LAMBDA
############################################################

# fit_hazard_model()
#
# PURPOSE
# -------
# Estimates a Poisson occurrence model for each hazard severity
# category at a given location.
#
# INPUT
# -----
# haz_loc : data frame
#   Categorized storm events for a single location
#   (output of categorize_events()).
#
# PROCESSING STEPS
# ----------------
# - Counts the number of events per year and severity category.
# - Ensures that all years in the observation period are represented,
#   filling missing years with zero events.
# - For each severity category, estimates:
#     * lambda:
#         Mean annual number of events (Poisson rate).
#     * n_years:
#         Length of the observation record (years).
#     * p_at_least_one:
#         Probability of at least one event in a year:
#           1 - exp(-lambda)
#     * p_zero:
#         Probability of zero events in a year:
#           exp(-lambda)
#
# OUTPUT
# ------
# A data frame (lambda table) with one row per severity category
# containing Poisson rate and annual probabilities.
#
# INTERPRETATION
# --------------
# The resulting parameters define a stationary Poisson process
# for event occurrence and directly answer questions such as:
#
#   "What is the probability that at least one Major TC affects
#    this location in any given year?"
#
# NOTES
# -----
# - Assumes independence between years and between severity classes.
# - Provides the statistical foundation for subsequent Monte Carlo
#   simulation or risk analysis.


fit_hazard_model <- function(haz_loc) {
  
  # annual counts
  haz_loc_counts <- haz_loc %>%
    filter(severity != "none") %>%
    count(year, severity, name = "n_events") %>%
    tidyr::complete(
      year = full_seq(range(year), 1),
      severity,
      fill = list(n_events = 0)
    )
  
  # lambdas
  haz_loc_counts %>%
    group_by(severity) %>%
    summarise(
      lambda = mean(n_events),
      n_years = n(),
      p_at_least_one = 1 - exp(-lambda),
      p_zero = exp(-lambda),
      .groups = "drop"
    )
}
